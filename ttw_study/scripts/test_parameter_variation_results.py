"""Scan-weight and end-to-end conditional-batch fixtures, never predictions."""
import contextlib
import hashlib
import io
import json
from pathlib import Path
import shutil
import tarfile
import tempfile
import unittest

import numpy as np

from campaign import ROOT, STUDY, digest, save
from load_results import POINTS, canonical_columns, load as load_main
from parameter_variation_cards import configure
from parameter_variation_results import columns, vectors, verify_inputs, audit, load, read_batches, load_batches
from replica_statistics import independent_jackknife_many
from read_splits import split_ensembles
import test_parameter_variation_inputs


def labels():
    return ['central value', 'dy']+[
        'dyn=3 muR=%g muF=%g d-6=%g d6=%g' % (r, f, a, t) for r, f, t, a in POINTS]


def synthetic_hwu(factor, error_factor=.01):
    lines = ['##& xmin & xmax & '+' & '.join(labels()), '']
    for config in ('R04_b25', 'R04_b30', 'R04_b40', 'R03_b25', 'R05_b25'):
        for charge in ('W+', 'W-'):
            for local in range(1, 22):
                title = config+' '+charge+(' rates: synthetic' if local == 1 else ' synthetic_%d' % local)
                base = [100., 80., 60., 40., 30., 20., 8., 2., 0.] if local == 1 else [4.]
                lines.append('<histogram> %d "%s"' % (len(base), title))
                for i, value in enumerate(base):
                    value *= factor if charge == 'W+' else 0.
                    varied = [value*(r+2*f+3*t+4*a)/10 for r, f, t, a in POINTS]
                    row = [i+.5, i+1.5, value, abs(value)*error_factor]+varied
                    lines.append(' '.join('%+.14e' % v for v in row))
                lines.extend(['<\\histogram>', ''])
    return ('\n'.join(lines)+'\n').encode()


class TestParameterVariationResults(unittest.TestCase):
    def test_signed_coordinates_ignore_order_and_reject_unsupported_columns(self):
        original = labels()
        reversed_labels = original[:2]+list(reversed(original[2:]))
        selected = columns(reversed_labels)
        self.assertEqual([reversed_labels[i] for i in selected], [original[0]]+original[2:])
        extra = ['delta_mu_'+x+' 3 @aux' for x in ('cen', 'min', 'max')]
        self.assertEqual(columns(original+extra), columns(original))
        for bad in (original[:-1], original+[original[2]], original+['PDF=14400'],
                    original+['unknown @aux'], original[:2]+[x.replace('d-6=', 'd6=') for x in original[2:]]):
            with self.assertRaises(ValueError): columns(bad)
        self.assertEqual(columns(original+['PDF=331700'], [331700]), columns(original))
        with self.assertRaises(ValueError): columns(original+['PDF=331700'], [331700, 331701])

    def test_real_control_headers_have_exactly_the_existing_scale_indices(self):
        queue = json.loads((STUDY/'inputs/rng_pilot_queue_alignment_controls_v1.json').read_text())
        for job in queue['jobs']:
            record = json.loads(Path(job['audit']).read_text())
            with open(record['path']) as stream:
                header = [x.strip() for x in next(stream).split('&')[3:]]
            self.assertEqual(columns(header, range(331700, 331801)), canonical_columns(header, job['variant'])[:82])

    def test_card_audit_checks_written_mass_energy_and_widths(self):
        cases = [('ct18', 'onshell'), ('mt171p5', 'all-bw'), ('as119', 'top-bw'), ('ecm13600', 'onshell')]
        with tempfile.TemporaryDirectory(prefix='scan_written_card_audit_') as folder:
            for index, (scenario, mode) in enumerate(cases):
                process = Path(folder)/str(index)
                test_parameter_variation_inputs.TestParameterVariationInputs.fixture(process, mode)
                with contextlib.redirect_stdout(io.StringIO()):
                    _, archive = configure(process, test_parameter_variation_inputs.INPUTS, scenario, 'Pi', 910001+index, mode)
                manifest = json.loads((archive/'manifest.json').read_text())
                verified = verify_inputs(archive, manifest)
                self.assertEqual(verified['name'], scenario)
                # Even updating the raw checksum cannot conceal a wrong actual
                # input: the physical scenario must still agree numerically.
                path = archive/'param_card.dat'
                from models.check_param_card import ParamCard
                card = ParamCard(str(path))
                card['yukawa'].get((6,)).value += 1.
                card.write(str(path), precision=16)
                manifest['hashes']['param_card.dat'] = digest(path)
                with self.assertRaises(ValueError): verify_inputs(archive, manifest)

    @staticmethod
    def archive_fixture(folder):
        process = folder/'synthetic_export'
        test_parameter_variation_inputs.TestParameterVariationInputs.fixture(process, 'onshell')
        (process/'Source').mkdir()
        for source in ('Source/ranmar.f90', 'SubProcesses/mint_module.f90'):
            shutil.copy2(ROOT/'Template/fNLO'/source, process/source)
        shutil.copy2(ROOT/'Template/fNLO/FixedOrderAnalysis/HwU.f90', process/'FixedOrderAnalysis/HwU.f90')
        seed = 920001
        with contextlib.redirect_stdout(io.StringIO()):
            name, cards = configure(process, test_parameter_variation_inputs.INPUTS, 'as117', 'Pi', seed)
        event = process/'Events'/name
        event.mkdir(parents=True)
        factors = np.linspace(.08, .12, 10)
        final = event/'MADatNLO.HwU'
        final.write_bytes(synthetic_hwu(factors.sum(), .01*np.linalg.norm(factors)/factors.sum()))
        run_log = folder/'synthetic_run.log'
        run_log.write_text('SYNTHETIC TEST\nSetting up grids\nRefining results, step 1\n')
        execution = dict(status='finished', outputs={str(final): digest(final)}, log=str(run_log),
                         data_scope='synthetic fixture only; no physics integration')
        save(cards/'execution.json', execution)
        audit_path = folder/'synthetic_audit.json'
        audit(final, audit_path)
        files, workers, stage_rows = {}, [], []
        for index, factor in enumerate(factors):
            stratum, split = index//5, index % 5+1
            directory = 'SubProcesses/P%d_test/all_G1_%d' % (stratum, split)
            log = ('with seed %d\nRanmar initialization seeds %d %d\n'
                   'Ranmar stream namespace refinement-v1 stage 1 split %d\n'
                   'channel %d fixed proposal\n------- iteration\nTime spent in Total : 1\n'
                   % (seed, 100+index, 200+index, split, stratum))
            raw = {'MADatNLO.HwU': synthetic_hwu(factor), 'log.txt': log.encode(),
                   'input_app.txt': b'synthetic batch\n', 'res.dat': b'synthetic result\n'}
            for filename, content in raw.items(): files[directory+'/'+filename] = content
            job = dict(mint_mode=-1, niters=1, niters_done=1, npoints=100, npoints_done=100,
                       p_dir='P%d_test' % stratum, channel=1, split=split, configs=[1], nchans=1,
                       accuracy=.015, run_mode='all', wgt_frac=1., wgt_mult=.2)
            workers.append(dict(directory=directory, job=job,
                                files={n: hashlib.sha256(v).hexdigest() for n, v in raw.items()}))
            stage_rows.append('<a name=%s></a><PRE>\n%s\n</PRE>' % (directory[len('SubProcesses'):], log))
        training = ('with seed %d\nRanmar initialization seeds 50 60\n'
                    'Ranmar stream namespace refinement-v1 stage 0 split 0\nTime spent in Total : 1\n' % seed)
        history = {source: (process/source).read_bytes() for source in ('Source/ranmar.f90', 'SubProcesses/mint_module.f90')}
        history['Events/'+name+'/alllogs_0.html'] = ('<a name=/P0_test/all_G1></a><PRE>\n'+training+'\n</PRE>').encode()
        history['Events/'+name+'/alllogs_1.html'] = '\n'.join(stage_rows).encode()
        files.update(history)
        archive = folder/'synthetic_workers.tar.gz'
        with tarfile.open(archive, 'w:gz') as output:
            for name, raw in files.items():
                info = tarfile.TarInfo(name)
                info.size = len(raw)
                output.addfile(info, io.BytesIO(raw))
        record = dict(audit=str(audit_path), audit_sha256=digest(audit_path), archive_sha256=digest(archive),
                      seed=seed, workers=workers, worker_count=len(workers),
                      history_files={n: hashlib.sha256(v).hexdigest() for n, v in history.items()})
        record_path = archive.with_suffix('.json')
        save(record_path, record)
        return record_path, factors, audit_path

    def test_complete_worker_chain_preserves_covariance_and_stage_history(self):
        with tempfile.TemporaryDirectory(prefix='synthetic_scan_batches_') as directory:
            folder = Path(directory)
            record, factors, audit_path = self.archive_fixture(folder)
            output = folder/'synthetic_batches.npz'
            read_batches(record, output)
            result = load_batches(output)
            self.assertEqual(result['contributions'].shape, (10, 290, 82))
            np.testing.assert_allclose(result['contributions'].sum(axis=0), result['values'], atol=1.e-12)
            metadata = result['batch_audit']
            self.assertEqual(metadata['refinement_rng_audit']['distinct_initialization_pairs'], 11)
            self.assertEqual(metadata['refinement_rng_audit']['final_workers_checked'], 10)
            populations = split_ensembles(result['contributions'][:, :9], result['strata'])
            estimate = independent_jackknife_many(populations, lambda *means: sum(means))
            expected_variance = sum(np.var(factors[start:start+5]*5, ddof=1)/5 for start in (0, 5))
            self.assertAlmostEqual(estimate['mc_error'][0, 0]**2, 100.**2*expected_variance, places=10)
            cov = estimate['covariance_factor'][:, :, 0]
            self.assertAlmostEqual(cov[:, 0]@cov[:, 3], 100.*40.*expected_variance, places=10)
            fraction = independent_jackknife_many(populations, lambda *means: sum(means)[3]/sum(means)[0])
            np.testing.assert_allclose(fraction['value'], .4, atol=1.e-14)
            self.assertLess(float(fraction['mc_error'].max()), 1.e-14)
            with self.assertRaises(ValueError): read_batches(record, output)
            # Physical input corruption must be caught again when loading the
            # batch archive, rather than trusting its earlier successful audit.
            report = json.loads(audit_path.read_text())
            card = Path(report['path']).parents[2]/'study_cards'/Path(report['path']).parent.name/'run_card.dat'
            card.write_text(card.read_text()+'\n# changed after the audit\n')
            with self.assertRaises(ValueError): load_batches(output)

    def test_changed_stage_stream_and_worker_sum_cannot_pass_with_new_checksums(self):
        for kind in ('reused_training_stream', 'missing_training_history', 'changed_worker_weight', 'unequal_batches'):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory(prefix='synthetic_scan_reject_') as directory:
                folder = Path(directory)
                record_path, factors, _ = self.archive_fixture(folder)
                record = json.loads(record_path.read_text())
                archive_path = record_path.with_suffix('.gz')
                with tarfile.open(archive_path, 'r:gz') as archive:
                    files = {member.name: archive.extractfile(member).read() for member in archive}
                stage0 = next(name for name in record['history_files'] if name.endswith('alllogs_0.html'))
                if kind == 'reused_training_stream':
                    files[stage0] = files[stage0].replace(b'initialization seeds 50 60', b'initialization seeds 100 200')
                    record['history_files'][stage0] = hashlib.sha256(files[stage0]).hexdigest()
                elif kind == 'missing_training_history':
                    del files[stage0]
                    del record['history_files'][stage0]
                elif kind == 'changed_worker_weight':
                    worker = record['workers'][0]
                    filename = worker['directory']+'/MADatNLO.HwU'
                    files[filename] = synthetic_hwu(1.1*factors[0])
                    worker['files']['MADatNLO.HwU'] = hashlib.sha256(files[filename]).hexdigest()
                else:
                    record['workers'][0]['job'].update(npoints=200, npoints_done=200)
                with tarfile.open(archive_path, 'w:gz') as archive:
                    for name, raw in files.items():
                        info = tarfile.TarInfo(name)
                        info.size = len(raw)
                        archive.addfile(info, io.BytesIO(raw))
                record['archive_sha256'] = digest(archive_path)
                save(record_path, record)
                expected = {'reused_training_stream': 'collision', 'missing_training_history': 'per-stage',
                            'changed_worker_weight': 'Worker sum', 'unequal_batches': 'Unequal batch'}[kind]
                output = folder/'rejected_batches.npz'
                with self.assertRaisesRegex(ValueError, expected): read_batches(record_path, output)
                self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main()
