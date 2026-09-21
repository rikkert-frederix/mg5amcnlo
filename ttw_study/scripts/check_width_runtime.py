#!/usr/bin/env python3
"""Physical signed-scale widths through an export's actual compiled objects."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import tempfile

from campaign import ROOT, STUDY, digest, now, save


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('process_dir', type=Path)
    args = parser.parse_args()
    process = args.process_dir.resolve()
    if process.parent != STUDY/'processes':
        raise ValueError('Not a study process')
    metadata = json.loads((process/'study_cards/benchmark_parameters.json').read_text())
    widths = json.loads((STUDY/'inputs/benchmark.json').read_text())['top_widths']
    matched = {row['scale_factor']:row for row in widths
               if row['mb_GeV'] == metadata['decay_mb'] and
               row['w_treatment'] == metadata['top_width']['w_treatment']}
    g0 = matched[1.]['gamma_lo']
    archive = STUDY/'inputs'/(process.name+'_runtime_widths.json')
    if archive.exists():
        raise ValueError('Check already archived')
    base = STUDY/'local/runtime_width_checks'/process.name
    base.mkdir(parents=True,exist_ok=True)
    build = Path(tempfile.mkdtemp(prefix='check_',dir=base))
    subdir = next(iter(sorted(process.glob('SubProcesses/P*'))))
    objects = [subdir/(name+'.o') for name in (
        'fnlo_process_common', 'process_dimensions', 'process_dimensions_bridge',
        'fks_metadata', 'fks_metadata_bridge', 'nlo_contribution_bundle',
        'nlo_decay_metadata', 'decay_chain_metadata', 'decay_chain_parameters')]
    libraries = [process/'lib/libgeneric.a',process/'lib/libpdf.a']
    source = STUDY/'scripts/check_runtime_widths.f90'
    command = ['gfortran','-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow',
               '-I'+str(subdir),'-I'+str(process/'lib'),'-I'+str(process/'Source'),
               str(source),*map(str,objects),'-Wl,--start-group',*map(str,libraries),
               '-Wl,--end-group','-L'+str(ROOT/'LHAPDF/lib'),'-lLHAPDF','-lstdc++',
               '-o',str(build/'check_widths')]
    compiled = subprocess.run(command,cwd=build,text=True,capture_output=True)
    (build/'build.log').write_text(compiled.stdout+compiled.stderr)
    compiled.check_returncode()
    result = subprocess.run(['taskset','-c','0-63',str(build/'check_widths'),str(subdir)],
                            cwd=process,text=True,capture_output=True)
    (build/'run.log').write_text(result.stdout+result.stderr)
    result.check_returncode()
    records = []
    for line in result.stdout.splitlines():
        if not line.startswith('WIDTH_CHECK'):
            continue
        t,a,wt,wa,ct,product = map(float,line.split()[1:])
        rt,ra = matched[t]['gamma_nlo_pdf'],matched[a]['gamma_nlo_pdf']
        expected = [rt,ra,-(rt-g0)/g0-(ra-g0)/g0,g0*g0/(rt*ra)]
        actual = [wt,wa,ct,product]
        assert all(math.isclose(x,y,rel_tol=1.e-12,abs_tol=1.e-13)
                   for x,y in zip(actual,expected)), (actual,expected)
        records.append(dict(top_factor=t,antitop_factor=a,width_top=wt,width_antitop=wa,
                            additive_counterterm=ct,product_denominator_rescaling=product))
    assert len(records) == 9
    save(archive,dict(created_utc=now(),status='passed',process=str(process),
                      compile_command=command,stdout=result.stdout,checks=records,
                      source_sha256=digest(source),
                      object_sha256={str(path):digest(path) for path in objects+libraries},
                      card_sha256=digest(subdir/'decay_card.dat'),
                      limitation='Validates running width and signed-axis normalization, not direct matrix-element scale reweighting.'))
    print('All nine signed top/antitop width and counterterm combinations passed.')


if __name__ == '__main__':
    main()
