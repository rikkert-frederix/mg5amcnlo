#!/usr/bin/env python3
"""Install frozen alternative PDF archives and match future variation widths."""
import json
import math
from pathlib import Path
import tarfile

import lhapdf
import yaml

from campaign import ROOT, STUDY, WIDTH, MW, MZ, calculator, digest, now, save

SETS = [('CT18NLO',14400,59,1), ('MSHT20nlo_as118',27100,65,3),
        ('NNPDF40_nlo_as_01170',333900,101,1),
        ('NNPDF40_nlo_as_01190',334100,101,1),
        ('NNPDF30_nlo_as_0118',260000,101,2)]


def main():
    benchmark_path = STUDY/'inputs/benchmark.json'
    benchmark = json.loads(benchmark_path.read_text())
    output = STUDY/'inputs/parameter_variations_preflight.json'
    if output.exists():
        raise ValueError('Variation preflight already archived')
    lhapdf.setVerbosity(0)
    records = []
    for name, pdf_id, members, version in SETS:
        archive = STUDY/'downloads'/(name+'.tar.gz')
        directory = ROOT/'LHAPDF/share/LHAPDF'/name
        if not directory.exists():
            with tarfile.open(archive) as package:
                for entry in package.getmembers():
                    path = Path(entry.name)
                    if (path.is_absolute() or '..' in path.parts or path.parts[0] != name
                            or not (entry.isfile() or entry.isdir())):
                        raise ValueError('Unsafe/unexpected PDF archive member: '+entry.name)
                package.extractall(path=directory.parent,filter='data')
        info = directory/(name+'.info')
        metadata = yaml.safe_load(info.read_text())
        assert [metadata[k] for k in ('SetIndex','NumMembers','DataVersion','OrderQCD','NumFlavors')] == [pdf_id,members,version,1,5]
        files = sorted(directory.glob(name+'_*.dat'))
        assert len(files) == members
        pdf = lhapdf.mkPDF(name,0)
        assert math.isfinite(pdf.xfxQ(21,.1,172.5))
        couplings = {str(q):pdf.alphasQ(q) for q in (MW,86.25,MZ,172.5,345.)}
        widths = []
        for base in benchmark['top_widths']:
            alpha = pdf.alphasQ(base['muR_GeV'])
            widths.append(dict(mb_GeV=base['mb_GeV'],w_treatment=base['w_treatment'],
                               scale_factor=base['scale_factor'],muR_GeV=base['muR_GeV'],
                               gamma_lo=base['gamma_lo'],qcd_coefficient=base['qcd_coefficient'],
                               alpha_s_pdf=alpha,
                               gamma_nlo_pdf=base['gamma_lo']+base['qcd_coefficient']*alpha))
        records.append(dict(name=name,id=pdf_id,members=members,DataVersion=version,
                            url='https://lhapdfsets.web.cern.ch/current/'+name+'.tar.gz',
                            archive_sha256=digest(archive),metadata_sha256=digest(info),
                            member_sha256={path.name:digest(path) for path in files},
                            metadata=metadata,couplings=couplings,top_widths=widths,
                            W_width_fixed=benchmark['W_width']['gamma_nlo_pdf']))
        print('Prepared',name,'ID',pdf_id,'members',members,flush=True)
    pdf = lhapdf.mkPDF('NNPDF40_nlo_as_01180',0)
    mass_scan = []
    for mt in (171.5,173.5):
        for treatment in ('onshell','bw'):
            ww = 0. if treatment == 'onshell' else benchmark['W_width']['gamma_nlo_pdf']
            for factor in (.5,1.,2.):
                scale = factor*mt
                row = calculator([WIDTH/'top_decay_width',mt,MW,0.,scale,ww,.118,1.,5])
                values = row['values']
                coefficient = values['Delta_Gamma_QCD^NLO [GeV]']/values['alpha_s(muR)']
                g0 = values['Gamma_LO [GeV]']
                row.update(mt_GeV=mt,mb_GeV=0.,w_treatment=treatment,
                           scale_factor=factor,muR_GeV=scale,gamma_lo=g0,
                           qcd_coefficient=coefficient,alpha_s_pdf=pdf.alphasQ(scale),
                           gamma_nlo_pdf=g0+coefficient*pdf.alphasQ(scale))
                mass_scan.append(row)
    save(output,dict(created_utc=now(),status='inputs prepared; no variation integrations',
                     benchmark_sha256=digest(benchmark_path),PDFs=records,top_mass_scan=mass_scan,
                     convention='The benchmark external W width remains fixed for every PDF/alpha-s/top-mass variation.',
                     reference_note='NNPDF30 widths here use the main weak inputs; literature matching still requires its separate weak inputs and cuts.'))
    print('Twelve top-mass/scale width rows prepared; no cross-section conclusions.')


if __name__ == '__main__':
    main()
