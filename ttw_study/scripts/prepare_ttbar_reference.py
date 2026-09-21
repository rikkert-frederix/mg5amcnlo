#!/usr/bin/env python3
"""Install/check the separate NNPDF3.1 and authors' fixed-scale reference data."""
from pathlib import Path
import tarfile

import lhapdf
import numpy as np
import yaml

from campaign import ROOT, STUDY, digest, now, save


def extract_new(archive, destination, root_name):
    target = destination/root_name
    if target.exists():
        raise ValueError('Refuse to replace an existing extracted archive: '+str(target))
    with tarfile.open(archive) as package:
        for entry in package.getmembers():
            path = Path(entry.name)
            if (path.is_absolute() or '..' in path.parts or path.parts[0] != root_name
                    or not (entry.isfile() or entry.isdir())):
                raise ValueError('Unsafe archive entry: '+entry.name)
        package.extractall(path=destination, filter='data')
    return target


def main():
    output = STUDY/'inputs/ttbar_reference_preflight.json'
    if output.exists():
        raise ValueError('ttbar reference preflight already exists')
    name = 'NNPDF31_nlo_as_0118'
    archive = STUDY/'downloads'/(name+'.tar.gz')
    directory = extract_new(archive, ROOT/'LHAPDF/share/LHAPDF', name)
    info = directory/(name+'.info')
    metadata = yaml.safe_load(info.read_text())
    if [metadata[k] for k in ('SetIndex','NumMembers','DataVersion','OrderQCD','NumFlavors')] != [303400,101,1,1,5]:
        raise ValueError('Unexpected NNPDF3.1 metadata')
    files = sorted(directory.glob(name+'_*.dat'))
    if len(files) != 101:
        raise ValueError('Incomplete PDF set')
    lhapdf.setVerbosity(0)
    for member in range(101):
        pdf = lhapdf.mkPDF(name, member)
        if not np.isfinite(pdf.xfxQ(21,.1,172.5)):
            raise ValueError('Nonfinite PDF interpolation')
    pdf = lhapdf.mkPDF(name,0)
    alpha = pdf.alphasQ(172.5)
    data_archive = STUDY/'downloads/ttbar-expanded-ratios.tar.gz'
    data = extract_new(data_archive, STUDY/'references', 'data-expanded-ratios')
    curves = {}
    for prefix in ('norm-incl-dPhill', 'norm-incl-exp-dPhill'):
        path = data/(prefix+'-NLO-NNPDF-mt.dat')
        values = np.loadtxt(path)
        if values.shape != (10,10) or not np.isfinite(values).all():
            raise ValueError('Unexpected fixed-scale angular data layout')
        if not np.allclose(values[:,0], np.arange(.05,1.,.1), atol=1.e-6, rtol=0.):
            raise ValueError('Unexpected angular bin centres')
        if not np.allclose(values[:,-1], .1, atol=4.e-6, rtol=0.):
            raise ValueError('Unexpected angular bin widths')
        integrals = (values[:,1:8]*values[:,-1,None]).sum(axis=0)
        if not np.allclose(integrals,1.,atol=1.e-4,rtol=0.):
            raise ValueError('Reference angular normalization failed')
        curves[path.name] = dict(sha256=digest(path), integrals=integrals.tolist(),
                                bin_centres=values[:,0].tolist(),
                                published_bin_widths=values[:,-1].tolist())
    save(output,dict(created_utc=now(),status='reference inputs checked; no ttbar integration yet',
                     PDF=dict(name=name,id=303400,DataVersion=1,members=101,
                              url='https://lhapdfsets.web.cern.ch/current/'+name+'.tar.gz',
                              archive_sha256=digest(archive),metadata_sha256=digest(info),
                              member_sha256={p.name:digest(p) for p in files},
                              alpha_s_at_mt=alpha,alpha_s_at_MZ=pdf.alphasQ(91.1876)),
                     physics=dict(mt=172.5,MW=80.385,MZ=91.1876,GF=1.166379e-5,
                                  W_width=2.0928,top_width_LO=1.48063,
                                  top_width_NLO=1.48063-1.18*alpha,
                                  top_width_coefficient=-1.18,
                                  width_convention='quoted rounded coefficient; fixed at mt for scale variations',
                                  production_and_decay_numerator_scale=172.5),
                     reference=dict(url='https://www.precision.hep.phy.cam.ac.uk/wp-content/results/ttbar-decay/data-expanded-ratios.tar.gz',
                                    archive_sha256=digest(data_archive),
                                    README_sha256=digest(data/'README'),curves=curves,
                                    source='arXiv:1901.05407 with authors June 2019 update; also discussed in arXiv:2008.11133 Sec. 3.1.3',
                                    observable='Delta phi(lepton+,lepton-)/pi; fully inclusive in all final-state momenta',
                                    scales=[[1,1],[2,2],[.5,.5],[1,.5],[.5,1],[2,1],[1,2]],
                                    normalized='non-exp files use NLO differential divided by NLO total; exp files expand this ratio',
                                    rounding='one printed width is 0.099997; use exact tenths in analysis, retain published widths in data audit')))
    print(output, flush=True)


if __name__ == '__main__':
    main()
