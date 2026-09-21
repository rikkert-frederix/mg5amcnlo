"""Member-wise PDF reduction, separate from scale envelopes and MC errors.

Call this only AFTER flavour/charge sums and nonlinear transformations.
The last input axis is the ordered member vector, including member zero.
Reference: https://www.lhapdf.org/classLHAPDF_1_1PDFSet.html
The installed, archived LHAPDF release supplies the uncertainty prescription.
"""
import math

import numpy as np


def uncertainty(values, set_name='NNPDF40_nlo_as_01180'):
    import lhapdf

    pdf_set = lhapdf.getPDFSet(set_name)
    values = np.asarray(values,dtype=float)
    if values.ndim < 1 or values.shape[-1] != pdf_set.size:
        raise ValueError('Require the complete ordered PDF-member axis')
    shape = values.shape[:-1]
    valid = np.isfinite(values).all(axis=-1)
    output = {key:np.full(shape,np.nan) for key in
              ('library_central','error_plus','error_minus','error_symmetric')}
    for i, row in enumerate(values.reshape(-1,pdf_set.size)):
        if not valid.flat[i]:
            continue
        result = pdf_set.uncertainty(row.tolist(),100*math.erf(1/math.sqrt(2)),False)
        for key, field in (('library_central','central'),('error_plus','errplus'),
                            ('error_minus','errminus'),('error_symmetric','errsymm')):
            output[key].flat[i] = getattr(result,field)
    output.update(nominal=values[...,0].copy(),valid=valid,
                  library_central_minus_nominal=output['library_central']-values[...,0],
                  set_name=set_name,error_type=pdf_set.errorType,
                  library_version=lhapdf.version(),confidence_level_percent=100*math.erf(1/math.sqrt(2)),
                  convention='Member zero remains nominal. PDF errors are not MC errors or scale-envelope points.')
    return output
