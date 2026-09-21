"""Complete six-prescription observables formed before scale/PDF reduction."""
import numpy as np

from full_flavour_statistics import VARIANTS
from replica_statistics import ratio

LABELS=VARIANTS+('Pi_minus_S','Pi_over_S','S_over_P','PiD_minus_D',
                 'production_related_product_shift')
ABSOLUTE_LABELS=LABELS+('S_minus_P_minus_D_plus_LO',)
CHARGE_LABELS=('plus','minus','sum','plus_over_minus','asymmetry')
CHARGE_SHAPE_LABELS=('plus','minus','combined','plus_shape_over_minus_shape','shape_asymmetry')


def prescription_observables(totals,charge,transform=None):
    """The strict linear identity is valid only before normalization."""
    rows={variant:(totals[charge,variant] if transform is None else transform(totals[charge,variant]))
          for variant in VARIANTS}
    output=[rows[variant] for variant in VARIANTS]
    output.extend([rows['Pi']-rows['S'],ratio(rows['Pi'],rows['S']),ratio(rows['S'],rows['P']),
        rows['PiD']-rows['D'],(rows['Pi']-rows['S'])-(rows['PiD']-rows['D'])])
    if transform is None:
        output.append(rows['S']-rows['P']-rows['D']+rows['LO'])
    return np.asarray(output)


def normalize_spectrum(total):
    """Input ends with its matched fiducial/jet-sector normalization rate."""
    if np.asarray(total).ndim<2 or len(total)<2:
        raise ValueError('Require spectrum bins followed by a parent rate')
    return ratio(total[:-1],total[-1])


def charge_observables(totals,transform=None):
    plus=np.asarray([totals['plus',variant] for variant in VARIANTS])
    minus=np.asarray([totals['minus',variant] for variant in VARIANTS])
    combined=plus+minus
    if transform is not None:
        combined=np.asarray([transform(row) for row in combined])
        plus=np.asarray([transform(row) for row in plus])
        minus=np.asarray([transform(row) for row in minus])
    total=plus+minus
    return np.asarray([plus,minus,combined,ratio(plus,minus),ratio(plus-minus,total)])
