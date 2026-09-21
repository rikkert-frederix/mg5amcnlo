"""Detach retained nominal-weight arrays from full multiweight estimates."""
import numpy as np


def nominal_copy(values):
    """Store a nominal slice without retaining its large multiweight parent."""
    return np.asarray(values)[..., 0].copy()
