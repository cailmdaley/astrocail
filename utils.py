import numpy as np

def chi(data, model, weights=1):
    return np.sum((data - model)**2 * weights)
