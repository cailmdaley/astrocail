import numpy as np

def chi(data, model, weights=1):
    return np.sum((np.array(data) - np.array(model))**2 * np.array(weights))
