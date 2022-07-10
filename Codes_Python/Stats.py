import numpy as np
from scipy.stats import linregress


def stats(X, Y):

    X2 = X[np.logical_not(np.isnan(X)) & np.logical_not(np.isnan(Y))]
    Y2 = Y[np.logical_not(np.isnan(X)) & np.logical_not(np.isnan(Y))]

    # coefficietn directeur et ordonnée à l origine
    coef = linregress(X2, Y2)
    y = coef[0] * X2 + coef[1]

    return coef, y, X2




