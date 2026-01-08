import numpy as np
import numpy as np


def bayesian_blocks_gaussian(x, y, sigma=None, penalty=3.0, fitrange=None):
    x = np.asarray(x)
    y = np.asarray(y)
    x = x[y != -np.inf]
    if sigma is None:
        sigma = np.ones_like(y)
    else:
        sigma = np.asarray(sigma)[y != -np.inf] + 1e-6  # avoid zero division

    y = y[y != -np.inf]

    if fitrange is not None:
        y = y[np.bitwise_and(x >= fitrange[0], x <= fitrange[1])]
        sigma = sigma[np.bitwise_and(x >= fitrange[0], x <= fitrange[1])]
        x = x[np.bitwise_and(x >= fitrange[0], x <= fitrange[1])]
    n = len(y)

    if sigma is None:
        sigma = np.ones_like(y)
    else:
        sigma = np.asarray(sigma)

    # sort by x
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    sigma = sigma[order]

    w = 1.0 / sigma**2

    # cumulative sums
    W = np.zeros(n + 1)
    Y = np.zeros(n + 1)
    Z = np.zeros(n + 1)

    W[1:] = np.cumsum(w)
    Y[1:] = np.cumsum(w * y)
    Z[1:] = np.cumsum(w * y * y)

    best = -np.inf * np.ones(n + 1)
    last = np.zeros(n + 1, dtype=int)
    mean = np.zeros(n + 1)

    best[0] = 0.0
    mean[0] = np.inf  # allows any first block

    for j in range(1, n + 1):
        for i in range(1, j + 1):
            Wij = W[j] - W[i - 1]
            Yij = Y[j] - Y[i - 1]
            Zij = Z[j] - Z[i - 1]

            mu = Yij / Wij
            chi2 = Zij - (Yij * Yij) / Wij
            fitness = -0.5 * chi2

            # monotonic constraint
            if mu > mean[i - 1]:
                continue

            score = best[i - 1] + fitness - penalty

            if score > best[j]:
                best[j] = score
                last[j] = i
                mean[j] = mu

    # backtrack
    change_points = []
    j = n
    while j > 0:
        change_points.append(j)
        j = last[j] - 1
    change_points.append(0)
    change_points = change_points[::-1]

    edges = x[change_points[:-1]]

    values = []
    for i0, i1 in zip(change_points[:-1], change_points[1:]):
        Wij = W[i1] - W[i0]
        Yij = Y[i1] - Y[i0]
        values.append(Yij / Wij)
    # build model values per point and compute chi2
    y_model = np.empty_like(y)
    for i0, i1, v in zip(change_points[:-1], change_points[1:], values):
        y_model[i0:i1] = v
    chi2 = np.sum(((y - y_model) / sigma) ** 2) / len(values)
    return edges, np.asarray(values), chi2
