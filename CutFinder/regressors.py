import numpy as np


def bayesian_blocks_gaussian(x, y, sigma=None, penalty=3.0, fitrange=None, **kwargs):
    x = np.asarray(x)
    y = np.asarray(y)
    x = x[y!=-np.inf]
    y = y[y!=-np.inf]

    if fitrange is not None:
        x = x[np.bitwise_and(x >= fitrange[0], x <= fitrange[1])]
        y = y[np.bitwise_and(x >= fitrange[0], x <= fitrange[1])]

    n = len(y)

    if sigma is None:
        sigma = np.ones_like(y)*0.1
    else:
        sigma = np.asarray(sigma)

    order = np.argsort(x)
    x = x[order]
    y = y[order]
    sigma = sigma[order]

    w = 1.0 / sigma**2

    W = np.zeros(n + 1)
    Y = np.zeros(n + 1)
    Z = np.zeros(n + 1)

    W[1:] = np.cumsum(w)
    Y[1:] = np.cumsum(w * y)
    Z[1:] = np.cumsum(w * y * y)

    best = np.zeros(n + 1)
    last = np.zeros(n + 1, dtype=int)

    for j in range(1, n + 1):
        i = np.arange(1, j + 1)

        Wij = W[j] - W[i - 1]
        Yij = Y[j] - Y[i - 1]
        Zij = Z[j] - Z[i - 1]

        chi2 = Zij - (Yij * Yij) / Wij
        fitness = -0.5 * chi2

        scores = best[i - 1] + fitness - penalty
        k = np.argmax(scores)

        best[j] = scores[k]
        last[j] = i[k]

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
    chi2 = np.sum(((y - y_model) / sigma) ** 2)/len(values)
    return edges, np.asarray(values), chi2
