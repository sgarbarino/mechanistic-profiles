import numpy as np

def norm(x):
    return (x-np.min(x))/(np.max(x)-np.min(x))

## ISRA
def isra(A, y, x, maxiter):
    iter = 0
    Aty = np.matmul(A.T,  y)
    while (iter < maxiter):
        x = np.multiply(x, Aty) / np.matmul(A.T, np.matmul(A, x))
        iter = iter + 1
    return x


def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())