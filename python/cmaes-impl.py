import numpy as np

from numpy.linalg import norm
from cmaes import CMA


def init_bounds() -> np.ndarray:
    bounds = []
    for i in range(2):
        bounds.append([-np.inf, 0])
    for i in range(8):
        bounds.append([0, np.inf])
    for i in range(3):
        bounds.append([-np.inf, np.inf])

    return np.array(bounds)


def regulator(a1: np.float64, a2: np.float64):
    A = np.zeros((20, 11), np.float64)
    b = np.zeros((20), np.float64)

    Ae = np.array([[1, 1], [-a1, -a2]])
    A1 = np.array([[1, -1, 0], [1, 0, -1]])
    A2 = np.array([[1, 1, 0], [2, 0, 1]])
    A3 = np.array([[1, -1, 0], [2, 0, -1]])
    A4 = np.array([[1, 1, 0], [1, 0, 1]])
    A5 = np.array([[0, 2, 0], [0, 1, 2], [0, 0, 1]])

    be = np.array([1, -(a1 + a2), a1 * a2])

    for i in range(4):
        A[i : i + 2, i : i + 2] = Ae
        A[8 + i * 2 : 11 + i * 2, 8:11] = A5
    A[0:2, 8:11] = A1
    A[2:4, 8:11] = A2
    A[4:6, 8:11] = A3
    A[6:8, 8:11] = A4

    for i in range(4):
        b[8 + i * 3 : 11 + i * 3] = be

    return A, b


def fitness(x: np.array):
    if len(x) != 13:
        raise RuntimeError("invalid x dimension, expected 13")

    A, b = regulator(x[0], x[1])
    y = x[2:]
    val = norm(np.matmul(A, y) - b)

    return val


if __name__ == "__main__":
    optimizer = CMA(mean=np.zeros(13), sigma=1.3, bounds=init_bounds())

    for generation in range(200):
        solutions = []
        for _ in range(optimizer.population_size):
            x = optimizer.ask()
            value = fitness(x)
            solutions.append((x, value))
            # print(f"#{generation} {value} (x={x})")
        print(f"#{generation} {value} (x={x})")
        optimizer.tell(solutions)
