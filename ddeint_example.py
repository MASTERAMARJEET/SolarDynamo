import numpy as np
import matplotlib.pyplot as plt
from ddeint import ddeint


def model(XY, t, d):
    """
    We solve the following system:
    X(t) = 1 (t < 0)
    Y(t) = 2 (t < 0)
    dX/dt = X * (1 - Y(t-d)) / 2
    dY/dt = -Y * (1 - X(t-d)) / 2
    """
    x, y = XY(t)
    xd, yd = XY(t - d)
    return np.array([0.5 * x * (1 - yd), -0.5 * y * (1 - xd)])


past = lambda t: np.array([1, 2])  # 'history' at t<0
tt = np.linspace(2, 30, 20000)
d = 0.2  # set parameter d
yy = ddeint(model, tt, past, fargs=(d,))  # solve the DDE !
plt.plot(yy[:, 0], yy[:, 1], lw=2, label=f"delay = {d:.01}")
plt.show()
