import numpy as np
import matplotlib.pyplot as plt
from ddeint import ddeint


def model(XY, t, d):
    x, y = XY(t)
    xd, yd = XY(t - d)
    return np.array([0.5 * x * (1 - yd), -0.5 * y * (1 - xd)])


past = lambda t: np.array([1 + t, 2 - t])  # 'history' at t<0
tt = np.linspace(0, 30, 20000)  # times for integration
d = 0.5  # set parameter d
yy = ddeint(model, tt, past, fargs=(d,))  # solve the DDE !
plt.plot(tt, yy)
plt.show()
