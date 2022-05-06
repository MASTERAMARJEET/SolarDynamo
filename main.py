import numpy as np
import matplotlib.pyplot as plt
from ddeint import ddeint
from scipy.special import erf

omega_by_L = -0.34
alpha = 0.17
tau = 15
T0 = 2
T1 = 0.5
Bmin = 1
Bmax = 7
N_D = alpha * omega_by_L * tau ** 2


def quenching(B):
    return (1 + erf(B ** 2 - Bmin ** 2)) * (1 - erf(B ** 2 - Bmax ** 2)) * 0.25


def model(Y, t):
    """
    `omega`: differential rotation
    `L`: length scale of the differential rotation
    `aplha`: coeffecient of alpha-quenching
    `tau`: diffusion time-scale
    """
    B_t, A_t = Y(t)
    _, A_T0 = Y(t - T0)
    B_T1, _ = Y(t - T1)
    dB = omega_by_L * (A_t - A_T0) - B_t / tau
    dA = alpha * quenching((B_t - B_T1)) * (B_t - B_T1) - A_t / tau
    return [dB, dA]


tt = np.linspace(0, 300, int(1e4))
yy = ddeint(model, tt, lambda t: np.array([(Bmax + Bmin) / 2, (Bmax + Bmin) / 2]))
plt.plot(tt, yy[:, 0], label="B")
plt.plot(tt, yy[:, 1], label="A")
plt.text(150, 4, f"{N_D=} {T0+T1=} {tau=}")
plt.legend()
plt.show()
