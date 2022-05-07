import numpy as np
import matplotlib.pyplot as plt
from ddeint import ddeint
from scipy.special import erf

omega_by_L = -0.34
alpha = 0.17
tau = 15
T0 = 9.5  # keep it above 9.5 for oscillations
T1 = 8.5  # keep it above 8.3 for oscillations
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
    dB = omega_by_L * A_T0 - (B_t / tau)
    dA = alpha * quenching(B_T1) * B_T1 - (A_t / tau)
    return [dB, dA]


time_steps = 1e4
cutoff = int(time_steps / 4)

tt = np.linspace(0, 600, int(time_steps))
yy = ddeint(model, tt, lambda t: np.array([(Bmax + Bmin) / 2, (Bmax + Bmin) / 2]))

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.plot(tt[cutoff:], yy[cutoff:, 0], "b")
ax1.set_ylabel("B")
ax1.set_xlabel("Time")
ax2.plot(tt[cutoff:], yy[cutoff:, 1], "r")
ax2.set_ylabel("A")
ax2.set_xlabel("Time")
fig.suptitle(f"{N_D=} {T0+T1=} {tau=}")
plt.show()
