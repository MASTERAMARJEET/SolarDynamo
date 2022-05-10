from random import uniform
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
N_D = alpha * omega_by_L * tau**2
epsilon_lim = np.sqrt(3) * 0.01 * alpha


def past(t):
    return np.array([(Bmax + Bmin) / 2, (Bmax + Bmin) / 2])


def quenching(B):
    return (1 + erf(B**2 - Bmin**2)) * (1 - erf(B**2 - Bmax**2)) * 0.25


def model(Y, t, epsilon_lim):
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
    dA = alpha * quenching(B_T1) * B_T1 - ( A_t / tau) + uniform(-epsilon_lim, epsilon_lim)
    return [dB, dA]


time_steps = 1e3
cutoff = 0
# cutoff = int(time_steps / 4)

tt = np.linspace(0, time_steps, int(time_steps + 1))
yy = ddeint(model, tt, past, modelargs=(0,))
yy_noise = ddeint(model, tt, past, modelargs=(epsilon_lim,))

fig, ((ax1,ax3),(ax2, ax4)) = plt.subplots(2, 2, sharex=True)
ax1.plot(tt[cutoff:], yy[cutoff:, 0], "b")
ax1.set_ylabel("B")
ax1.set_xlabel("Time")
ax2.plot(tt[cutoff:], yy[cutoff:, 1], "r")
ax2.set_ylabel("A")
ax2.set_xlabel("Time")
ax3.plot(tt[cutoff:], yy_noise[cutoff:, 0], "b")
ax3.set_ylabel("B")
ax3.set_xlabel("Time")
ax4.plot(tt[cutoff:], yy_noise[cutoff:, 1], "r")
ax4.set_ylabel("A")
ax4.set_xlabel("Time")
fig.suptitle(f"{N_D=} {T0+T1=} {tau=}")
plt.show()
