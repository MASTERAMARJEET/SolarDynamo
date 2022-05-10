import os.path
from random import uniform
import numpy as np
import matplotlib.pyplot as plt
from ddeint import ddeint
from scipy.special import erf
# from numpy.random import uniform, standard_normal

omega_by_L = -0.34
# alpha = 0.17
tau = 15
T0 = 2
T1 = 0.5
Bmin = 1
Bmax = 7
# N_D = alpha * omega_by_L * tau**2
# epsilon_lim = 1


def past(t):
    return np.array([(Bmax + Bmin) / 2, (Bmax + Bmin) / 2])


def quenching(B):
    return (1 + erf(B**2 - Bmin**2)) * (1 - erf(B**2 - Bmax**2)) * 0.25


def model(Y, t, alpha, epsilon_lim):
    """
    `omega`: differential rotation
    `L`: length scale of the differential rotation
    `aplha`: coeffecient of alpha-quenching
    `tau`: diffusion time-scale
    """
    # epsilon = gauss(0, epsilon_lim)
    epsilon = uniform(-epsilon_lim, epsilon_lim)
    B_t, A_t = Y(t)
    _, A_T0 = Y(t - T0)
    B_T1, _ = Y(t - T1)
    dB = omega_by_L * A_T0 - (B_t / tau)
    dA = alpha * quenching(B_T1) * B_T1 - ( A_t / tau) + epsilon
    return [dB, dA]


time_range = 1e4
time_steps = 1
cutoff = 0

tt = np.linspace(0, time_range, int(time_range / time_steps + 1))


for alpha in [0.24]:
    fig, ax = plt.subplots(1, 1, sharex=True)
    for percent, color in zip([0.00001, 0.00005, 0.0001], ["b","g", "r", "k"]):
        data_file = f"data/a{alpha}_p{percent}_{time_range}_{time_steps}.txt"
        if os.path.isfile(data_file):
            print("file found.. Nice!")
            yy = np.loadtxt(data_file)
        else:
            epsilon_lim = np.sqrt(3) * percent * alpha
            yy = ddeint(model, tt, past, modelargs=(alpha, epsilon_lim,))
            np.savetxt(data_file, yy)
        B = yy[:,0]
        dB = np.diff(B)

        ax.plot(B[1:], dB, color, label=f"{alpha=} {percent=}")
        ax.set_xlabel("$B_\phi$")
        ax.set_ylabel("$dB_\phi$")
    plt.legend()
    plt.show()
# fig.suptitle(f"{N_D=:0.4f} {T0+T1=} {tau=}")
# plt.show()
