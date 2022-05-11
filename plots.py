from random import uniform, gauss
import numpy as np
import matplotlib.pyplot as plt
from ddeint import ddeint
from scipy.special import erf

#omega_by_L = -0.34
#alpha = 0.17
tau = 1
T0 =10
T1 = 4
Bmin = 1
Bmax = 7



def past(t):
    return np.array([(Bmax + Bmin) / 2, (Bmax + Bmin) / 2])


def quenching(B):
    return (1 + erf(B**2 - Bmin**2)) * (1 - erf(B**2 - Bmax**2)) * 0.25


def model(Y, t,alpha):
    """
    `omega`: differential rotation
    `L`: length scale of the differential rotation
    `aplha`: coeffecient of alpha-quenching
    `tau`: diffusion time-scale
    """
    B_t, A_t = Y(t)
    _, A_T0 = Y(t - T0)
    B_T1, _ = Y(t - T1)
    dB = -2*alpha* A_T0 - (B_t / tau)
    dA = alpha * quenching(B_T1) * B_T1 - ( A_t / tau)
    return [dB, dA]


time_range = 1e3
time_steps = 1
cutoff = 100

for i in range(70):

    #tau = tau1 + i*0.1
    alpha = -1+0.05*i
    N_D = -2*alpha **2 *tau**2 
    tt = np.linspace(0, time_range, int(time_range / time_steps + 1))
    yy = ddeint(model, tt, past, modelargs=(alpha,))

    fig, (ax1,ax2, ax3) = plt.subplots(3, 1, sharex=True)
    ax1.plot(tt[cutoff:], yy[cutoff:, 1], "r")
    ax1.set_ylabel("A", fontsize=14)
    ax2.plot(tt[cutoff:], yy[cutoff:, 0], "b")
    ax2.set_ylabel("$B_\phi$", fontsize=14)
    ax3.plot(tt[cutoff:], np.square(yy[cutoff:, 0]), "b")
    ax3.set_ylabel("$B_\phi^2$", fontsize=14)
    ax3.set_xlabel("Time", fontsize=14)
    fig.suptitle(f"$N_D$={N_D:0.2f}   $T_0$={T0}   $T_1$={T1}   $\omega$/L={-2*alpha}   $\\tau$={tau}")
    plt.tight_layout()
    #plt.show()
    plt.savefig(str(N_D)+'.png')
    plt.clf
