"""
This module implements ddeint, a simple Differential Delay Equation
solver built on top of Scipy's odeint """

# REQUIRES Numpy and Scipy.
import numpy as np
from scipy.integrate import ode
import scipy.interpolate
from tqdm import tqdm


class ddeVar:
    """
    The instances of this class are special function-like
    variables which store their past values in an interpolator and
    can be called for any past time: Y(t), Y(t-d).
    Very convenient for the integration of DDEs.
    """

    def __init__(self, past, tc=0):
        """ past(t) = expression of Y(t) for t<tc """

        self.past = past
        self.tc = tc
        # We must fill the interpolator with 2 points minimum

        self.interpolator = scipy.interpolate.interp1d(
            np.array([tc - 1, tc]),  # X
            np.array([self.past(tc), self.past(tc)]).T,  # Y
            kind="linear",
            bounds_error=False,
            fill_value=self.past(tc),
        )

    def update(self, t, Y):
        """ Add one new (ti,yi) to the interpolator """
        Y2 = Y if (Y.size == 1) else np.array([Y]).T
        self.interpolator = scipy.interpolate.interp1d(
            np.hstack([self.interpolator.x, [t]]),  # X
            np.hstack([self.interpolator.y, Y2]),  # Y
            kind="linear",
            bounds_error=False,
            fill_value=Y,
        )

    def __call__(self, t=0):
        """ Y(t) will return the instance's value at time t """

        return self.past(t) if (t <= self.tc) else self.interpolator(t)


class dde(ode):
    """
    This class overwrites a few functions of ``scipy.integrate.ode``
    to allow for updates of the pseudo-variable Y between each
    integration step.
    """

    def __init__(self, model, jac=None):
        def f2(t, y, args):
            return model(self.Y, t, *args)

        ode.__init__(self, f2, jac)
        self.set_f_params(None)

    def integrate(self, t):

        ode.integrate(self, t)
        self.Y.update(self.t, self.y)
        return self.y

    def set_initial_value(self, Y):

        self.Y = Y  # !!! Y will be modified during integration
        ode.set_initial_value(self, Y(Y.tc), Y.tc)


def ddeint(model, tt, past, modelargs=None):
    """ Solves Delay Differential Equations

    Similar to scipy.integrate.odeint. Solves a Delay differential
    Equation system (DDE) defined by

        Y(t) = past(t) for t<0

        Y'(t) = model(Y,t) for t>= 0

    Where model can involve past values of Y, like Y(t-d).


    Parameters
    -----------

    model
      a function Y,t,args -> Y'(t), where args is optional.
      The variable Y is an instance of class ddeVar, which means that
      it is called like a function: Y(t), Y(t-d), etc. Y(t) returns
      either a number or a numpy array (for multivariate systems).

    tt
      The vector of times [t0, t1, ...] at which the system must
      be solved.

    past
      The 'history function'. A function past(t)=Y(t) for t<0, past(t)
      returns either a number or a numpy array (for multivariate
      systems).

    modelargs
      Additional arguments to be passed to parameter ``model``, if any.


    Examples
    ---------

    We will solve the delayed Lotka-Volterra system defined as

        For t < 0:
        x(t) = 1+t
        y(t) = 2-t

        For t >= 0:
        dx/dt =  0.5* ( 1- y(t-d) )
        dy/dt = -0.5* ( 1- x(t-d) )

    The delay ``d`` is a tunable parameter of the model.

    >>> import numpy as np
    >>> from ddeint import ddeint
    >>> 
    >>> def model(XY,t,d):
    >>>     x, y = XY(t)
    >>>     xd, yd = XY(t-d)
    >>>     return np.array([0.5*x*(1-yd), -0.5*y*(1-xd)])
    >>> 
    >>> past = lambda t : np.array([1+t,2-t]) # 'history' at t<0
    >>> tt = np.linspace(0,30,20000) # times for integration
    >>> d = 0.5 # set parameter d 
    >>> yy = ddeint(model,tt,past,modelargs=(d,)) # solve the DDE !

    """

    dde_ = dde(model)
    dde_.set_initial_value(ddeVar(past, tt[0]))
    dde_.set_f_params(modelargs if modelargs else [])
    results = []
    for dt in tqdm(np.diff(tt)):
        results.append(dde_.integrate(dde_.t + dt))
    # results = [dde_.integrate(dde_.t + dt) for dt in np.diff(tt)]
    return np.array([past(tt[0])] + results)
