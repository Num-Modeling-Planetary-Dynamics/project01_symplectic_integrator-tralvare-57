# astro 591 symplectic integrator
import numpy as np

### this reads data ###
def read_data():
    # use nasa horizons to look up data initialy, impement automated method later
    x
    y
    z
    vx
    vy
    vz

    # ecentric enomilly from book using the e and semimajor axis

### moves simulation by dt ###
def step(self, dt):
    dth = 0.5 * dt
    self.lindrift(Q, P, M_sun, dt)
    self.kick(dth)
    self.drift(dt)
    self.kick(dt)
    self.lindrift(dth)

    return


### computes the SUN DRIFT ###
def lindrift(Q, P, M_sun, dt):
    Q_dot = P / N_sum
    Q_tmp = Q + Q_dot * dt / 2


### find INTERATION accelerations ###
def interaction(Q, P, m, dt)
    sum.


### KEPLER DRIFT ###
def kepler_drift(dt):
    E = [M]
    e = ecc

    def dandby(M, ecc):

        def f(E):
            return E - e * np.sin(E) - M

        def fP(E):
            return 1 - e * np.cos(E)

        def fP2(E):
            return e * np.sin(E)

        def fP3(E):
            return e * np.cos(E)

        def D1(E):
            return -(f(E)) / (fP(E))

        def D2(E):
            return -(f(E)) / (fP(E) + 0.5 * (D1(E) * fP2(E)))

        def D3(E):
            return -(f(E)) / (fP(E) + 0.5 * (D2(E) * fP2(E)) + (1 / 6) * ((D2(E) ** 2) * fP3(E)))

        while True:
            Enew = E[-1] + D3(E[-1])
            E.append(Enew)
            if np.abs((E[-1] - E[-2]) / E[-2]) < accuracy:
                break

    # this is the danby method
    if ecc < np.finfo(np.float64).tiny:
        E0 = 0.0
    else:
        E0 = np.arccos(-(rmag0 - a) / a * ecc))

        if np.sign(np.vdot(r0, v0)) < 0.0:
            E0 = 2 * np.pi - E0

        # solves keplers equations over a single step
        M0 = E0 - ecc * np.sin(E0)

        M = M0 + n * dt
        E = dandby(M, ecc)
        dE = E - E0

        # this is the f and g functions
        f = a / rmag0 * (np.cos(dE) - 1.0) + 1.0
        g = dt + 1.0 / n * (np.sin(dE) - dE)

        # advance position vector
        rvec = f * rvec0 + g * vvec0
        rmag = np.linalg.norm(rvec)

        fdot = -a ** 2 / (rmag * rmag0) * n * np.sin(dE)
        gdot = a / rmag * (np.cos(dE) - 1.0) + 1.0

        vvec = fdot * rvec0 + gdot * vvec0

    return rvec, vvec