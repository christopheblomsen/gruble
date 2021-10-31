"""
Author: Christophe Blomsen

This code is for
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as c
from mpl_toolkits.mplot3d import Axes3D


class ODE:
    def __init__(self, M, B):
        """
        Constructor, figure it out
        B are how many bodies
        M must be an array of the masses
        """
        self.G = c.G  # Graviational constant
        self.M = M    # Array of masses
        self.B = B    # Number of bodies


    def newton(self, G, M, dist):
        """
        Calculates newton gravitational forces
        """
        num = -G*M*dist
        den = np.linalg.norm(dist)**3
        ans = num/den
        return ans

    def euler_cromer(self, r0, v0, T, dt, t0=0):
        """
        euler cromer intergration
        The arrays are to be understod as the first one is the timestep
        The second row is the one for the body
        The last row is the position part in x, y and z

        Example
        x = np.zeros(N)
        y = np.zeros(N)
        r = np.array(list(zip(x, y))

        Would have given the same result as
        r = np.zeros((N, 2))
        """
        B = self.B
        G = self.G                          # For less writing
        N = int(T/dt)                       # Length of all our vectors
        DIM = len(r0[0])
        #print(f'{N} = {T}/{dt}')
        t = np.zeros(N, float)              # time array
        M = self.M                          # Star mass
        r = np.zeros((N, B, DIM), float)    # Position vector
        v = np.zeros((N, B, DIM), float)    # Velocity vector
        a = np.zeros((N, B, DIM), float)    # Acceleration vector

        t[0] = t0                           # incase you don't start at 0
        r[0, :, :] = r0                     # sets initial conditions
        v[0, :, :] = v0

        for i in range(N-1):                # timestep
            for j in range(len(r0)):        # body number
                for k in range(B):          # forces from the other bodies
                    if j != k:              # to make sure it isn't double counted
                        temp = r[i, j, :] - r[i, k, :]
                        a[i, j, :] += self.newton(G, M[k], temp)
            v[i + 1, :, :] = v[i, :, :] + a[i, :, :]*dt
            r[i + 1, :, :] = r[i, :, :] + v[i + 1, :, :]*dt

        self.r = r
        self.v = v
        self.a = a
        self.t = t
        return r, v, a, t

    def plotyplot(self, bodies_names, figs, lim, r):
        """
        For plotting
        """
        DIM = len(r[0][0])
        if figs == 0:
            fig = plt.figure()   # incase you don't know what size you want
                                 # and go by default
        else:
            fig = plt.figure(figsize=figs)

        if lim != 0:             # change x and y limits from default
            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)

        if DIM > 2:
            ax = fig.add_subplot(111, projection='3d')
            for i in range(self.B):  # goes through all the bodies
                ax.plot(r[:, i, 0], r[:, i, 1], zs=r[:, i, 2], label=f"{bodies_names[i]}")
        else:
            ax = fig.add_subplot(111)
            for i in range(self.B):  # goes through all the bodies
                ax.plot(r[:, i, 0], r[:, i, 1], label=f"{bodies_names[i]}")
                ax.plot(r[0, i, 0], r[0, i, 1], "ro", markersize=3)     # plots start
                ax.plot(r[-1, i, 0], r[-1, i, 1], "ko", markersize=3)   # plots end

        ax.legend(loc=1)  # upper right
        plt.show()

M = np.array([1e12, 1e12])  # in kg
dt = 1e-5
B = 2

r0 = np.array([[-1, 1, 1], [1, -1, 1]])
v0 = np.array([[0, -1, 1], [0, 1, 1]])
T = 3
names = ['A', 'B']

three_body = ODE(M, B)
r, v, a, t = three_body.euler_cromer(r0, v0, T, dt)
figs = (16,16)
three_body.plotyplot(names, figs, 0, r)
