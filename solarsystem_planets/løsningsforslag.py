# Se koden din fra første semester og AST2000
import numpy as np
import matplotlib.pyplot as plt

class Solsystem:
    def __init__(self, number_of_planets, duration, dt, end = 0, dim = 3):
        n_time_steps_tmp = (duration - end) / dt
        n_time_steps = int(n_time_steps_tmp)
        self.a = np.zeros((n_time_steps, number_of_planets, dim))
        self.v = np.zeros((n_time_steps, number_of_planets, dim))
        self.r = np.zeros((n_time_steps, number_of_planets, dim))
        self.t = np.zeros(n_time_steps)
        self.number_of_planets = number_of_planets
        self.dim = dim
        self.dt = dt
        self.N = n_time_steps
        self.G = 6.67384e-11 # Graviational constant

    def initial_values(self, array_pos, array_vel, masses):
        # Kan bruke en dicitionary for både number_of_planets og arrayene
        # med deres possisjoner for å koble en planet til dens verdier
        number_of_planets = self.number_of_planets
        dim = self.dim
        v = self.v
        r = self.r
        self.masses = masses

        for i in range(number_of_planets):
            v[0, i, :] = array_vel[i]
            r[0, i, :] = array_pos[i]

    def aks_func(self, r_pos):
        M = self.masses
        G = self.G

        r = np.linalg.norm(r_pos)
        a = - G * M / (r**3)
        # import pdb; pdb.set_trace()
        a_vec = a * r_pos
        return a_vec

    def Euler_Cromers(self):
        t = self.t
        a = self.a
        v = self.v
        r = self.r
        N = self.N
        number_of_planets = self.number_of_planets

        for i in range(N - 1):
            t[i + 1] = (i + 1) * dt
            for k in range(number_of_planets):
                a[i, k, :] = self.aks_func(r[i, k, :])
            v[i + 1, :, :] = v[i, :, :] + a[i, :, :] * dt
            r[i + 1, :, :] = r[i, :, :] + v[i + 1, :, :] * dt
        return r, v, a, t

    def plot(self, name_of_planets, save = False, filename = "example_plot.png"):
        x = self.r[:, :, 0]
        y = self.r[:, :, 1]

        for i in range(number_of_planets):
            plt.plot(x[:, i], y[:, i], label = f"{name_of_planets[i]}")
            plt.axis("equal")
            plt.legend()
        if save == True:
            plt.savefig(filename)
        plt.show()

if __name__ == '__main__':
    # Dette er 13. november 2019 midtnatten
    duration = (1.5 * 365.25636) * 24 * 60 * 60 #I sekunder
    dt = 23668 # Seconds
    number_of_planets = 2
    testing = Solsystem(number_of_planets, duration, dt)

    dim = 3
    planet_pos = np.zeros((number_of_planets, dim))
    planet_vel = np.zeros((number_of_planets, dim))
    masses = np.zeros(number_of_planets)
    planet_names = ["Venus", "Earth"]

    # Venus
    planet_pos[0, :] = np.array([3.043847272216039E+07 * 1000,
                                -1.032026891160324E+08 * 1000,
                                -3.215531269702129E+06 * 1000])
    planet_vel[0, :] = np.array([3.332530175747475E+01 * 1000,
                                9.837914285608647E+00 * 1000,
                                -1.788529521785596E+00 * 1000])

    # Earth
    planet_pos[1, :] = np.array([9.460337608621770E+07 * 1000,
                                1.145968643799976E+08 * 1000,
                                -3.391763799369335E+03 * 1000])
    planet_vel[1, :] = np.array([-2.332042112652936E+01 * 1000,
                                1.901369366752253E+01 * 1000,
                                -1.227219686033010E-03 * 1000])
    solar_mass = 1.98855e30
    testing.initial_values(planet_pos, planet_vel, solar_mass)
    testing.Euler_Cromers()
    testing.plot(planet_names, save = True)



# TODO: Kommentarer og enheter på alt, (i dokumentasjonsstilen fra IN1910)
# prøv deg på PEP8 også

# TODO: Lag tester og en vektorisert akselerasjons funksjon
