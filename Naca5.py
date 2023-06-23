import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
import scipy


class Naca5Creator:
    def __init__(self, NACA=39012, n_points=500):
        self.NACA = NACA  # this is name of naca before processing
        self.n_points = n_points
        self.NACA_designator = 5

        self.name = 'NACA ' + str(NACA)
        self.Cl = None
        self.Xf = None
        self.Q = None
        self.thickness = None

        self.m = None
        self.K1 = None

        self.x_camberline = None
        self.y_camberline = None

    def check_name(self):
        self.NACA = str(self.NACA)
        if len(self.NACA) == self.NACA_designator:
            self.Cl = float(self.NACA[0]) * 3 / 20  # design lift coefficient
            self.Xf = float(self.NACA[1]) / 20  # the designator for the position of maximum camber
            self.Q = float(self.NACA[2])  # 0/1 standard or ''reflexed
            self.thickness = float(self.NACA[3:]) / 100  # thickness
        if len(self.NACA) != self.NACA_designator:
            print(f'This NACA only asks for a {self.NACA_designator} number')

    def chev_grid_creator(self):
        xi = np.cos(np.pi/(2*(self.n_points+1))*(2*np.linspace(0,self.n_points+1,self.n_points+1)-1))
        return ((-xi +1)/2)*(1-0) + 0

    def reflex(self):
        pass

    def get_m(self):
        xf = np.array([0.05, 0.1, 0.15, 0.2, 0.25])
        m = np.array([0.0580, 0.1260, 0.2025, 0.2900, 0.3910])

        z = np.polyfit(xf, m, 3)
        f = np.poly1d(z)

        x_new = self.Xf
        y_new = round(f(x_new), 5)
        self.m = y_new
        print(self.m)

    def get_K1(self):
        xf = np.array([0.05,0.1,0.15,0.2,0.25])
        K1 = np.array([361.4,51.64,15.957,6.643,3.230])
        K1 = K1*self.Cl/0.3

        def asymptotic_func(x, a, b, c):
            return a * np.exp(-b * x) + c / x

        # Fit the data points using the asymptotic function
        params, _ = curve_fit(asymptotic_func, xf, K1)


        Xf_new = self.Xf
        self.K1 = asymptotic_func(Xf_new, *params)


    def normal_camberline_function(self,x):
        yc = np.zeros(x.shape)
        yc[x < self.m] = self.K1/6 * ((x[x < self.m])**3 - 3 * self.m * (x[x < self.m])**2 + (self.m)**2 * (3 - self.m) * x[x < self.m])
        yc[x >= self.m] = self.K1/6 * (self.m)**3 * (1-x[x >= self.m])
        return yc

    def get_contour_thickness(self,x):
        # getting the thickness of the contour
        a0 = 0.2969
        a1 = -0.1260
        a2 = -0.3516
        a3 = 0.2843
        a4 = -0.1015
        yt = self.thickness * 5 * (a0 * np.sqrt(x) + a1 * x + a2 * x**2 +a3 * x**3 + a4 * x**4)
        return yt

    def make_airfoil(self):
        self.check_name()
        self.x_camberline = self.chev_grid_creator()
        if self.Q ==0:
            self.get_m()
            self.get_K1()
            self.y_camberline = self.normal_camberline_function(self.chev_grid_creator())


            fig, ax = plt.subplots(figsize=(8, 6))
            ax.axis([-0.25, 1.25, -1, 1])
            ax.scatter(self.x_camberline,self.y_camberline)
            plt.show()

        elif self.Q ==1:
            pass
        else:
            pass




x = Naca5Creator(n_points=30)
x.make_airfoil()