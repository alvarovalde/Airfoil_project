import math
import os
import numpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from random import choice
from math import sqrt
from pathlib import Path
import glob

matplotlib.use('Qt5Agg')
import pandas as pd
from directory_management import clean_directory
import json


class Naca4Creator:
    def __init__(self, NACA=2309, n_points=500,JSON= True, export_camberline = False,accuracy = 4):
        """
        This class creates the naca 4 series outline as a series of (X,Y) points
         and exports it:
            - As a JSON file if JSON = True
            - As a .dat file (not recomended)
        It also accepts:
         + the NACA denomination number (defaults to NACA2309)
         + the number of points in the final outline (defaults to 500)
         + whether it exports the (X,Y) coordinates of the camberline
         + the significant figures of the

        """
        self.accuracy = accuracy
        self.JSON = JSON

        self.NACA = NACA  # this is name of naca before processing
        self.n_points = math.ceil(n_points / 2)
        self.NACA_designator = 4

        self.name = 'NACA ' + str(NACA)
        self.max_camber = None  # m
        self.loc_mc = None  # p
        self.thickness = None  # t

        # point distributions
        self.y_camberline = None
        self.theta_contour = None
        self.y_thickness = None

        # final distributions
        self.camberline = export_camberline
        self.coordinates = None  # coordinates of the thickness distribution
        self.cor_camb= None

    def check_name(self):
        """
        Checks the given string and if correct retrieves the characteristics
        """
        self.NACA = str(self.NACA)
        try:
            if len(self.NACA) == self.NACA_designator:
                self.max_camber = float(self.NACA[0]) / 100
                self.loc_mc = float(self.NACA[1]) / 10
                self.thickness = float(self.NACA[2:]) / 100
            if len(self.NACA) != self.NACA_designator:
                print(f'This NACA only asks for a {self.NACA_designator} number')
        except SyntaxError as err:
            print('Check the name',err)
            quit()

    def grid_uniform(self):
        return np.round(np.linspace(0, 1, self.n_points),self.accuracy)

    def grid_chevychev(self):
        xi = np.cos(np.pi / (2 * (self.n_points + 1)) * (2 * np.linspace(0, self.n_points + 1, self.n_points) - 1))
        return ((-xi + 1) / 2) * (1 - 0) + 0

    def camberline_function(self, x):
        yc = np.zeros(x.shape)
        yc[x < self.loc_mc] = self.max_camber * x[x < self.loc_mc] / self.loc_mc ** 2 * (
                    2 * self.loc_mc - x[x < self.loc_mc])
        yc[x >= self.loc_mc] = self.max_camber * (1 - x[x >= self.loc_mc]) / (1 - self.loc_mc) ** 2 * (
                    1 + x[x >= self.loc_mc] - 2 * self.loc_mc)
        return yc

        # getting theta function for the plotting of the contour

    def theta_funct(self, x):
        if self.max_camber > 0:
            dydx = np.zeros(x.shape)
            dydx[x < self.loc_mc] = 2 * self.max_camber / self.loc_mc ** 2 * (self.loc_mc - x[x < self.loc_mc])
            dydx[x >= self.loc_mc] = 2 * self.max_camber / (1 - self.loc_mc) ** 2 * (self.loc_mc - x[x >= self.loc_mc])
            theta = np.arctan(dydx)
        else:
            theta = 0
        return theta

    def get_contour_thickness(self, x):
        # getting the thickness of the contour
        a0 = 0.2969
        a1 = -0.1260
        a2 = -0.3516
        a3 = 0.2843
        a4 = -0.1015
        yt = self.thickness * 5 * (a0 * np.sqrt(x) + a1 * x + a2 * x ** 2 + a3 * x ** 3 + a4 * x ** 4)
        return yt

    def create_point_distribution(self):
        x_chev = self.grid_chevychev()

        xu = x_chev - self.y_thickness * np.sin(self.theta)
        xl = x_chev + self.y_thickness * np.sin(self.theta)

        yu = self.y_camberline + self.y_thickness * np.cos(self.theta)
        yl = self.y_camberline - self.y_thickness * np.cos(self.theta)

        cor_up = np.array((xu, yu)).transpose()
        cor_down = np.array((xl, yl)).transpose()

        self.coordinates = np.round(np.concatenate((cor_down, cor_up)),self.accuracy)  # (102,1)
        json_list = self.coordinates.tolist()

        if self.camberline:
            self.cor_camb = np.round(np.array((self.grid_chevychev(), self.y_camberline)).transpose(),self.accuracy).tolist()  # (54,1)


        self.foilJSON = {
            "name": f"{self.name}",
            "camberline": self.camberline,
            "points": json_list,
            "camber points": self.cor_camb

        }

        return self.coordinates

    def send_array(self):
        # clean directory
        clean_directory()
        # create new file
        if self.JSON:
            with open('FoilToAnalize/foil.json', 'w') as file:
                json.dump(self.foilJSON, file, indent=4)
        else:
            np.savetxt(f"FoilToAnalize/{self.name}.dat", self.coordinates, delimiter=' ')
    def make_foil_structure(self):
        self.y_camberline = self.camberline_function(self.grid_chevychev())
        self.theta = self.theta_funct(self.grid_chevychev())
        self.y_thickness = self.get_contour_thickness(self.grid_chevychev())
        self.coordinates = self.create_point_distribution()

    def generate_airfoil(self):
        self.check_name()
        self.make_foil_structure()
        self.send_array()

