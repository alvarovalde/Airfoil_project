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
from multiprocessing import Process,cpu_count

class Naca4Creator:
    def __init__(self,NACA=2309,n_points=5000):
        self.NACA = NACA  #this is name of naca before processing
        self.n_points = n_points
        self.NACA_designator = 4

        self.name = 'NACA ' + str(NACA)
        self.max_camber = None #m
        self.loc_mc = None #p
        self.thickness = None #t

        #point distributions
        self.y_camberline = None
        self.theta_contour = None
        self.y_thickness = None

        # final distributions
        self.coordinates = None #coordinates of the thickness distribution
        self.camberline_array =None
        self.file_array = None



    def check_name(self):
        self.NACA = str(self.NACA)
        if len(self.NACA) == self.NACA_designator:
                self.max_camber = float(self.NACA[0])/100
                self.loc_mc = float(self.NACA[1])/10
                self.thickness = float(self.NACA[2:]) / 100
        if len(self.NACA) != self.NACA_designator:
              print(f'This NACA only asks for a {self.NACA_designator} number')

    def grid_uniform(self):
        return np.linspace(0,1,self.n_points)

    def grid_chevychev(self):
        xi = np.cos(np.pi/(2*(self.n_points+1))*(2*np.linspace(0,self.n_points+1,self.n_points+1)-1))
        return ((-xi +1)/2)*(1-0) + 0

    def camberline_function(self,x):
        yc = np.zeros(x.shape)
        yc[x < self.loc_mc] = self.max_camber * x[x < self.loc_mc]/self.loc_mc**2 * (2*self.loc_mc - x[x < self.loc_mc])
        yc[x >= self.loc_mc] = self.max_camber * (1-x[x >= self.loc_mc]) / (1 - self.loc_mc)**2 * (1 + x[x >= self.loc_mc] - 2*self.loc_mc)
        return yc

        # getting theta function for the plotting of the contour
    def theta_funct(self,x):
        if self.max_camber > 0:
            dydx = np.zeros(x.shape)
            dydx[x < self.loc_mc] = 2 * self.max_camber / self.loc_mc ** 2 * (self.loc_mc - x[x < self.loc_mc])
            dydx[x >= self.loc_mc] = 2 * self.max_camber / (1 - self.loc_mc) ** 2 * (self.loc_mc - x[x >= self.loc_mc])
            theta = np.arctan(dydx)
        else:
            theta = 0
        return theta
    def get_contour_thickness(self,x):
        # getting the thickness of the contour
        a0 = 0.2969
        a1 = -0.1260
        a2 = -0.3516
        a3 = 0.2843
        a4 = -0.1015
        yt = self.thickness * 5 * (a0 * np.sqrt(x) + a1 * x + a2 * x**2 +a3 * x**3 + a4 * x**4)
        return yt

    def create_point_distribution(self):
        x_chev = self.grid_chevychev()

        xu = x_chev - self.y_thickness * np.sin(self.theta)
        xl = x_chev + self.y_thickness * np.sin(self.theta)

        yu = self.y_camberline + self.y_thickness * np.cos(self.theta)
        yl = self.y_camberline - self.y_thickness * np.cos(self.theta)

        cor_up = np.array((xu, yu)).transpose()
        cor_down = np.array((xl, yl)).transpose()
        self.coordinates = np.concatenate((cor_down, cor_up))  # (102,1)
        return self.coordinates

    def create_final_array(self):
        cor_camb = np.array((self.grid_chevychev(),self.y_camberline)).transpose()  # (54,1)
        n_elem = self.coordinates.shape[0]
        nan_array = np.full((n_elem - cor_camb.shape[0], 2), np.nan)
        self.camberline_array = np.concatenate((cor_camb, nan_array))
        self.file_array = np.concatenate((self.coordinates, self.camberline_array), axis=1)

    def send_array(self):

        # clean directory
        clean_directory()

        # create new file
        np.savetxt(f"FoilToAnalize/{self.name}.dat", self.file_array, delimiter=' ')


    def make_airfoil(self):
        self.check_name()
        self.y_camberline = self.camberline_function(self.grid_chevychev())
        self.theta = self.theta_funct(self.grid_chevychev())
        self.y_thickness = self.get_contour_thickness(self.grid_chevychev())
        self.coordinates = self.create_point_distribution()
        self.create_final_array()
        self.send_array()

