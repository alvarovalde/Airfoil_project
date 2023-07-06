import os
import math
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import glob
import matplotlib
import json
matplotlib.use('Qt5Agg')





class Panel:
    """
    Contains information related to a panel.
    """

    def __init__(self, xa, ya, xb, yb):
        """
        Initializes the panel.

        Sets the end-points and calculates the center, length,
        and angle (with the x-axis) of the panel.
        Defines if the panel is on the lower or upper surface of the geometry.
        Initializes the source-sheet strength, tangential velocity,
        and pressure coefficient to zero.

        Parameters
        ----------
        xa: float
            x-coordinate of the first end-point.
        ya: float
            y-coordinate of the first end-point.
        xb: float
            x-coordinate of the second end-point.
        yb: float
            y-coordinate of the second end-point.
        """
        self.xa, self.ya = xa, ya
        self.xb, self.yb = xb, yb

        self.xc, self.yc = (xa + xb) / 2, (ya + yb) / 2  # control-point (center-point)
        self.length = math.sqrt((xb - xa) ** 2 + (yb - ya) ** 2)  # length of the panel

        # orientation of the panel (angle between x-axis and panel's normal)
        if xb - xa <= 0.0:
            self.beta = math.acos((yb - ya) / self.length)
        elif xb - xa > 0.0:
            self.beta = math.pi + math.acos(-(yb - ya) / self.length)

        # location of the panel
        if self.beta <= math.pi:
            self.loc = 'upper'
        else:
            self.loc = 'lower'

        self.sigma = 0.0  # source strength
        self.vt = 0.0  # tangential velocity
        self.cp = 0.0  # pressure coefficient




class Panel_method:
    def __init__(self):
        # aiforil main data
        self.name = None  # NACA####
        self.data = None  # np array with (2,n) points
        self.raw_coordinates = None  # TRANSPOSED DATA FOR EASIER VISUALIZATION
        self.x = None  # np array with (1,n) points
        self.y = None

    def get_airfoil_metrics(self):
        '''
        getting the metrics for the airfoil with the built-in generator, after generating ethe data.
        :return: safe formated information into the object
        '''

        try:
            with open('FoilToAnalize/foil.json', 'r') as file:
                json_data = json.load(file)

                self.name = json_data['name']
                self.camber_builder = json_data['camberline']
                self.n_points = json_data['number of points']

                if self.camber_builder:
                    self.camber_line = np.array(json_data['camber points']).transpose()
                    #print(self.camber_line)
                self.data = np.array(json_data['points']).transpose()

        except SyntaxError as err:
            print(err)
            pass

        self.x = self.data[0]
        self.y = self.data[1]
        self.raw_coordinates = self.data.transpose()
        #np.array((x, y)) #array of shape (2,n): two down and n right

    def define_panels(self):


        x_ends = np.copy(self.x)  # projection of the x-coord on the surface
        y_ends = np.copy(self.y)  # initialization of the y-coord Numpy array


        self.panels = np.empty(self.n_points - 1, dtype=object)
        for i in range(self.n_points -1):
            self.panels[i] = Panel(x_ends[i], y_ends[i], x_ends[i + 1], y_ends[i + 1])

        #print(self.panels)
        return self.panels



    def plot(self):
      # plot the geometry
        width = 10
        plt.figure(figsize=(width, width))
        plt.grid()
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.scatter(self.x, self.y, color='k', linestyle='-', s=0.5)
        plt.axis('scaled')
        plt.xlim(-0.1, 1.1)
        plt.ylim(-0.1, 0.3)
        plt.show()

    def plot_panels(self):

        # plot the geometry and the panels
        width = 10
        plt.figure(figsize=(width, width))
        plt.grid()
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.plot(self.x, self.y, color='k', linestyle='-', linewidth=2)
        plt.plot(np.append([panel.xa for panel in self.panels], self.panels[0].xa),
                 np.append([panel.ya for panel in self.panels], self.panels[0].ya),
                 linestyle='-', linewidth=1, marker='o', markersize=6, color='#CD2305')
        plt.axis('scaled')
        plt.xlim(-0.1, 1.1)
        plt.ylim(-0.1, 0.3)
        plt.show()

pm = Panel_method()
pm.get_airfoil_metrics()
pm.plot()
pm.define_panels()
pm.plot_panels()





class Freestream:
    """
    Freestream conditions.
    """

    def __init__(self, u_inf=1.0, alpha=0.0):
        """
        Sets the freestream speed and angle (in degrees).

        Parameters
        ----------
        u_inf: float, optional
            Freestream speed;
            default: 1.0.
        alpha: float, optional
            Angle of attack in degrees;
            default 0.0.
        """
        self.u_inf = u_inf
        self.alpha = numpy.radians(alpha)  # degrees to radians