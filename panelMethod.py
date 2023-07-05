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
                    print(self.camber_line)
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

        print(self.panels)
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


quit()

























class Freestream:
    """
    Freestream conditions.
    """

    def __init__(self, u_inf=1.0, alpha=0.0):
        """
        Sets the freestream speed and angle (with the x-axis).

        Parameters
        ----------
        u_inf: float, optional
            Freestream speed;
            default: 1.0.
        alpha: float, optional
            Angle of attack in degrees;
            default: 0.0.
        """
        self.u_inf = u_inf
        self.alpha = np.radians(alpha)  # degrees --> radians


# define and creates the object freestream
u_inf = 1.0                            # freestream spee
alpha = 0.0                            # angle of attack (in degrees)
freestream = Freestream(u_inf, alpha)  # instantiation of the object freestream


def integral(x, y, panel, dxdz, dydz):
    """
    Evaluates the contribution of a panel at one point.

    Parameters
    ----------
    x: float
        x-coordinate of the target point.
    y: float
        y-coordinate of the target point.
    panel: Panel object
        Source panel which contribution is evaluated.
    dxdz: float
        Derivative of x in the z-direction.
    dydz: float
        Derivative of y in the z-direction.

    Returns
    -------
    Integral over the panel of the influence at the given target point.
    """

    def integrand(s):
        return (((x - (panel.xa - math.sin(panel.beta) * s)) * dxdz +
                 (y - (panel.ya + math.cos(panel.beta) * s)) * dydz) /
                ((x - (panel.xa - math.sin(panel.beta) * s)) ** 2 +
                 (y - (panel.ya + math.cos(panel.beta) * s)) ** 2))

    return integrate.quad(integrand, 0.0, panel.length)[0]


def build_matrix(panels):
    """
    Builds the source matrix.

    Parameters
    ----------
    panels: 1D array of Panel object
        The source panels.

    Returns
    -------
    A: 2D Numpy array of floats
        The source matrix (NxN matrix; N is the number of panels).
    """
    N = len(panels)
    A = np.empty((N, N), dtype=float)
    np.fill_diagonal(A, 0.5)

    for i, p_i in enumerate(panels):
        for j, p_j in enumerate(panels):
            if i != j:
                A[i, j] = 0.5 / math.pi * integral(p_i.xc, p_i.yc, p_j,
                                                   math.cos(p_i.beta),
                                                   math.sin(p_i.beta))

    return A


def build_rhs(panels, freestream):
    """
    Builds the RHS of the linear system.

    Parameters
    ----------
    panels: 1D array of Panel objects
        The source panels.
    freestream: Freestream object
        The freestream conditions.

    Returns
    -------
    b: 1D Numpy array of floats
        RHS of the linear system.
    """
    b = np.empty(len(panels), dtype=float)

    for i, panel in enumerate(panels):
        b[i] = -freestream.u_inf * math.cos(freestream.alpha - panel.beta)

    return b



A = build_matrix(panels)           # compute the singularity matrix
b = build_rhs(panels, freestream)  # compute the freestream RHS


# solve the linear system
sigma = np.linalg.solve(A, b)

for i, panel in enumerate(panels):
    panel.sigma = sigma[i]


def get_tangential_velocity(panels, freestream):
    """
    Computes the tangential velocity on the surface of the panels.

    Parameters
    ---------
    panels: 1D array of Panel objects
        The source panels.
    freestream: Freestream object
        The freestream conditions.
    """
    N = len(panels)
    A = np.empty((N, N), dtype=float)
    np.fill_diagonal(A, 0.0)

    for i, p_i in enumerate(panels):
        for j, p_j in enumerate(panels):
            if i != j:
                A[i, j] = 0.5 / math.pi * integral(p_i.xc, p_i.yc, p_j,
                                                   -math.sin(p_i.beta),
                                                   math.cos(p_i.beta))

    b = freestream.u_inf * np.sin([freestream.alpha - panel.beta
                                      for panel in panels])

    sigma = np.array([panel.sigma for panel in panels])

    vt = np.dot(A, sigma) + b

    for i, panel in enumerate(panels):
        panel.vt = vt[i]

# compute the tangential velocity at the center-point of each panel
get_tangential_velocity(panels, freestream)


def get_pressure_coefficient(panels, freestream):
    """
    Computes the surface pressure coefficients on the panels.

    Parameters
    ---------
    panels: 1D array of Panel objects
        The source panels.
    freestream: Freestream object
        The freestream conditions.
    """
    for panel in panels:
        panel.cp = 1.0 - (panel.vt / freestream.u_inf) ** 2


# computes the surface pressure coefficients
get_pressure_coefficient(panels, freestream)


def get_velocity_field(panels, freestream, X, Y):
    """
    Computes the velocity field on a given 2D mesh.

    Parameters
    ---------
    panels: 1D array of Panel objects
        The source panels.
    freestream: Freestream object
        The freestream conditions.
    X: 2D Numpy array of floats
        x-coordinates of the mesh points.
    Y: 2D Numpy array of floats
        y-coordinate of the mesh points.

    Returns
    -------
    u: 2D Numpy array of floats
        x-component of the velocity vector field.
    v: 2D Numpy array of floats
        y-component of the velocity vector field.
    """
    # freestream contribution
    u = freestream.u_inf * math.cos(freestream.alpha) * np.ones_like(X, dtype=float)
    v = freestream.u_inf * math.sin(freestream.alpha) * np.ones_like(X, dtype=float)
    # add the contribution from each source (superposition powers!!!)
    vec_intregral = np.vectorize(integral)
    for panel in panels:
        u += panel.sigma / (2.0 * math.pi) * vec_intregral(X, Y, panel, 1.0, 0.0)
        v += panel.sigma / (2.0 * math.pi) * vec_intregral(X, Y, panel, 0.0, 1.0)

    return u, v

# define a mesh grid
nx, ny = 20, 20  # number of points in the x and y directions
x_start, x_end = -1.0, 2.0
y_start, y_end = -0.3, 0.3
X, Y = np.meshgrid(np.linspace(x_start, x_end, nx),
                      np.linspace(y_start, y_end, ny))

# compute the velocity field on the mesh grid
u, v = get_velocity_field(panels, freestream, X, Y)


# plot the velocity field
width = 10
plt.figure(figsize=(width, width))
plt.xlabel('x', fontsize=16)
plt.ylabel('y', fontsize=16)
plt.streamplot(X, Y, u, v,
                  density=1, linewidth=1, arrowsize=1, arrowstyle='->')
plt.fill([panel.xc for panel in panels],
            [panel.yc for panel in panels],
            color='k', linestyle='solid', linewidth=2, zorder=2)
plt.axis('scaled')
plt.xlim(x_start, x_end)
plt.ylim(y_start, y_end)
plt.title('Streamlines around a NACA 0012 airfoil (AoA = ${}^o$)'.format(alpha),fontsize=16)
plt.show()