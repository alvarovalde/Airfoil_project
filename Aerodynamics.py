import numpy as np
import matplotlib.pyplot as plt
import matplotlib, cmath, multiprocessing
matplotlib.use('TkAgg')






class SpaceDomain:
    """
    This class creates a virtual 2D space to plot areodynamics,
    it is based on linear, non-viscous, airfoil theory.
    you can create freestreams, sources,sinks, doublets and vortexes and plot them.
    """
    def __init__(self, num_of_points_per_direction=200):
        # grid parameters
        self.N = num_of_points_per_direction    # resolution of grid will be (N,N). Unstable below N = 50
        self.X = None                           # X position of points in grid: (N,N) array of x coordinates
        self.Y = None                            # Y position of points in grid: (N,N) array of y coordinates
        self.u = None                           # X comp of speed vector grid: (N,N) array of u component
        self.v = None                           # Y comp of speed vector grid: (N,N) array of v component
        self.cp = None                          # pressure coefficient: (N,N) array of C_p at every point
        self.psi = 0                            # streamline function,  (N,N) array of values

        # flow parameters and Points Of Interest
        self.u_inf = None                       # freestream speed, float
        self.s_s_pos = []                       # positions of sources,sinks,... , list of tuples
        self.st_points = []                     # positions of stagnation points, list of tuples

        # plotting parameters
        self.x_start = -4                       # start of view window in x direction
        self.x_end = 4                          # end of view window in x direction
        self.y_start = -4                       # start of view window in y direction
        self.y_end = 4                          # end of view window in y direction


        # variable values (depending on current construction)
        self.gamma = 0
        self.kappa = 0


        # plotting values
        self.vec_density = 4
        self.R = 0

    def create_grid(self):
        """
        This function creates the grid first generating a linear space (N points),
        and then creates a mesh of these points (N,N) ands divides it into its
         X and Y coordinates. (+/- 1 for plotting neatness). Initializes the shape of u and v.
        """
        x = np.linspace(self.x_start - 1, self.x_end + 1, self.N)
        y = np.linspace(self.y_start - 1, self.y_end + 1, self.N)
        self.X, self.Y = np.meshgrid(x, y)
        self.u = np.zeros_like(self.X, dtype=float)
        self.v = np.zeros_like(self.Y, dtype=float)

    def create_freestream(self, freestream_speed=1.0, angle=0.0):
        """
        Creates linear flow behaviour and ads it to current flow state.

        @param freestream_speed: the speed of the created flow [-]
        @param angle: angle of attack, degrees.
        """

        self.u_inf = freestream_speed
        a = angle * np.pi / 180
        self.psi = self.psi + self.u_inf * self.Y
        self.u += self.u_inf*np.cos(a)
        self.v += self.u_inf * np.sin(a)



    def create_source_flow(self, strength=1.0, x_sr=0.0, y_sr=0.0):
        """
        Creates a source behaviour and ads it to current flow state.
        Append the source position to s_s list.

        @param strength: strength of the source
        @param x_sr: x position of source flow
        @param y_sr: y position of source flow
        """
        strength = abs(strength)
        self.u += (strength/(2*np.pi))*((self.X - x_sr)/((self.X - x_sr)**2 + (self.Y - y_sr)**2))
        self.v += self.v + (strength / (2 * np.pi) * ((self.Y - y_sr) / ((self.X - x_sr) ** 2 + (self.Y - y_sr) ** 2)))
        self.psi = self.psi + strength / (2 * np.pi) * np.arctan2((self.Y - y_sr), (self.X - x_sr))
        self.s_s_pos.append((x_sr, y_sr))

    def create_sink_flow(self, strength=-1.0, x_snk=0.0, y_snk=0.0):
        """
        Creates a sink behaviour and ads it to current flow state.
        Append the sink position to s_s list.
        @param strength: strength of the source
        @param x_snk: x position of sink flow
        @param y_snk: y position of sink flow
        """
        strength = - abs(strength)
        self.u = self.u + (strength/(2*np.pi))*((self.X - x_snk)/((self.X - x_snk)**2 + (self.Y - y_snk)**2))
        self.v = self.v + (strength / (2 * np.pi) *
                           ((self.Y - y_snk) / ((self.X - x_snk) ** 2 + (self.Y - y_snk) ** 2)))
        self.psi = self.psi + strength / (2 * np.pi) * np.arctan2((self.Y - y_snk), (self.X - x_snk))
        self.s_s_pos.append((x_snk, y_snk))

    def create_vortex(self, gamma=5.0, x_vor=0.0, y_vor=0.0):
        """
         Creates a vortex behaviour and ads it to current flow state.
        Append the vortex position to s_s list.
        @param gamma: strength of vortex. (positive rotation is clockwise)
        @param x_vor: x position of the vortex
        @param y_vor: y position of the vortex
        """
        self.gamma = gamma
        self.u = self.u + +gamma / (2 * np.pi) * (self.Y - y_vor) / ((self.X - x_vor) ** 2 + (self.Y - y_vor) ** 2)
        self.v = self.v + -gamma / (2 * np.pi) * (self.X - x_vor) / ((self.X - x_vor) ** 2 + (self.Y - y_vor) ** 2)
        self.psi = self.psi + gamma / (4 * np.pi) * np.log((self.X - x_vor) ** 2 + (self.Y - y_vor) ** 2)
        self.s_s_pos.append((x_vor, y_vor))


    def create_doublet(self, kappa=1, x_db=0.0, y_db=0.0):
        """
        Creates a doublet behaviour and ads it to current flow state.
        Append the doublet position to s_s list.
        Create the radius of the doublet for plotting.
        @param kappa: strength of the doublet
        @param x_db: x position of center of doublet
        @param y_db: y position of center of doublet
        """
        self.kappa = kappa

        self.u = self.u + (-self.kappa / (2 * np.pi) *
            ((self.X - x_db) ** 2 - (self.Y - y_db) ** 2) /
            ((self.X - x_db) ** 2 + (self.Y - y_db) ** 2) ** 2)

        self.v = self.v + (-self.kappa / (2 * np.pi) *
            2 * (self.X - x_db) * (self.Y - y_db) /
            ((self.X - x_db) ** 2 + (self.Y - y_db) ** 2) ** 2)

        self.psi = -self.kappa / (2 * np.pi) * (self.Y - y_db) / ((self.X - x_db) ** 2 + (self.Y - y_db) ** 2)

        self.s_s_pos.append((x_db, y_db))

        # calculate the cylinder radius
        self.R = np.sqrt(kappa/(2*np.pi*self.u_inf))

    def create_doublet_with_vortex(self, kappa=1.0, gamma=1):
        """
        Creates a doublet (with a vortex inside) behaviour, esentialy a rotating solid body,
        and ads it to current flow state.

        Append the doublet position to s_s list.
        Create the radius of the doublet for plotting.
        @param kappa: strength of the doublet
        @param gamma: strength of vortex, normally same as kappa. (positive rotation is clockwise)
                    - gamma = 0: no rotation; two stagnation points (0,pi)
                    - gamma < 4 pi * R * V : two real stagnation points
                    - gamma = 4 pi * R * V : one stagnation point
                    - gamma > 4 pi * R * V : two imaginary stagnation points
        @param x_vdb: (not-variable) x position of ,center of circular object
        @param y_vdb: (not-variable) x position of ,center of circular object
        """

        x_vdb = 0.0
        y_vdb = 0.0
        self.kappa = kappa
        self.gamma = gamma

        self.u = self.u + (-self.kappa / (2 * np.pi) *
             ((self.X - x_vdb) ** 2 - (self.Y - y_vdb) ** 2) /
             ((self.X - x_vdb) ** 2 + (self.Y - y_vdb) ** 2) ** 2)


        self.v = self.v + (-self.kappa / (2 * np.pi) * 2 * (self.X - x_vdb) * (self.Y - y_vdb) /
             ((self.X - x_vdb) ** 2 + (self.Y - y_vdb) ** 2) ** 2)

        self.psi = -kappa / (2 * np.pi) * (self.Y - y_vdb) / ((self.X - x_vdb) ** 2 + (self.Y - y_vdb) ** 2)

        # create vortex
        self.create_vortex(self.gamma, x_vdb, y_vdb)

        # calculate the cylinder radius
        self.R = np.sqrt(kappa/(2*np.pi*self.u_inf))
        print(self.R*4*np.pi*self.u_inf)

    def get_cp(self):
        """
        Creates a scalar field of the preasure coefficients, which can then be plotted.
        @return:returns the Cp object
        """
        print(np.shape(self.u))
        try:
            self.cp = 1.0 - (self.u ** 2 + self.v ** 2) / self.u_inf ** 2
        except TypeError as err:
            print('You need to create a freestream to get a Cp, try creating a uniform flow!',err)
            quit()

        # calculate the surface tangential velocity on the cylinder
        theta = np.linspace(0.0, 2 * np.pi, 100)
        u_theta = -2 * self.u_inf * np.sin(theta) - self.gamma / (2 * np.pi * self.R)

        # compute the surface pressure coefficient
        cp = 1.0 - (u_theta / self.u_inf) ** 2

        return cp

    def get_stagnation_points(self):
        """
        gets the positions of the stagnation points and appends them to the list of points.
        """
        try:
            point1 = (+cmath.sqrt(self.R ** 2 - (self.gamma / (4 * cmath.pi * self.u_inf)) ** 2) * 1j,
                      -self.gamma / (4 * cmath.pi * self.u_inf))
            self.st_points.append(point1)
        except RuntimeWarning as err:
            print('gamma is too large')
            pass
        try:
            point2 = (-cmath.sqrt(self.R ** 2 - (self.gamma / (4 * cmath.pi * self.u_inf)) ** 2) * 1j,
                      -self.gamma / (4 * cmath.pi * self.u_inf))
            self.st_points.append(point2)
        except RuntimeWarning as err:
            print('gamma is too large')
            pass
        print(self.st_points)



    def plot_field(self,pressure_field = True):
        """
        plotting the field with Matplotlib.
        """

        # use for external window
        matplotlib.use('Qt5Agg')

        #start figure
        fig,ax = plt.subplots()
        ax.set_aspect('equal')
        ax.axis([self.x_start, self.x_end, self.y_start, self.y_end])

        #plot the vector field
        ax.streamplot(self.X, self.Y, self.u, self.v, density=self.vec_density, linewidth=1, arrowsize=1,arrowstyle='->')

        #plot the circle
        circle = plt.Circle((0, 0), radius=self.R, color='black', alpha=0.5)

        # Try to plot the source,sink,... positions
        try:
            x,y = np.array(self.s_s_pos).T
            ax.scatter(x, y, color='red')
            print(self.s_s_pos)
        except ValueError as err:
            print('No points to plot', err)

        #choose wether to plot the pressure field
        if pressure_field:
            contf = plt.contourf(self.X, self.Y, self.cp, levels=np.linspace(-2.0, 1.0, 100), extend='both')
            cbar = plt.colorbar(contf)
            cbar.set_label('$C_p$', fontsize=16)

        #try to plot the stagnation points
        try:
            x,y = np.array(self.st_points).T
            ax.scatter(x,y,color='g')
        except:
            pass

        #add the circle to the plot
        ax.add_patch(circle)

        #show the plot
        plt.show()

    def plot_cp(self,cp):
        size = 6
        plt.figure(figsize= (5,5))
        plt.grid(True)
        plt.xlabel(r'$\theta$', fontsize=18)
        plt.ylabel('$C_p$', fontsize=18)
        plt.xlim(0, 2*np.pi)
        theta = np.linspace(0.0, 2 * np.pi, 100)
        plt.plot(theta, cp, color='#CD2305', linewidth=2, linestyle='-')
        plt.show()

        matplotlib.use('Qt5Agg')
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        ax.axis([self.x_start, self.x_end, self.y_start, self.y_end])
        contf = plt.contourf(self.X, self.Y, self.psi, levels=np.linspace(-2.0, 1.0, 100), extend='both')
        cbar = plt.colorbar(contf)
        ax.contour(self.X,self.Y,self.psi,levels=[0.], colors='#CD2305', linewidths=2, linestyles='solid')
        plt.show()

if __name__ == '__main__':

    #create an instance of the Space_domain and generate the grid
    sd =SpaceDomain()
    sd.create_grid()

    #create a freestream
    sd.create_freestream(0.5,0)

    #create a doublet with a vortex
    sd.create_doublet_with_vortex()



    #sd.get_stagnation_points()
    cp = sd.get_cp()

    p1 = multiprocessing.Process(target=sd.plot_cp(cp))
    p2 = multiprocessing.Process(target=sd.plot_field(pressure_field=True))

    # Start both processes
    p1.start()
    p2.start()

    # Wait for both processes to finish
    p1.join()
    p2.join()