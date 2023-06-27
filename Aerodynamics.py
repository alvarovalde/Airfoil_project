import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Qt5Agg')
class Space_domain:
    ''''
    This class creates a virtual 2D space to plot areodynamics,
    it is based on linear, non-viscous, airfoil theory.
    you can create freestreams, sources,sinks, douvlets and vortexes and plot them.
    '''
    def __init__(self,num_of_points_per_direction=200):
        self.N = num_of_points_per_direction    # resolution of grid will be (N,N). Unstable below N = 50
        self.X = None                           # X position of points in grid: (N,N) array of x coordinates
        self.Y =None                            # Y position of points in grid: (N,N) array of y coordinates
        self.u = None                           # X comp of speed vector grid: (N,N) array of u component
        self.v = None                           # Y comp of speed vector grid: (N,N) array of v component
        self.u_inf =None                        # freestream speed, float
        self.cp = None                          # pressure coefficient: (N,N) array of C_p at every point
        self.s_s_pos = []                       # positions of sources,sinks,... , list of tuples
        self.psi = 0                            #

        #plotting parameters
        self.x_start = -2                       #
        self.x_end = 2                          #
        self.y_start = -1                       #
        self.y_end = 1                          #


    def create_grid(self):
        x = np.linspace(self.x_start - 1,self.x_end +11,self.N)
        y = np.linspace(self.y_start - 1,self.y_end +1,self.N)
        self.X,self.Y = np.meshgrid(x,y)
        self.u = np.zeros_like(self.X,dtype=float)
        self.v = np.zeros_like(self.Y,dtype=float)

    def create_freestream(self, freestream_speed=1.0,angle=0.0):
        self.u_inf = freestream_speed
        self._inf = self.u_inf * np.ones(self.X.shape,dtype=float)
        a = angle * np.pi /180
        self.psi = self.psi + self.u_inf * self.Y
        self.u = self.u + self.u_inf*np.cos(a)
        self.v = self.v + self.u_inf * np.sin(a)

    def create_source_flow(self,strength=1.0,x_sr= 0.0,y_sr= 0.0):
        # compute the velocity field on the mesh grid
        self.u = self.u + (strength/(2*np.pi))*((self.X -x_sr)/((self.X -x_sr)**2 + (self.Y -y_sr)**2))
        self.v = self.v + (strength / (2 * np.pi) *((self.Y - y_sr) / ((self.X - x_sr) ** 2 + (self.Y - y_sr) ** 2)))
        self.s_s_pos.append((x_sr,y_sr))
        self.psi = self.psi + strength / (2 * np.pi) * np.arctan2((self.Y - y_sr), (self.X - x_sr))

    def create_sink_flow(self,strength=-1.0,x_snk= 0.0,y_snk= 0.0):
        # compute the velocity field on the mesh grid
        self.u = self.u + (strength/(2*np.pi))*((self.X -x_snk)/((self.X -x_snk)**2 + (self.Y -y_snk)**2))
        self.v = self.v + (strength / (2 * np.pi) *((self.Y - y_snk) / ((self.X - x_snk) ** 2 + (self.Y - y_snk) ** 2)))
        self.s_s_pos.append((x_snk, y_snk))
        self.psi = self.psi + strength / (2 * np.pi) * np.arctan2((self.Y - y_snk), (self.X - x_snk))

    def create_vortex(self,gamma=5.0,x_vor=0.0,y_vor=0.0):
        self.u = self.u + +gamma / (2 * math.pi) * (self.Y - y_vor) / ((self.X - x_vor) ** 2 + (self.Y - y_vor) ** 2)
        self.v = self.v + -gamma / (2 * math.pi) * (self.X - x_vor) / ((self.X - x_vor) ** 2 + (self.Y - y_vor) ** 2)
        self.psi = self.psi + gamma / (4 * math.pi) * np.log((self.X - x_vor) ** 2 + (self.Y - y_vor) ** 2)
        self.s_s_pos.append((x_vor, y_vor))


    def create_doublet(self,kappa=1, xpos=0.0,ypos=0.0 ):
        self.u = self.u +(-kappa / (2 * math.pi) *
             ((self.X - xpos) ** 2 - (self.Y - ypos) ** 2) /
             ((self.X - xpos) ** 2 + (self.Y - ypos) ** 2) ** 2)
        self.v = self.v + (-kappa / (2 * math.pi) *
             2 * (self.X - xpos) * (self.Y - ypos) /
             ((self.X - xpos) ** 2 + (self.Y - ypos) ** 2) ** 2)

        self.psi = -kappa / (2 * math.pi) * (self.Y - ypos) / ((self.X - xpos) ** 2 + (self.Y - ypos) ** 2)
        self.s_s_pos.append((xpos, ypos))
        # calculate the cylinder radius
        self.R = np.sqrt(kappa/(2*np.pi*self.u_inf))

    def get_cp(self):
        try:
            self.cp = 1.0 - (self.u ** 2 + self.v ** 2) / self.u_inf ** 2
        except TypeError as err:
            print('You need to create a freestream to get a Cp, try creating a uniform flow!',err)
            quit()

    def plot_field(self):
        ax = plt.subplot()
        ax.axis([self.x_start, self.x_end, self.y_start, self.y_end])
        ax.streamplot(self.X, self.Y, self.u, self.v, density=5, linewidth=1, arrowsize=1,arrowstyle='->')
        ax.contour(self.X,self.Y,self.psi,levels=[0.], colors='#CD2305', linewidths=2, linestyles='solid')

        circle = plt.Circle((0, 0), radius=self.R, color='#CD2305', alpha=0.5)
        plt.gca().add_patch(circle)
        try:
            x,y = np.array(self.s_s_pos).T
            ax.scatter(x,y,color='red')
            print(so.s_s_pos)
        except ValueError as err:
            print('No points to plot',err)





    def plot_cp(self):
        ax = plt.subplot()
        ax.axis([self.x_start, self.x_end, self.y_start, self.y_end])
        contf = plt.contourf(self.X, self.Y, self.cp,levels=np.linspace(-2.0, 1.0, 100), extend='both')
        cbar =plt.colorbar(contf)
        cbar.set_label('$C_p$', fontsize=16)
        cbar.set_ticks([-2.0, -1.0, 0.0, 1])
        plt.show()

if __name__ == '__main__':
    so =Space_domain()
    so.create_grid()
    so.create_freestream(0.5,0)
    #so.create_vortex()
    so.create_doublet()
    so.get_cp()

    so.plot_field()
    so.plot_cp()
