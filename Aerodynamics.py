import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Qt5Agg')
class Space_domain:
    def __init__(self,num_of_points_per_direction=200):
        self.N = num_of_points_per_direction
        self.x_start = -4
        self.x_end = 4
        self.y_start = -2
        self.y_end = 2
        self.X = None
        self.Y =None # Y spoition of point (
        self.u = None   # X comp of speed vector
        self.v =None # Y comp of speed vector
        self.s_s_pos = []
        #self.psi = None


    def create_grid(self):
        x = np.linspace(self.x_start,self.x_end,self.N)
        y = np.linspace(self.y_start,self.y_end,self.N)
        self.X,self.Y = np.meshgrid(x,y)
        self.u = np.zeros_like(self.X,dtype=float)
        self.v = np.zeros_like(self.Y,dtype=float)

    def create_freestream(self, freestream_speed=1,angle=0):
        u_inf = freestream_speed
        u_inf = u_inf * np.ones(self.X.shape,dtype=float)
        a = angle * np.pi /180
        #self.psi= u_inf * self.Y
        self.u = self.u + u_inf*np.cos(a)
        self.v = self.v + u_inf * np.sin(a)

    def create_source_flow(self,strength=1,x_sr= -1,y_sr= 0):
        # compute the velocity field on the mesh grid
        self.u = self.u + (strength/(2*np.pi))*((self.X -x_sr)/((self.X -x_sr)**2 + (self.Y -y_sr)**2))
        self.v = self.v + (strength / (2 * np.pi) *((self.Y - y_sr) / ((self.X - x_sr) ** 2 + (self.Y - y_sr) ** 2)))
        self.s_s_pos.append((x_sr,y_sr))
        #self.psi = strength / (2 * np.pi) * np.arctan2((self.Y - y_sr), (self.X - x_sr))

    def create_sink_flow(self,strength=-1,x_snk= 1,y_snk= 0):
        # compute the velocity field on the mesh grid
        self.u = self.u + (strength/(2*np.pi))*((self.X -x_snk)/((self.X -x_snk)**2 + (self.Y -y_snk)**2))
        self.v = self.v + (strength / (2 * np.pi) *((self.Y - y_snk) / ((self.X - x_snk) ** 2 + (self.Y - y_snk) ** 2)))
        self.s_s_pos.append((x_snk, y_snk))

    def plot_field(self):
        fig, ax = plt.subplots()
        ax.axis([self.x_start, self.x_end, self.y_start, self.y_end])
        ax.streamplot(self.X, self.Y, self.u, self.v, density=5, linewidth=1, arrowsize=1,arrowstyle='->')
        #ax.contour(self.X,self.Y,self.psi,levels=[-1 / 2,1 / 2], colors='#CD2305', linewidths=2, linestyles='solid')
        plt.show()

so =Space_domain()
so.create_grid()
so.create_freestream(0.1,0)
so.create_source_flow(1,-1,0)
so.create_sink_flow(-1,1,0)
so.plot_field()
