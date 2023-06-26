import math
import numpy as np
import matplotlib.pyplot as plt
class Space_domain:
    def __init__(self,num_of_points_per_direction=100):
        self.N = num_of_points_per_direction
        self.x_start = -3
        self.x_end = 3
        self.y_start = -3
        self.y_end = 3
        self.X = None
        self.Y =None # Y spoition of point (
        self.u = None   # X comp of speed vector
        self.v =None # Y comp of speed vector
        self.s_s_pos = []


    def create_grid(self):
        x = np.linspace(self.x_start,self.x_end,self.N)
        y = np.linspace(self.y_start,self.y_end,self.N)
        self.X,self.Y = np.meshgrid(x,y)
        self.u,self.v = np.zeros_like(np.meshgrid(self.N,self.N))

    def create_freestream(self, freestream_speed=1):
        u_inf = freetream_speed

    def create_source_flow(self,strength=1,x_sr= -1,y_sr= 0):
        # compute the velocity field on the mesh grid
        self.u = self.u + (strength/(2*np.pi))*((self.X -x_sr)/((self.X -x_sr)**2 + (self.Y -y_sr)**2))
        self.v = self.v + (strength / (2 * np.pi) *((self.Y - y_sr) / ((self.X - x_sr) ** 2 + (self.Y - y_sr) ** 2)))
        self.s_s_pos.append((x_sr,y_sr))
    def create_sink_flow(self,strength=-1,x_snk= 1,y_snk= 0):
        # compute the velocity field on the mesh grid
        self.u = self.u + (strength/(2*np.pi))*((self.X -x_snk)/((self.X -x_snk)**2 + (self.Y -y_snk)**2))
        self.v = self.v + (strength / (2 * np.pi) *((self.Y - y_snk) / ((self.X - x_snk) ** 2 + (self.Y - y_snk) ** 2)))
        self.s_s_pos.append((x_snk, y_snk))

    def plot_field(self):
        fig, ax = plt.subplots()
        ax.axis([self.x_start, self.x_end, self.y_start, self.y_end])
        ax.streamplot(self.X, self.Y, self.u, self.v, density=2, linewidth=1, arrowsize=2)
        plt.show()

so =Space_domain()
so.create_grid()
so.create_source_flow(1,1,0)
so.create_sink_flow(-1,-1,0)
so.plot_field()