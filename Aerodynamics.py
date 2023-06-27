import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Qt5Agg')
class Space_domain:
    def __init__(self,num_of_points_per_direction=200):
        self.N = num_of_points_per_direction
        self.x_start = -2
        self.x_end = 2
        self.y_start = -1
        self.y_end = 1
        self.X = None
        self.Y =None # Y spoition of point (
        self.u = None   # X comp of speed vector
        self.v = None # Y comp of speed vector
        self.u_inf =None
        self.cp = None
        self.s_s_pos = []
        self.psi = 0


    def create_grid(self):
        x = np.linspace(self.x_start,self.x_end,self.N)
        y = np.linspace(self.y_start,self.y_end,self.N)
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

        x,y = np.array(self.s_s_pos).T
        ax.scatter(x,y,color='red')
    def create_doublet(self,strength=1, xpos=0.0,ypos=0.0 ):
        self.create_source_flow(strength,xpos -0.5, ypos)
        self.create_sink_flow(-strength,xpos+0.5,ypos)


    def plot_cp(self):
        ax = plt.subplot()
        ax.axis([self.x_start, self.x_end, self.y_start, self.y_end])
        contf = plt.contourf(self.X, self.Y, self.cp,levels=np.linspace(-2.0, 1.0, 100), extend='both')
        cbar =plt.colorbar(contf)
        cbar.set_label('$C_p$', fontsize=16)
        cbar.set_ticks([-2.0, -1.0, 0.0, 1])
        plt.show()

so =Space_domain()
so.create_grid()
so.create_freestream(1,0)
so.create_vortex(4)
so.create_doublet()
so.get_cp()
print(so.s_s_pos)
so.plot_field()
so.plot_cp()
