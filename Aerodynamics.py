import math
import numpy as np
import matplotlib.pyplot as plt
class Space_domain:
    def __init__(self,num_of_points_per_direction=50):
        self.N = num_of_points_per_direction
        self.x_start = -2
        self.x_end = 2
        self.y_start = -2
        self.y_end = 2
        self.X = None
        self.Y =None # Y spoition of point (
        self.u = None   # X comp of speed vector
        self.v =None # Y comp of speed vector


    def create_grid(self):
        x = np.linspace(self.x_start,self.x_end,self.N)
        y = np.linspace(self.y_start,self.y_end,self.N)
        self.X,self.Y = np.meshgrid(x,y)
        self.u,self.v = np.zeros_like(np.meshgrid(self.N,self.N))


    def create_source_flow(self,strength=5,x_sr= -1,y_sr= 0):
        # compute the velocity field on the mesh grid
        self.u = self.u +(strength / (2 * math.pi) * ((self.X - x_sr) / ((self.X - x_sr) ** 2 + (self.Y - y_sr) ** 2)))
        self.v = self.v + (strength / (2 * math.pi) *((self.Y - y_sr) / ((self.X - y_sr) ** 2 + (self.Y - y_sr) ** 2)))

        width = 10
        fig, ax = plt.subplots(figsize=(width, (self.y_end - self.y_start) / (self.x_end - self.x_start) * width))
        # ax.axis([-0.25, 1.25, -1, 1])
        ax.streamplot(self.X, self.Y, self.u, self.v, density=2, linewidth=1, arrowsize=2, arrowstyle='->')
        ax.scatter(x_sr,y_sr, color='#CD2305', s=40, marker='o')
        plt.show()

so =Space_domain()
so.create_grid()
so.create_source_flow(1,1,0)
print(so.u,so.v)
