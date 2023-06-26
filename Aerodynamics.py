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
        self.Y =None


    def create_grid(self):
        x = np.linspace(self.x_start,self.x_end,self.N)
        y = np.linspace(self.y_start,self.y_end,self.N)
        self.X,self.Y = np.meshgrid(x,y)


class source_flow(Space_domain):
    def __init__(self,strength=5,x_sr= 0,y_sr= 0):
        super().__init__()
        self.strength_source = strength
        self.x_sr = x_sr
        self.y_sr = y_sr
