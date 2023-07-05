from directory_management import clean_directory
import json
import numpy as np
import os
import matplotlib
import math
from scipy import interpolate
matplotlib.use('Qt5Agg')  #Try different backends (such as 'TkAgg', 'Agg', or 'Qt5Agg')
import matplotlib.pyplot as plt
import glob
import PyQt5

class Airfoil_manager:
    """Class to read airfoil data from file (or use predefined airfoil)"""
    def __init__(self):
        #aiforil main data
        self.name = None             # NACA####
        self.data =None              # np array with (2,n) points
        self.raw_coordinates = None  #TRANSPOSED DATA FOR EASIER VISUALIZATION
        self.x = None                # np array with (1,n) points
        self.y = None                # np array with (1,n) points

        #Chord stuff (... a class maybe here in the future...?)
        self.camber_builder = True
        self.leading_edge = None     #tuple of coordinates
        self.trailing_edge = None
        self.chord_length = None
        self.camber_line = None



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

                if self.camber_builder:
                    self.camber_line = np.array(json_data['camber points']).transpose()
                    print(self.camber_line)
                self.data = np.array(json_data['points']).transpose()

        except SyntaxError as err:
            print(err)
            pass



        # finally:
        #     print('There is no foil to analize in the folder')
        #     quit()


        self.x = self.data[0]
        self.y = self.data[1]
        self.raw_coordinates = self.data.transpose()                  #np.array((x, y)) #array of shape (2,n): two down and n right

    def closest_point_to_origin(self):
        distances = np.linalg.norm(self.raw_coordinates, axis=1)     # Calculate Euclidean distances
        closest_index = np.argmin(distances)                         # Find the index of the closest point
        closest_point = self.raw_coordinates[closest_index]          # Get the closest point coordinates
        return closest_point

    def format_airfoil(self,closest_point):
        #normalize airfoil to unit chord
        self.data[0] -= np.min(self.data[0])                         #This step shifts the x-coordinates such that the minimum x-coordinate becomes 0.
        divisor = np.max(self.data[0])                               #he purpose of this divisor is to normalize the x-coordinates such that the maximum x-coordinate becomes 1.
        self.data[0] /= divisor                                      #This step scales the x-coordinates to the range [0, 1].
        self.data[1] /= divisor                                      #This step scales the y-coordinates to maintain their relative proportions to the x-coordinates.

        #getting the leading and trailing edge
        if self.camber_builder:
            try:
                max_x_index = np.nanargmax(self.raw_coordinates[:, 0])
                self.trailing_edge = tuple(self.raw_coordinates[max_x_index])
                print(self.trailing_edge)



                min_x_index = np.nanargmin(self.raw_coordinates[:, 0])
                self.leading_edge = tuple(closest_point)
                print(min_x_index,self.leading_edge)

                self.chord_length = math.dist(self.leading_edge, self.trailing_edge)
            except TypeError as e:
                print('There is no camberline to get the LE')


    def rotate_airfoil(self,angle=0):
        '''
        A function to rotate the gridpoints of the object to create an AoA
        :param angle: angle to rotate the airfoil grid [deg]
        :return: grid of the airfoil updated to new AoA
        '''
        angle = -angle * np.pi / 180

        def rotateMatrix(angle):
            return np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        if self.camber_builder:
            try:
                self.camber_line_transposed = self.camber_line.transpose() @ rotateMatrix(angle).T
                self.camber_line = self.camber_line_transposed.transpose()
                self.trailing_edge = self.trailing_edge @ rotateMatrix(angle).T
            except:
                print("this dataset has no camberline defined")

        self.raw_coordinates = self.raw_coordinates @ rotateMatrix(angle).T

        self.data =self.raw_coordinates.transpose()



    def plot_airfoil(self):
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(self.data[0],self.data[1],s=3, color='black')
        if self.camber_builder:
            ax.scatter(self.trailing_edge[0], self.trailing_edge[1], color='blue')
            ax.scatter(self.camber_line[0], self.camber_line[1],s=2,color='0.8')
            ax.scatter(self.leading_edge[0], self.leading_edge[1], color='red')





        # setting axis and title
        ax.axis([-0.25, 1.25, -1, 1])
        ax.set_title(self.name, fontsize=14)
        ax.set_xlabel('x-pos', fontsize=10)
        ax.set_ylabel('y-pos', fontsize=10)

        plt.show() #doesnt work in this backend config

        #plt.savefig(f'airfoil point representation{self.name}1.PNG', bbox_inches='tight', dpi=300)



















'''
    def get_airfoil_Javafoil(self):
        filepath = glob.glob('FoilToAnalize\\*')[0]
        with open(filepath, "r+") as f:
            file = f.readlines()
            fileData = file[1:]
            name = file[:1]
            self.data = np.loadtxt(fileData, delimiter='\t', unpack=True)

            self.name = ''.join(name)[0:-1:1]
        self.x = self.data[0]
        self.y = self.data[1]
        self.raw_coordinates = self.data.transpose()  # np.array((x, y)) #array of shape (2,n): two down and n right
'''