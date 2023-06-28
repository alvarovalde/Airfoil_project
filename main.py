import time
from airfoil import Airfoil_manager
from directory_management import clean_directory
from grid_creator import Naca4Creator
import threading
from PySide6 import QtGui, QtCore, QtWidgets


def main():
    # airfoil creator
    t1_start = time.perf_counter()
    x = Naca4Creator(NACA=2112, n_points=500)
    x.make_airfoil()

    # airfoil asignation
    foil = Airfoil_manager()
    foil.get_airfoil_metrics()
    foil.format_airfoil(foil.closest_point_to_origin())
    foil.rotate_airfoil(30)
    t1_stop = time.perf_counter()
    foil.plot_airfoil()
    print("Elapsed time during the whole program in seconds:",
          t1_stop - t1_start)


if __name__ == '__main__':
    clean_directory()
    main()
