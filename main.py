import time
from airfoil import Airfoil_manager
from directory_management import clean_directory
from Naca4 import Naca4Creator
import threading
from PySide6 import QtGui, QtCore, QtWidgets


def main():
    # airfoil creator
    t1_start = time.perf_counter()
    x = Naca4Creator(NACA=2112, n_points=40,export_camberline=True)
    x.generate_airfoil()


    # airfoil asignation
    foil = Airfoil_manager()
    foil.manager(5)

    pm = Panel_method()
    pm.get_airfoil_metrics()
    pm.plot()
    pm.define_panels()
    pm.plot_panels()


    t1_stop = time.perf_counter()

    print("Elapsed time during the whole program in seconds:",
          t1_stop - t1_start)


if __name__ == '__main__':
    clean_directory()
    main()
