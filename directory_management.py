import glob
import os
def clean_directory():
    #clean directory before
        filepath = glob.glob('FoilToAnalize\\*')
        for f in filepath:
            os.remove(f)