#import the required libraries
import numpy as np
import os
from pathlib import Path

# go through all files in the gridfiles folder
# So, no unnecessary file should exist in this folder
# change the directory based on where you care running the script
directory = "./run/Tutorials/2DBump/CreateBlocks/grid"
gridfiles = [os.path.join(root,file) for root,dir,f in os.walk(directory) for file in f]

# sort the name of the gridfile 
gridfiles.sort()

# loop through the gridfiles and process the VOf expression 
# and write it in a separate file
for index, gridfile in enumerate(gridfiles):
    # Info to output
    #print("Working on file no. {0:d}: {1}".format(index, gridfile))
    # open the file to read
    with open(gridfile, 'r') as f:
        # load header of the file
        nx, ny, nz = [int(dim) for dim in f.readline().split()]
        #load the grid points of the file
        x,y,z = np.loadtxt(f, unpack=True)
        
    # Now we will evaluate VOF expression for 
    vof = np.zeros_like(x)
    vof[x<0] = 1
    # end of vof expression

    # write the output in a file
    ## Directory where the vof files will be stored
    vofDirectory="./run/Tutorials/2DBump/CreateBlocks/vof"
    # create the vof directory if does not exist
    Path(vofDirectory).mkdir(exist_ok=True)
    # name of the vof file to write in vofdirectory 
    # corresponding to (Index) block number
    filename = "{0}/block_{1:02d}.vof".format(vofDirectory, index)
    # header of the file which contains information about 
    # number of points in i,j,k direction
    header = "{0:d} {1:} {2:d}".format(nx, ny, nz)
    np.savetxt(filename, vof, header=header)
    #Info
    print("Written file no. {0:d}: {1}".format(index, filename))