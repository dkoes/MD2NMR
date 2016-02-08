#!/usr/bin/python

# Script that analyzes a simulation and outputs the 
# frame number and backbone distances.

# Print format
# <frame_number>|<backbone distances>|
import sys, MDAnalysis
import math
import numpy as np
import time
from os.path import splitext
from MDAnalysis.analysis.rms import *
from MDAnalysis.analysis.distances import *

print '========================================================================================'
print '=============== HELLO ================================= FRAME_BB_DIST.PY ==============='
print '========================================================================================'

start_time=time.time()

# Set up file names
top_name = sys.argv[1]
traj_name = sys.argv[2]

# Set up universe, selections and data structures
u = MDAnalysis.Universe(top_name, traj_name)
prot = u.select_atoms("protein")
res = prot.atoms.residues

# Set print options to only display 3 decimal places
np.set_printoptions(threshold=np.inf,precision=3)

# Computes distances function sqrt( (x2 - x1)^2 + (y2-y1)^2 + (z2-z1)^2 )
def distance(atom1, atom2):
	x = atom1[0] - atom2[0]
	y = atom1[1] - atom2[1]
	z = atom1[2] - atom2[2]
	return math.sqrt((x * x) + (y * y) + (z * z))

# Iterates through every residue in each frame, computing the distance between N's and O's
for ts in u.trajectory:
	for i in range(2, len(res)-2):
		#Print into the open files
		file_name_O = 'DUMP.O.' + str(i) + '.txt'
		file_name_N = 'DUMP.N.' + str(i) + '.txt'

		with open(file_name_O, "a") as file_O:
			file_O.write('O|' + str(ts.frame) + '|' + \
					   str("{0:.3f}".format(distance(res[ i ].N.pos, res[i].O.pos))) + '|' + \
					   str("{0:.3f}".format(distance(res[i-1].O.pos, res[i].O.pos))) + '|' + \
					   str("{0:.3f}".format(distance(res[i+1].N.pos, res[i].O.pos))) + '|' + \
					   str("{0:.3f}".format(distance(res[i+1].O.pos, res[i].O.pos))) + '|' + \
					   str("{0:.3f}".format(distance(res[i+2].N.pos, res[i].O.pos))) + '\n')
		
		
		with open(file_name_N, "a") as file_N:
			file_N.write('N|' + str(ts.frame) + '|' + \
					   str("{0:.3f}".format(distance(res[ i ].O.pos, res[i].N.pos))) + '|' + \
					   str("{0:.3f}".format(distance(res[i-1].N.pos, res[i].N.pos))) + '|' + \
					   str("{0:.3f}".format(distance(res[i-2].O.pos, res[i].N.pos))) + '|' + \
					   str("{0:.3f}".format(distance(res[i-1].O.pos, res[i].N.pos))) + '|' + \
					   str("{0:.3f}".format(distance(res[i+1].N.pos, res[i].N.pos))) + '\n')

print "Program time: " + str("{0:.3f}".format(time.time() - start_time)) + " seconds."
