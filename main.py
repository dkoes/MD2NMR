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

# Holds calculations that are used more than once
next_res_dist = []

# Holds all open files
open_O_files = []
open_N_files = []

# Open all files at the beginning
for i in range(2, len(res)-2):
	file_name_O = 'DUMP.O.' + str(i+1) + '.txt'
	file_name_N = 'DUMP.N.' + str(i+1) + '.txt'
	
	f_O = open(file_name_O,'w')
	f_N = open(file_name_N,'w')
	
	open_O_files.append(f_O)
	open_N_files.append(f_N)

# Iterates through every residue in each frame, computing the distance between N's and O's
for ts in u.trajectory:
	for i, f in zip(range(2, len(res)-2), range(0, len(open_O_files))):
		if i == 2 and res[i].name not in "PRO":
			# d3 becomes d1 in the next residue
			d3 = str("{0:.3f}".format(distance(res[i+1].O.pos, res[i].O.pos)))
			next_res_dist.append(d3) 

			# Compute calculations as distances are written to files
			open_O_files[f].write('O|' + str(ts.frame) + '|' + \
						 str("{0:.3f}".format(distance(res[ i ].N.pos, res[i].O.pos))) + '|' + \
						 str("{0:.3f}".format(distance(res[i-1].O.pos, res[i].O.pos))) + '|' + \
						 str("{0:.3f}".format(distance(res[i+1].N.pos, res[i].O.pos))) + '|' + \
						 d3 + '|' + str("{0:.3f}".format(distance(res[i+2].N.pos, res[i].O.pos))) + '\n')
		
			open_N_files[f].write('N|' + str(ts.frame) + '|' + \
						 str("{0:.3f}".format(distance(res[i-1].N.pos, res[i].H.pos))) + '|' + \
						 str("{0:.3f}".format(distance(res[i-2].O.pos, res[i].H.pos))) + '|' + \
						 str("{0:.3f}".format(distance(res[i-1].O.pos, res[i].H.pos))) + '|' + \
						 str("{0:.3f}".format(distance(res[ i ].O.pos, res[i].H.pos))) + '|' + \
						 str("{0:.3f}".format(distance(res[i+1].N.pos, res[i].H.pos))) + '\n')
		elif i > 2 and res[i].name not in "PRO":
			# d3 becomes d1 in the next residue
			d3 = str("{0:.3f}".format(distance(res[i+1].O.pos, res[i].O.pos)))

			# Compute calculations as distances are written to files
			open_O_files[f].write('O|' + str(ts.frame) + '|' + \
						 str("{0:.3f}".format(distance(res[ i ].N.pos, res[i].O.pos))) + '|' + next_res_dist[0] + '|' + \
						 str("{0:.3f}".format(distance(res[i+1].N.pos, res[i].O.pos))) + '|' + d3 + '|' + \
						 str("{0:.3f}".format(distance(res[i+2].N.pos, res[i].O.pos))) + '\n')
		
			open_N_files[f].write('N|' + str(ts.frame) + '|' + \
						 str("{0:.3f}".format(distance(res[i-1].N.pos, res[i].H.pos))) + '|' + \
						 str("{0:.3f}".format(distance(res[i-2].O.pos, res[i].H.pos))) + '|' + \
						 str("{0:.3f}".format(distance(res[i-1].O.pos, res[i].H.pos))) + '|' + \
						 str("{0:.3f}".format(distance(res[ i ].O.pos, res[i].H.pos))) + '|' + 
						 str("{0:.3f}".format(distance(res[i+1].N.pos, res[i].H.pos))) + '\n')

			next_res_dist[0] = d3
		else:
			continue

print "Program time: " + str("{0:.3f}".format(time.time() - start_time)) + " seconds."
