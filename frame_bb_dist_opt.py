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

next_res_dist = []
open_O_files = []
open_N_files = []

# Open all files at the beginning
for i in range(2, len(res)-2):
	file_name_O = 'DUMP.O.' + str(i+2) + '.txt'
	f_O = open(file_name_O,'w')
	open_O_files.append(f_O)
	file_name_N = 'DUMP.N.' + str(i+2) + '.txt'
	f_N = open(file_name_N,'w')
	open_N_files.append(f_N)

# Iterates through every residue in each frame, computing the distance between N's and O's
for ts in u.trajectory:
	for i, f in zip(range(2, len(res)-2), range(0, len(open_O_files))):
		if i == 2:
			# Calculations that will be used more than once
			d0 = str("{0:.3f}".format(distance(res[ i ].N.pos, res[i].O.pos))) # Same as d5
			d2 = str("{0:.3f}".format(distance(res[i+1].N.pos, res[i].O.pos))) # Next res d8
			d3 = str("{0:.3f}".format(distance(res[i+1].O.pos, res[i].O.pos))) # Next res d1
			d9 = str("{0:.3f}".format(distance(res[i+1].N.pos, res[i].N.pos))) # Next res d6
			
			# Storing distances for use in next residue
			next_res_dist.append(d2) # d8
			next_res_dist.append(d3) # d1
			next_res_dist.append(d9) # d6

			# Compute calculations as distances are written to files
			open_O_files[f].write('O|' + str(ts.frame) + '|' + d0 + '|' + \
						 str("{0:.3f}".format(distance(res[i-1].O.pos, res[i].O.pos))) + '|' + \
						 d2 + '|' + d3 + '|' + str("{0:.3f}".format(distance(res[i+2].N.pos, res[i].O.pos))) + '\n')
		
			open_N_files[f].write('N|' + str(ts.frame) + '|' + d0 + '|' + \
						 str("{0:.3f}".format(distance(res[i-1].N.pos, res[i].N.pos))) + '|' + \
						 str("{0:.3f}".format(distance(res[i-2].O.pos, res[i].N.pos))) + '|' + \
						 str("{0:.3f}".format(distance(res[i-1].O.pos, res[i].N.pos))) + '|' + d9 + '\n')
		else:
			# New calculations
			d0 = str("{0:.3f}".format(distance(res[ i ].N.pos, res[i].O.pos)))
			d2 = str("{0:.3f}".format(distance(res[i+1].N.pos, res[i].O.pos))) # Next res d8
			d3 = str("{0:.3f}".format(distance(res[i+1].O.pos, res[i].O.pos))) # Next res d1
			d9 = str("{0:.3f}".format(distance(res[i+1].N.pos, res[i].N.pos))) # Next res d6 

			# Compute calculations as distances are written to files
			open_O_files[f].write('O|' + str(ts.frame) + '|' + d0 + '|' + next_res_dist[1] + '|' + \
						 d2 + '|' + d3 + '|' + \
						 str("{0:.3f}".format(distance(res[i+2].N.pos, res[i].O.pos))) + '\n')
		
			open_N_files[f].write('N|' + str(ts.frame) + '|' + d0 + '|' + next_res_dist[2] + '|' + \
						 str("{0:.3f}".format(distance(res[i-2].O.pos, res[i].N.pos))) + '|' + \
						 next_res_dist[0] + '|' + d9 + '\n')

			next_res_dist[0] = d2
			next_res_dist[1] = d3
			next_res_dist[2] = d9

print "Program time: " + str("{0:.3f}".format(time.time() - start_time)) + " seconds."
