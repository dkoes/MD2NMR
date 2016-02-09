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
oxygen_distances = []
nitrog_distances = []

# Set print options to only display 3 decimal places
np.set_printoptions(threshold=np.inf,precision=3)

# Computes distances function sqrt( (x2 - x1)^2 + (y2-y1)^2 + (z2-z1)^2 )
def distance(atom1, atom2):
	x = atom1[0] - atom2[0]
	y = atom1[1] - atom2[1]
	z = atom1[2] - atom2[2]
	return math.sqrt((x * x) + (y * y) + (z * z))

'''# Open all files at the beginning
for i in range(2, len(res)-2):
	file_name_O = 'DUMP.O.' + str(i+2) + '.txt'
	f_O = open(file_name_O,'w')
	open_O_files.append(f_O)
	file_name_N = 'DUMP.N.' + str(i+2) + '.txt'
	f_N = open(file_name_N,'w')
	open_N_files.append(f_N)'''

next_res_dist = []

# Iterates through every residue in each frame, computing the distance between N's and O's
for ts in u.trajectory:
	
	for i in range(2, len(res)-2):
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

			oxy_res_distances = [] # O distances for THIS residue
			nit_res_distances = [] # N distances for THIS residue

			oxy_res_distances.append(d0)	# d0
			oxy_res_distances.append(str("{0:.3f}".format(distance(res[i-1].O.pos, res[i].O.pos))))	# d1
			oxy_res_distances.append(d2)	# d2
			oxy_res_distances.append(d3)	# d3
			oxy_res_distances.append(str("{0:.3f}".format(distance(res[i+2].N.pos, res[i].O.pos))))	# d4
			oxygen_distances.append(oxy_res_distances)

			nit_res_distances.append(d0)	# d5
			nit_res_distances.append(str("{0:.3f}".format(distance(res[i-1].N.pos, res[i].N.pos))))	# d6
			nit_res_distances.append(str("{0:.3f}".format(distance(res[i-2].O.pos, res[i].N.pos))))	# d7
			nit_res_distances.append(str("{0:.3f}".format(distance(res[i-1].O.pos, res[i].N.pos))))	# d8
			nit_res_distances.append(d9)	# d9
			nitrog_distances.append(nit_res_distances)
		else:
			# New calculations
			d0 = str("{0:.3f}".format(distance(res[ i ].N.pos, res[i].O.pos)))
			d2 = str("{0:.3f}".format(distance(res[i+1].N.pos, res[i].O.pos))) # Next res d8
			d3 = str("{0:.3f}".format(distance(res[i+1].O.pos, res[i].O.pos))) # Next res d1
			d9 = str("{0:.3f}".format(distance(res[i+1].N.pos, res[i].N.pos))) # Next res d6 

			oxy_res_distances = [] # O distances for THIS residue
			nit_res_distances = [] # N distances for THIS residue

			oxy_res_distances.append(d0)	# d0
			oxy_res_distances.append(next_res_dist[1])	# d1
			oxy_res_distances.append(d2)	# d2
			oxy_res_distances.append(d3)	# d3
			oxy_res_distances.append(str("{0:.3f}".format(distance(res[i+2].N.pos, res[i].O.pos))))	# d4
			oxygen_distances.append(oxy_res_distances)

			nit_res_distances.append(d0)	# d5
			nit_res_distances.append(next_res_dist[2])	# d6
			nit_res_distances.append(str("{0:.3f}".format(distance(res[i-2].O.pos, res[i].N.pos))))	# d7
			nit_res_distances.append(next_res_dist[0])	# d8
			nit_res_distances.append(d9)	# d9
			nitrog_distances.append(nit_res_distances)

			next_res_dist[0] = d2
			next_res_dist[1] = d3
			next_res_dist[2] = d9

result = ''
all_ox_dist = []
all_ni_dist = []

for j in range (0, len(range(2,len(res)-2))):
	frame_num = 0
	result = ''
	dist_range = range(j, len(oxygen_distances), len(range(2,len(res)-2)))
	res_range = range(2, len(res)-2)	
	for i in dist_range:
		result += 'O|' + str(frame_num)
		frame_num += 1
		for dist in oxygen_distances[i]:
			result += '|' + dist
		result += '\n'
	all_ox_dist.append(result)

for i in range(0, len(all_ox_dist)):
	file_name = 'DUMP.O.' + str(i+4) + '.txt'
	with open(file_name, 'w') as f:
		f.write(all_ox_dist[i])

for j in range (0, len(range(2,len(res)-2))):
	frame_num = 0
	result = ''
	dist_range = range(j, len(nitrog_distances), len(range(2,len(res)-2)))
	res_range = range(2, len(res)-2)	
	for i in dist_range:
		result += 'N|' + str(frame_num)
		frame_num += 1
		for dist in nitrog_distances[i]:
			result += '|' + dist
		result += '\n'
	all_ni_dist.append(result)

for i in range(0, len(all_ni_dist)):
	file_name = 'DUMP.N.' + str(i+4) + '.txt'
	with open(file_name, 'w') as f:
		f.write(all_ni_dist[i])

print "Program time: " + str("{0:.3f}".format(time.time() - start_time)) + " seconds."
