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

# Iterates through every residue in each frame, computing the distance between N's and O's
for ts in u.trajectory:
	
	for i in range(2, len(res)-2):
		'''print '\n'
		print '========================    Frame # is ' + str(ts.frame) + '    ========================'  + '\n'
		print 'Residue #: \t' + str(res[i].id)
		print 'i = \t' + str(i)
		print '\n' '''

		oxy_res_distances = [] # O distances for THIS residue
		nit_res_distances = [] # N distances for THIS residue
		
		oxy_res_distances.append(distance(res[ i ].N.pos, res[i].O.pos))	# d0
		oxy_res_distances.append(distance(res[i-1].O.pos, res[i].O.pos))	# d1
		oxy_res_distances.append(distance(res[i+1].N.pos, res[i].O.pos))	# d2
		oxy_res_distances.append(distance(res[i+1].O.pos, res[i].O.pos))	# d3
		oxy_res_distances.append(distance(res[i+2].N.pos, res[i].O.pos))	# d4
		oxygen_distances.append(oxy_res_distances)

		nit_res_distances.append(distance(res[ i ].O.pos, res[i].N.pos))	# d5
		nit_res_distances.append(distance(res[i-1].N.pos, res[i].N.pos))	# d6
		nit_res_distances.append(distance(res[i-2].O.pos, res[i].N.pos))	# d7
		nit_res_distances.append(distance(res[i-1].O.pos, res[i].N.pos))	# d8
		nit_res_distances.append(distance(res[i+1].N.pos, res[i].N.pos))	# d9
		nitrog_distances.append(nit_res_distances)
		
		'''print 'res[i].N.pos:\t' + str(res[i].N)
		print 'res[i-1].O.pos\t' + str(res[i-1].O)
		print 'res[i+1].N.pos\t' + str(res[i+1].N)
		print 'res[i+1].O.pos\t' + str(res[i+1].O)
		print 'res[i+2].N.pos\t' + str(res[i+2].N)
		print '\n'
		print oxy_res_distances
		print '\n\n'
		
		print 'res[i].O.pos\t' + str(res[i].O)
		print 'res[i-1].N.pos\t' + str(res[i-1].N)
		print 'res[i-2].O.pos\t' + str(res[i-2].O)
		print 'res[i-1].O.pos\t' + str(res[i-1].O)
		print 'res[i+1].N.pos\t' + str(res[i+1].N)
		print '\n'
		print nit_res_distances
		print '\n\n' '''

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
			result += '|' + str("{0:.3f}".format(dist))
		result += '\n'
	all_ox_dist.append(result)

for i in range(-1, len(all_ox_dist)):
	file_name = 'DUMP.O.' + str(i+4) + '.txt'
	with open(file_name, 'w') as f:
		f.write(all_ox_dist[i])

for j in range (0, len(range(2,len(res)-2))):
	frame_num = 0
	result = ''
	dist_range = range(j, len(nitrog_distances), len(range(2,len(res)-2)))
	res_range = range(2, len(res)-2)	
	for i in dist_range:
		result += 'O|' + str(frame_num)
		frame_num += 1
		for dist in nitrog_distances[i]:
			result += '|' + str("{0:.3f}".format(dist))
		result += '\n'
	all_ni_dist.append(result)

for i in range(-1, len(all_ni_dist)):
	file_name = 'DUMP.N.' + str(i+4) + '.txt'
	with open(file_name, 'w') as f:
		f.write(all_ni_dist[i])

print "Program time: " + str("{0:.3f}".format(time.time() - start_time)) + " seconds."
