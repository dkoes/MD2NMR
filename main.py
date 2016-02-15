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

# Dictionary that maps atom names to groups
atom_reference = {
"NALA": "N", "NARG": "N", "NASN": "N", "NASP": "N", "NCYS": "N", "NGLN": "N", "NGLY": "N", "NGLU": "N",
"NHSD": "N", "NHSE": "N", "NHSP": "N", "NILE": "N", "NLEU": "N", "NLYS": "N", "NMET": "N", "NPHE": "N",
"NPRO": "N", "NSER": "N", "NTHR": "N", "NTRP": "N", "NTYR": "N", "NVAL": "N",
"OALA": "O", "OARG": "O", "OASN": "O", "OASP": "O", "OCYS": "O", "OGLN": "O", "OGLY": "O", "OGLU": "O",
"OHSD": "O", "OHSE": "O", "OHSP": "O", "OILE": "O", "OLEU": "O", "OLYS": "O", "OMET": "O", "OPHE": "O",
"OPRO": "O", "OSER": "O", "OTHR": "O", "OTRP": "O", "OTYR": "O", "OVAL": "O", "OD1ASN": "O", "OE1GLN": "O",
"OH2TIP3": "W", "CD1PHE": "R", "CD1TRP": "R", "CD1TYR": "R", "CD2HSD": "R", "CD2HSE": "R", "CD2HSP": "R",
"CD2PHE": "R", "CD2TRP": "R", "CD2TYR": "R", "CE1HSD": "R", "CE1HSE": "R", "CE1HSP": "R",
"CE1PHE": "R", "CE1TYR": "R", "CE2PHE": "R", "CE2TRP": "R", "CE2TYR": "R", "CE3TRP": "R",
"CGHSD": "R", "CGHSE": "R", "CGHSP": "R", "CGPHE": "R", "CGTRP": "R", "CGTYR": "R",
"CH2TRP": "R", "CZ2TRP": "R", "CZ3TRP": "R", "CZPHE": "R", "CZTYR": "R", "ND1HSD": "R",
"ND1HSE": "R", "ND1HSP": "R", "NE1TRP": "R", "NE2HSD": "R", "NE2HSE": "R", "NE2HSP": "R",	
"OG1THR": "L", "OGSER":  "L", "OHTYR":  "L", "OD2ASPH": "L", "OE2GLUH": "L", "ND2ASN": "D", "NE2GLN": "D", 
"SDMET":  "S", "SGCYS":  "S", "OD1ASP": "A", "OD2ASP": "A", "OE1GLU": "A", "OE2GLU": "A", "OT1ALA": "A", 
"OT1ARG": "A", "OT1ASN": "A", "OT1ASP": "A", "OT1CYS": "A", "OT1GLN": "A", "OT1GLY": "A", "OT1GLU": "A", 
"OT1HSD": "A", "OT1HSE": "A", "OT1HSP": "A", "OT1ILE": "A", "OT1LEU": "A",  "OT1LYS": "A", "OT1MET": "A", 
"OT1PHE": "A", "OT1PRO": "A", "OT1SER": "A", "OT1THR": "A", "OT1TRP": "A", "OT1TYR": "A", "OT1VAL": "A",
"NZLYS": "B", "NTALA": "B", "NEARG": "G", "NH1ARG": "G", "NH2ARG": "G", "SODSOD": "P", "CLACLA": "M",
"Z": "Z", "OD1ASPH": "U", "OE1GLUH": "U", "CALA": "C", "CARG": "C", "CASN": "C", "CASP": "C", "CCYS": "C", 
"CGLN": "C", "CGLY": "C", "CGLU": "C", "CHSD": "C", "CHSE": "C", "CHSP": "C", "CILE": "C", "CLEU": "C", 
"CLYS": "C", "CMET": "C", "CPHE": "C", "CPRO": "C", "CSER": "C", "CTHR": "C", "CTRP": "C","CTYR": "C", "CVAL": "C"
}

#for key, value in atom_reference.iteritems():
	#print key, value

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

# Iterates through every residue in each frame, computing and writing the distance between N's and O's
for ts in u.trajectory:
	for i, f in zip(range(2, len(res)-2), range(0, len(open_O_files))):
		if i == 2 and res[i].name not in "PRO":
			#Grab all atoms within 5.0 A of target
			print "== FRAME NUMBER: " + str(ts.frame) + " =="
			print str(res[i].H.number + 1), res[i].H
			cloud = u.select_atoms("around 5.0 atom SYSTEM " + str(res[i].id) + " H")
			print "Length of cloud atoms: \t" + str(len(cloud))
			for at in cloud:
				print "Atom name: \t" + at.name
				if at.name in atom_reference.keys():
					print "Atom name found in dictionary:"
					print atom_reference[at.name]
					
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
