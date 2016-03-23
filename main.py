#!/usr/bin/python

# Script that analyzes a simulation and outputs the 
# frame number and atom distances.

import sys, MDAnalysis
import math
import numpy as np
import time
from os.path import splitext
from MDAnalysis.analysis.rms import *
from MDAnalysis.analysis.distances import *
import __builtin__

# Retrieve input files
exc = getattr(__builtin__,"IOError","FileNotFoundError")
	
while True:
	top_name = raw_input("Hello, please type in the name of the topology file (.prmtop) and press 'Enter': ")
	try:
		top = open(top_name)
	except exc:
		print("Wrong file or file path")
	else:
		break

while True:
	traj_name = raw_input("Hello, please type in the name of the trajectory file (.dcd) and press 'Enter': ")
	try:
		traj = open(traj_name)
	except exc:
		print("Wrong file or file path")
	else:
		break

start_time=time.time()

# Set up universe, selections and data structures
u = MDAnalysis.Universe(top_name, traj_name)
prot = u.select_atoms("protein")

res = prot.atoms.residues
bb_atoms = ["N", "CA", "C", "O"]

# Dictionary that maps atom names to groups
atom_reference = {
"NALA": "N", "NARG": "N", "NASN": "N", "NASP": "N", "NCYS": "N", "NGLN": "N", "NGLY": "N", "NGLU": "N",
"NHSD": "N", "NHSE": "N", "NHSP": "N", "NILE": "N", "NLEU": "N", "NLYS": "N", "NMET": "N", "NPHE": "N",
"NPRO": "N", "NSER": "N", "NTHR": "N", "NTRP": "N", "NTYR": "N", "NVAL": "N",

"OALA": "O", "OARG": "O", "OASN": "O", "OASP": "O", "OCYS": "O", "OGLN": "O", "OGLY":   "O", "OGLU":   "O",
"OHSD": "O", "OHSE": "O", "OHSP": "O", "OILE": "O", "OLEU": "O", "OLYS": "O", "OMET":   "O", "OPHE":   "O",
"OPRO": "O", "OSER": "O", "OTHR": "O", "OTRP": "O", "OTYR": "O", "OVAL": "O", "OD1ASN": "O", "OE1GLN": "O",

"OH2TIP3": "W", 

"CD1PHE": "R", "CD1TRP": "R", "CD1TYR": "R", "CD2HSD": "R", "CD2HSE": "R", "CD2HSP": "R",
"CD2PHE": "R", "CD2TRP": "R", "CD2TYR": "R", "CE1HSD": "R", "CE1HSE": "R", "CE1HSP": "R",
"CE1PHE": "R", "CE1TYR": "R", "CE2PHE": "R", "CE2TRP": "R", "CE2TYR": "R", "CE3TRP": "R",
"CGHSD":  "R", "CGHSE":  "R", "CGHSP":  "R", "CGPHE":  "R", "CGTRP":  "R", "CGTYR":  "R",
"CH2TRP": "R", "CZ2TRP": "R", "CZ3TRP": "R", "CZPHE":  "R", "CZTYR":  "R", "ND1HSD": "R",
"ND1HSE": "R", "ND1HSP": "R", "NE1TRP": "R", "NE2HSD": "R", "NE2HSE": "R", "NE2HSP": "R",

"OG1THR": "L", "OGSER":  "L", "OHTYR":  "L", "OD2ASPH": "L", "OE2GLUH": "L", 

"ND2ASN": "D", "NE2GLN": "D",

"SDMET":  "S", "SGCYS":  "S", 

"OD1ASP": "A", "OD2ASP": "A", "OE1GLU": "A", "OE2GLU": "A", "OT1ALA": "A", "OT1ARG": "A", "OT1ASN": "A", 
"OT1ASP": "A", "OT1CYS": "A", "OT1GLN": "A", "OT1GLY": "A", "OT1GLU": "A", "OT1HSD": "A", "OT1HSE": "A", 
"OT1HSP": "A", "OT1ILE": "A", "OT1LEU": "A", "OT1LYS": "A", "OT1MET": "A", "OT1PHE": "A", "OT1PRO": "A", 
"OT1SER": "A", "OT1THR": "A", "OT1TRP": "A", "OT1TYR": "A", "OT1VAL": "A",

"NZLYS": "B", "NTALA": "B", 

"NEARG": "G", "NH1ARG": "G", "NH2ARG": "G", 

"SODSOD": "P", 

"CLACLA": "M",

"Z": "Z", 

"OD1ASPH": "U", "OE1GLUH": "U", 

"CALA": "C", "CARG": "C", "CASN": "C", "CASP": "C", "CCYS": "C", "CGLN": "C", "CGLY": "C", "CGLU": "C", 
"CHSD": "C", "CHSE": "C", "CHSP": "C", "CILE": "C", "CLEU": "C", "CLYS": "C", "CMET": "C", "CPHE": "C", 
"CPRO": "C", "CSER": "C", "CTHR": "C", "CTRP": "C", "CTYR": "C", "CVAL": "C"
}

# Dictionary that maps atom type to ranking preference
atom_ranking = {
	"N": 0, "O": 1, "W": 2, "R":  3, "L":  4, "D":  5, "S":  6, "A":  7,
	"B": 8, "G": 9, "P":10, "M": 11, "Z": 12, "U": 13, "C": 14
}

######################################################################################################################
## Helper functions - calculate distance function, find closest atoms, etc.                                         ##
######################################################################################################################

# Computes distances function sqrt( (x2 - x1)^2 + (y2-y1)^2 + (z2-z1)^2 )
def distance(atom1, atom2):
	x = atom1[0] - atom2[0]
	y = atom1[1] - atom2[1]
	z = atom1[2] - atom2[2]
	return math.sqrt((x * x) + (y * y) + (z * z))

# Calculates the angle between the two vectors x and v
def vangle(x, v):
    theta = 0.0
    
    dot = x[0] * v[0] + x[1] * v[1] + x[2] * v[2]
    
    udx = math.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])
    vdx = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])

    udxvdx = dot/(udx*vdx);

    if(udxvdx > 1.00000000):
    	udxvdx = 1.00000000

    theta = math.degrees(math.acos(udxvdx))
    
    if(theta > 180.000):
    	theta = 180.000

    return theta


# Determines nearest carbon to current target residue atom, based upon its location of oxygen
def find_closest_atom(cloud, target_atom):
	minC = 10000.0
	minO = 10000.0
	
	# Organize carbon and oxygen atom types
	carbons = []
	oxygens = []
	other_atoms = []
	nearestC = None
	nearestO = None
	
	ref_name = ''
	nearest_atoms = []
	
	# Classifies carbon and oxygen atom types
	for atom in cloud:
		ref_name = atom[0].name + atom[0].resname
		if atom_reference[ref_name] == 'C':
			carbons.append(atom)
		elif atom_reference[ref_name] == 'O':
			oxygens.append(atom)
		else:
			other_atoms.append(atom)

	# Finds closest oxygen to current residue target atom
	for atom in oxygens:
		dist = distance(atom[0].pos, target_atom.pos)
		if dist < minO:
			minC = dist
			nearestO = atom
	
	# If there are no O atoms or any other atoms, return no atoms (None)
	if nearestO == None and len(other_atoms) == 0:
		return None
	# If there are no O atoms and there are other atoms, return the other atoms
	elif nearestO == None and len(other_atoms) > 0:
		return other_atoms
	# If there is an O, add it first to the list of closest atoms, then find closest C to nearestO
	else:
		if nearestO:
			nearest_atoms.append(nearestO)
	
	# Finds closest C to previously found nearestO
	if len(carbons) > 0:
		for atom in carbons:
			dist = distance(atom[0].pos, nearestO[0].pos)
			if dist < minC:
				minC = dist
				nearestC = atom
	
	if nearestC:
		nearest_atoms.append(nearestC)

	if len(other_atoms) > 0:
		for atom in other_atoms:
			nearest_atoms.append(atom)
	
	return nearest_atoms


######################################################################################################################
## Main code - opens all files, iterates through each timestep once while iterating through residues.               ##
######################################################################################################################

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
	print "\nProcessing frame #: " + str(ts.frame)
	for i, f in zip(range(2, len(res)-2), range(0, len(open_O_files))):
		if i == 2 and res[i].name not in "PRO":
			######################################################
			##  Process the N dome								##
			######################################################
			
			# Select all atoms within 5.0 A of target
			cloud = u.select_atoms("around 5.0 atom SYSTEM " + str(res[i].id) + " H")
			processed_cloud = []
			processed_cloud_02 = []
			atom_pattern_N = ''
			atom_pattern_dist = ''
			
			# First round of processing
			for atom in cloud:
				ref_name = atom.name + atom.resname
				
				# Check to ensure backbone atoms within (+/-) 2 residues are not included
				if atom.name in bb_atoms and abs(res[i].id - atom.resid) < 3:
					continue
				else:
					if ref_name in atom_reference:
						x = res[i].N.pos - res[i].H.pos
						v = atom.pos - res[i].H.pos

						# Eliminate any atoms that have vector angles less than 90.0
						if abs(vangle(x, v)) < 90.0:
							continue
						else:
							######################################################
							##  Compute distance between atoms based on types	##
							##  Aromatic atoms: within 5.0 A 					##
							##  Carbon atoms: 	within 4.0 A 					##
							##  Other atoms: 	wtihin 2.5 A 					##
							######################################################
							if atom_reference[ref_name] in "R":
								processed_cloud.append([atom, distance(atom.pos, res[i].H.pos)])
							elif atom_reference[ref_name] in "C":
								C_dist = distance(atom.pos, res[i].H.pos)
								if C_dist <= 4.0:
									processed_cloud.append([atom, C_dist])
							else:
								dist = distance(atom.pos, res[i].H.pos)
								if dist <= 2.5:
									processed_cloud.append([atom, dist])
			
			# Second round of processing					
			if processed_cloud is None:
				atom_pattern_N = "Z:"
			else:
				processed_cloud_02 = find_closest_atom(processed_cloud, res[i].H)
			
				if processed_cloud_02 == None:
					atom_pattern_N = "Z:"
				else:
					# Rank atoms based on the atom type preferences noted above
					processed_cloud_02.sort(key = lambda atom: atom_ranking[atom_reference[atom[0].name + atom[0].resname]])
					
					# Generate the atom pattern with their corresponding distances
					for atom in processed_cloud_02:
						atom_pattern_N+= atom_reference[atom[0].name + atom[0].resname] + ':'
						atom_pattern_dist += str("{0:.3f}".format(distance(res[i].H.pos, atom[0].pos))) + '|' + \
							str("{0:.3f}".format(distance(res[i-1].N.pos, atom[0].pos))) + '|' + \
							str("{0:.3f}".format(distance(res[i-2].O.pos, atom[0].pos))) + '|' + \
							str("{0:.3f}".format(distance(res[ i ].O.pos, atom[0].pos))) + '|' + \
							str("{0:.3f}".format(distance(res[i+1].N.pos, atom[0].pos))) + '|'

			# Write the output to the open file
			open_N_files[f].write(atom_pattern_N + '|' + str(ts.frame) + '|' + \
				str("{0:.3f}".format(distance(res[i-1].N.pos, res[i].H.pos))) + '|' + \
				str("{0:.3f}".format(distance(res[i-2].O.pos, res[i].H.pos))) + '|' + \
				str("{0:.3f}".format(distance(res[i-1].O.pos, res[i].H.pos))) + '|' + \
				str("{0:.3f}".format(distance(res[ i ].O.pos, res[i].H.pos))) + '|' + \
				str("{0:.3f}".format(distance(res[i+1].N.pos, res[i].H.pos))) + '|' + \
				atom_pattern_dist + '\n')

			# Intentionally offset the oxygen output	
			o = i - 1
			if (o + 2) < len(res)-1:
				######################################################
				##  Process the O dome								##
				######################################################

				# Select all atoms within 3.9 A of target
				cloud = u.select_atoms("around 3.9 atom SYSTEM " + str(res[o].id) + " O")
				processed_cloud = []
				processed_cloud_02 = []
				atom_pattern_O = ''
				atom_pattern_dist = ''

				# First round of processing
				for atom in cloud:
					ref_name = atom.name + atom.resname

					# Check to ensure backbone atoms within (+/-) 2 residues are not included
					if atom.name in bb_atoms and abs(res[o].id - atom.resid) < 3:
						continue
					else:
						if ref_name in atom_reference:
							# Eliminate all Carbon atoms
							if "C" in atom_reference[ref_name]:
								continue
							else:
								x = res[o].N.pos - res[o].O.pos
								v = atom.pos - res[o].O.pos

								# Eliminate any atoms that have vector angles less than 90.0
								if abs(vangle(x, v)) < 90.0:
									continue
								else:
									if ref_name in atom_reference:
										processed_cloud.append([atom, distance(atom.pos, res[o].O.pos)])

				# Second round of processing
				if processed_cloud is None:
					atom_pattern_O = "Z:"
				else:
					processed_cloud_02 = find_closest_atom(processed_cloud, res[o].O)
			
					if processed_cloud_02 == None or len(processed_cloud_02) < 1:
						atom_pattern_O = "Z:"
					else:
						# Rank atoms based on the atom type preferences noted above
						processed_cloud_02.sort(key = lambda atom: atom_ranking[atom_reference[atom[0].name + atom[0].resname]])
						
						# Generate the atom pattern with their corresponding distances
						for atom in processed_cloud_02:
							atom_pattern_O += atom_reference[atom[0].name + atom[0].resname] + ":"
							atom_pattern_dist += str("{0:.3f}".format(distance(res[o].O.pos, atom[0].pos))) + '|' + \
								str("{0:.3f}".format(distance(res[o].N.pos, atom[0].pos))) + '|' + \
								str("{0:.3f}".format(distance(res[o-1].O.pos, atom[0].pos))) + '|' + \
								str("{0:.3f}".format(distance(res[o+1].O.pos, atom[0].pos))) + '|' + \
								str("{0:.3f}".format(distance(res[o+2].N.pos, atom[0].pos))) + '|'
				
				# d3 becomes d1 in the next residues
				d3 = str("{0:.3f}".format(distance(res[o+1].O.pos, res[o].O.pos)))
				next_res_dist.append(d3) 

				# Write the output to the open file
				open_O_files[f].write(atom_pattern_O + '|' + str(ts.frame) + '|' + \
					str("{0:.3f}".format(distance(res[ o ].N.pos, res[o].O.pos))) + '|' + \
					str("{0:.3f}".format(distance(res[o-1].O.pos, res[o].O.pos))) + '|' + \
					str("{0:.3f}".format(distance(res[o+1].N.pos, res[o].O.pos))) + '|' + \
					d3 + '|' + str("{0:.3f}".format(distance(res[o+2].N.pos, res[o].O.pos))) + '|' + \
					atom_pattern_dist + '\n')
		elif i > 2 and res[i].name not in "PRO":
			######################################################
			##  Process the N dome								##
			######################################################
			
			# Select all atoms within 5.0 A of target
			cloud = u.select_atoms("around 5.0 atom SYSTEM " + str(res[i].id) + " H")
			processed_cloud = []
			processed_cloud_02 = []
			atom_pattern_N = ''
			atom_pattern_dist = ''

			# First round of processing
			for atom in cloud:
				ref_name = atom.name + atom.resname
				
				# Check to ensure backbone atoms within (+/-) 2 residues are not included
				if atom.name in bb_atoms and abs(res[i].id - atom.resid) < 3:
					continue
				else:
					if ref_name in atom_reference:
						x = res[i].N.pos - res[i].H.pos
						v = atom.pos - res[i].H.pos
						
						# Eliminate any atoms that have vector angles less than 90.0
						if abs(vangle(x, v)) < 90.0:
							continue
						else:
							######################################################
							##  Compute distance between atoms based on types	##
							##  Aromatic atoms: within 5.0 A 					##
							##  Carbon atoms: 	within 4.0 A 					##
							##  Other atoms: 	wtihin 2.5 A 					##
							######################################################
							if atom_reference[ref_name] in "R":
								processed_cloud.append([atom, distance(atom.pos, res[i].H.pos)])
							elif atom_reference[ref_name] in "C":
								C_dist = distance(atom.pos, res[i].H.pos)
								if C_dist <= 4.0:
									processed_cloud.append([atom, C_dist])
							else:
								dist = distance(atom.pos, res[i].H.pos)
								if dist <= 2.5:
									processed_cloud.append([atom, dist])

			# Second round of processing
			if processed_cloud is None:
				atom_pattern_N = "Z:"
			else:
				processed_cloud_02 = find_closest_atom(processed_cloud, res[i].H)
			
				if processed_cloud_02 == None:
					atom_pattern_N = "Z:"
				else:
					# Rank atoms based on the atom type preferences noted above
					processed_cloud_02.sort(key = lambda atom: atom_ranking[atom_reference[atom[0].name + atom[0].resname]])
					
					# Generate the atom pattern with their corresponding distances
					for atom in processed_cloud_02:
						atom_pattern_N+= atom_reference[atom[0].name + atom[0].resname] + ':'
						atom_pattern_dist += str("{0:.3f}".format(distance(res[i].H.pos, atom[0].pos))) + '|' + \
							str("{0:.3f}".format(distance(res[i-1].N.pos, atom[0].pos))) + '|' + \
							str("{0:.3f}".format(distance(res[i ].CA.pos, atom[0].pos))) + '|' + \
							str("{0:.3f}".format(distance(res[ i ].O.pos, atom[0].pos))) + '|' + \
							str("{0:.3f}".format(distance(res[i+1].N.pos, atom[0].pos))) + '|'

			# Write the output to the open file
			open_N_files[f].write(atom_pattern_N + '|' + str(ts.frame) + '|' + \
				str("{0:.3f}".format(distance(res[i-1].N.pos, res[i].H.pos))) + '|' + \
				str("{0:.3f}".format(distance(res[i-2].O.pos, res[i].H.pos))) + '|' + \
				str("{0:.3f}".format(distance(res[i-1].O.pos, res[i].H.pos))) + '|' + \
				str("{0:.3f}".format(distance(res[ i ].O.pos, res[i].H.pos))) + '|' + 
				str("{0:.3f}".format(distance(res[i+1].N.pos, res[i].H.pos))) + '|' + \
				atom_pattern_dist + '\n')

			# Intentionally offset the oxygen output
			o = i - 1
			if o+2 < len(res)-1:
				######################################################
				##  Process the O dome								##
				######################################################

				# Select all atoms within 3.9 A of target
				cloud = u.select_atoms("around 3.9 atom SYSTEM " + str(res[o].id) + " O")
				processed_cloud = []
				processed_cloud_02 = []
				atom_pattern_O = ''
				atom_pattern_dist = ''
			
				# First round of processing
				for atom in cloud:
					ref_name = atom.name + atom.resname

					# Check to ensure backbone atoms within (+/-) 2 residues are not included
					if atom.name in bb_atoms and abs(res[o].id - atom.resid) < 3:
						continue
					else:
						if ref_name in atom_reference:
							# Eliminate all Carbon atoms
							if "C" in atom_reference[ref_name]:
								continue
							else:
								x = res[o].N.pos - res[o].O.pos
								v = atom.pos - res[o].O.pos

								# Eliminate any atoms that have vector angles less than 90.0
								if abs(vangle(x, v)) < 90.0:
									continue
								else:
									if ref_name in atom_reference:
										processed_cloud.append([atom, distance(atom.pos, res[o].O.pos)])
				
				# Second round of processing
				if processed_cloud is None:
					atom_pattern_O = "Z:"
				else:
					processed_cloud_02 = find_closest_atom(processed_cloud, res[o].O)
			
					if processed_cloud_02 == None or len(processed_cloud_02) < 1:
						atom_pattern_O = "Z:"
					else:
						# Rank atoms based on the atom type preferences noted above
						processed_cloud_02.sort(key = lambda atom: atom_ranking[atom_reference[atom[0].name + atom[0].resname]])
						
						# Generate the atom pattern with their corresponding distances
						for atom in processed_cloud_02:
							atom_pattern_O += atom_reference[atom[0].name + atom[0].resname] + ":"
							atom_pattern_dist += str("{0:.3f}".format(distance(res[o].O.pos, atom[0].pos))) + '|' + \
								str("{0:.3f}".format(distance(res[o].N.pos, atom[0].pos))) + '|' + \
								str("{0:.3f}".format(distance(res[o-1].O.pos, atom[0].pos))) + '|' + \
								str("{0:.3f}".format(distance(res[o+1].O.pos, atom[0].pos))) + '|' + \
								str("{0:.3f}".format(distance(res[o+2].N.pos, atom[0].pos))) + '|'

				# d3 becomes d1 in the next residue
				d3 = str("{0:.3f}".format(distance(res[o+1].O.pos, res[o].O.pos)))

				# Write the output to the open file
				open_O_files[f].write(atom_pattern_O + '|' + str(ts.frame) + '|' + \
					str("{0:.3f}".format(distance(res[ o ].N.pos, res[o].O.pos))) + '|' + next_res_dist[0] + '|' + \
					str("{0:.3f}".format(distance(res[o+1].N.pos, res[o].O.pos))) + '|' + d3 + '|' + \
					str("{0:.3f}".format(distance(res[o+2].N.pos, res[o].O.pos))) + '|' + \
					atom_pattern_dist + '\n')
				next_res_dist[0] = d3
		else:
			continue

print "Program time: " + str("{0:.3f}".format(time.time() - start_time)) + " seconds."
