#!/usr/bin/env python3

# Script that analyzes a molecular dynamics simulation for experimental validation.

import cProfile
import pstats, os
import io
import sys, MDAnalysis
import math
import numpy as np
import time
import argparse
from os.path import splitext
from MDAnalysis.analysis.rms import *
from MDAnalysis.analysis.distances import *
import gzip
from scipy import spatial


def addindex(s, r, aname):
	if aname in r.atoms.names:
		i = np.argmax(r.atoms.names == aname)
		a = r.atoms[i]
		s.add(a.index)

			
# Dictionary that maps atom names to groups
atom_reference = {
"NALA": "N", "NARG": "N", "NASN": "N", "NASP": "N", "NCYS": "N","NCYX": "N", "NGLN": "N", "NGLY": "N", "NGLU": "N", "NHIS": "N","NHID": "N", "NHIE": "N", "NHIP": "N",
"NHSD": "N", "NHSE": "N", "NHSP": "N", "NILE": "N", "NLEU": "N", "NLYS": "N", "NMET": "N", "NPHE": "N",
"NPRO": "N", "NSER": "N", "NTHR": "N", "NTRP": "N", "NTYR": "N", "NVAL": "N",

"OALA": "O", "OARG": "O", "OASN": "O", "OASP": "O", "OCYS": "O", "OCYX": "O", "OGLN": "O", "OGLY":   "O", "OGLU":   "O", "OHIS": "O","OHID": "O", "OHIE": "O", "OHIP": "O",
"OHSD": "O", "OHSE": "O", "OHSP": "O", "OILE": "O", "OLEU": "O", "OLYS": "O", "OMET":   "O", "OPHE":   "O",
"OPRO": "O", "OSER": "O", "OTHR": "O", "OTRP": "O", "OTYR": "O", "OVAL": "O", "OD1ASN": "O", "OE1GLN": "O",

"OH2TIP3": "W", "OWAT": "W",

"CD1PHE": "R", "CD1TRP": "R", "CD1TYR": "R", "CD2HSD": "R", "CD2HSE": "R", "CD2HSP": "R",  "CD2HIS": "R","CD2HID": "R", "CD2HIE": "R", "CD2HIP": "R",
"CD2PHE": "R", "CD2TRP": "R", "CD2TYR": "R", "CE1HSD": "R", "CE1HSE": "R", "CE1HSP": "R", "CE1HIS": "R","CE1HID": "R", "CE1HIE": "R", "CE1HIP": "R",
"CE1PHE": "R", "CE1TYR": "R", "CE2PHE": "R", "CE2TRP": "R", "CE2TYR": "R", "CE3TRP": "R",
"CGHSD":  "R", "CGHSE":  "R", "CGHSP":  "R", "CGHIS":  "R","CGHID":  "R", "CGHIE":  "R", "CGHIP":  "R", "CGPHE":  "R", "CGTRP":  "R", "CGTYR":  "R",
"CH2TRP": "R", "CZ2TRP": "R", "CZ3TRP": "R", "CZPHE":  "R", "CZTYR":  "R", "ND1HSD": "R",
"ND1HSE": "R", "ND1HSP": "R", "ND1HIS": "R","ND1HID": "R", "ND1HIE": "R", "ND1HIP": "R",
"NE1TRP": "R", "NE2HSD": "R", "NE2HSE": "R", "NE2HSP": "R", "NE2HIS": "R","NE2HID": "R", "NE2HIE": "R", "NE2HIP": "R",

"OG1THR": "L", "OGSER":  "L", "OHTYR":  "L", "OD2ASPH": "L", "OE2GLUH": "L", 

"ND2ASN": "D", "NE2GLN": "D",

"SDMET":  "S", "SGCYS":  "S", "SGCYX": "S",

"OD1ASP": "A", "OD2ASP": "A", "OE1GLU": "A", "OE2GLU": "A", "OT1ALA": "A", "OT1ARG": "A", "OT1ASN": "A", 
"OT1ASP": "A", "OT1CYS": "A","OT1CYX": "A", "OT1GLN": "A", "OT1GLY": "A", "OT1GLU": "A", "OT1HSD": "A", "OT1HSE": "A", 
"OT1HSP": "A", "OT1HIS": "A", "OT1HID": "A", "OT1HIE": "A", "OT1HIP": "A",
 "OT1ILE": "A", "OT1LEU": "A", "OT1LYS": "A", "OT1MET": "A", "OT1PHE": "A", "OT1PRO": "A", 
"OT1SER": "A", "OT1THR": "A", "OT1TRP": "A", "OT1TYR": "A", "OT1VAL": "A",

"NZLYS": "B", "NTALA": "B", 

"NEARG": "G", "NH1ARG": "G", "NH2ARG": "G", 

#we intentionally ignore ions!
#"SODSOD": "P", 
#"Na+Na+": "P", 

#"CLACLA": "M",
#"Cl-Cl-": "M",

"Z": "Z", 

"OD1ASPH": "U", "OE1GLUH": "U", 
"OD1ASH": "U", "OE1GLH": "U", 

"CALA": "C", "CARG": "C", "CASN": "C", "CASP": "C", "CCYS": "C","CCYX": "C", "CGLN": "C", "CGLY": "C", "CGLU": "C", 
"CHSD": "C", "CHSE": "C", "CHSP": "C","CHIS": "C", "CHID": "C", "CHIE": "C", "CHIP": "C",
"CILE": "C", "CLEU": "C", "CLYS": "C", "CMET": "C", "CPHE": "C", 
"CPRO": "C", "CSER": "C", "CTHR": "C", "CTRP": "C", "CTYR": "C", "CVAL": "C"
}

# Dictionary that maps atom type to ranking preference
atom_ranking = {
	"N": 0, "O": 1, "W": 2, "R":  3, "L":  4, "D":  5, "S":  6, "A":  7,
	"B": 8, "G": 9, "P":10, "M": 11, "Z": 12, "U": 13, "C": 14
}

#the following are global sets
isNT = None 
isCT = None
def atom_to_pattern(atom):
	#treat N-term and C-term specially
	if atom.index in isNT:
		return 'B'
	elif atom.index in isCT:
		return 'A'
	name = atom.name + atom.resname
	if name in atom_reference:
		return atom_reference[name]
	else:
		return None
	

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

	if udxvdx > 1.0:
		udxvdx = 1.0
	if udxvdx < -1.0:
		udxvdx = -1.0
	
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
		ref_name = atom_to_pattern(atom[0])
		if ref_name == 'C':
			carbons.append(atom)
		elif atom[0].name == 'O' and ref_name == 'O':
			oxygens.append(atom)
		else:
			other_atoms.append(atom)

	# Finds closest oxygen to current residue target atom
	for atom in oxygens:
		dist = distance(atom[0].position, target_atom.position)
		if dist < minO:
			minO = dist
			nearestO = atom
	
	# If there are no O atoms or any other atoms, return no atoms (None)
	if nearestO == None and len(other_atoms) == 0:
		return None
	# If there are no O atoms and there are other atoms, return the other atoms
	elif nearestO == None and len(other_atoms) > 0:
		return other_atoms
	
	# Add all remaining oxygens
	nearest_atoms += oxygens
	
	# Finds closest C to previously found nearestO
	if len(carbons) > 0:
		for atom in carbons:
			dist = distance(atom[0].position, nearestO[0].position)
			if dist < minC:
				minC = dist
				nearestC = atom
	
	if nearestC:
		nearest_atoms.append(nearestC)

	if len(other_atoms) > 0:
		for atom in other_atoms:
			nearest_atoms.append(atom)
	
	return nearest_atoms

def H(r):
	i = np.argmax(r.atoms.names == 'H')
	return r.atoms[i]
	
def O(r):
	i = np.argmax(r.atoms.names == 'O')
	return r.atoms[i]
		
def N(r):
	i = np.argmax(r.atoms.names == 'N')
	return r.atoms[i]
	
def process_frame(u, residues, prot, notH):
	'''Calculate distance descriptors for all requested residues of the current
	frame of a molecular dynamics simulation. Provide selections of the protein atoms and all notH atoms'''
	
	bb_atoms = ["N", "O"]
	res = prot.atoms.residues
	# Load this frame into a spatial cKDTree for fast nearest-neighbor search
	tree = spatial.cKDTree(np.array([x.position for x in notH.atoms]))
	
	ret = dict() #indexed by resid
	# Simultaneously iterate through each residue and corresponding open file
	for r in residues:
		i = r.resid-1;
		assert res[i] == r
		######################################################
		##  Process the N dome								
		######################################################
		
		# Original atom selection
		# orig_cloud = u.select_atoms("around 5.0 atom SYSTEM " + str(res[i].id) + " H")

		cloud = []
		processed_cloud = []
		atom_pattern_N = ''
		atom_pattern_dist = []
		current_res_H_pos = H(r).position
		x = N(r).position - current_res_H_pos
		res_i_O_pos = O(res[ i ]).position
		res_i_m1_N_pos = N(res[i-1]).position
		res_i_p1_N_pos = N(res[i+1]).position
		res_i_m2_O_pos = O(res[i-2]).position
		
		# First round of processing from enhanced atom selection
		for j in tree.query_ball_point(current_res_H_pos, 5.0):
			atom = notH.atoms[j]
			atom_pos = atom.position
			ref_name = atom_to_pattern(atom) 
			# Check to ensure backbone atoms within (+/-) 2 residues are not included
			if atom.name in bb_atoms and abs(r.resid - atom.resid) < 3:
				continue 
			else:
				if ref_name:
					# Variables for calculating the vector angles of each atom
					v = atom_pos - current_res_H_pos

					# Eliminate any atoms that have vector angles less than 90.0
					if abs(vangle(x, v)) < 90.0:
						continue
					else:
						######################################################
						##  Compute distance between atoms based on types	
						##  Aromatic atoms: within 5.0 A 					
						##  Carbon atoms: 	within 4.0 A 					
						##  Other atoms: 	wtihin 2.5 A 					
						######################################################
						dist = distance(atom_pos, current_res_H_pos)
						
						if ref_name == "R":
							cloud.append([atom, dist])
						elif ref_name == "C":
							if dist <= 4.0:
								cloud.append([atom, dist])
						elif dist <= 2.5:
							cloud.append([atom, dist])
		
		# Second round of processing					
		if cloud is None:
			atom_pattern_N = "Z:"
			atom_pattern_dist += [0.0,0.0,0.0,0.0,0.0]
		else:
			processed_cloud = find_closest_atom(cloud, H(res[i]))
	
			if processed_cloud == None:
				atom_pattern_N = "Z:"
				atom_pattern_dist += [0.0,0.0,0.0,0.0,0.0]
			else:
				# Rank atoms based on the atom type preferences noted above
				processed_cloud.sort(key = lambda atom: (atom_ranking[atom_to_pattern(atom[0])],atom[1]))				
				
				# Generate the atom pattern with their corresponding distances
				for atom in processed_cloud:
					atom_pattern_N += atom_to_pattern(atom[0]) + ':'
					apos = atom[0].position
					atom_pattern_dist += [
						distance(current_res_H_pos, apos),distance(res_i_m1_N_pos, apos),
						distance(res_i_m2_O_pos, apos),distance(res_i_O_pos, apos),
						distance(res_i_p1_N_pos, apos)]

		# Write the output to the open file
		nvals = [atom_pattern_N,u.trajectory.frame,
				distance(res_i_m1_N_pos, current_res_H_pos),distance(res_i_m2_O_pos, current_res_H_pos),
				distance(O(res[i-1]).position, current_res_H_pos),distance(res_i_O_pos, current_res_H_pos),
				distance(N(res[i+1]).position, current_res_H_pos)] + atom_pattern_dist
		
		# Intentionally offset the oxygen output	
		o = i - 1
		######################################################
		##  Process the O dome								
		######################################################

		# Original atom selection
		# orig_cloud = u.select_atoms("around 3.9 atom SYSTEM " + str(res[o].id) + " O")
		
		cloud = []
		processed_cloud = []
		atom_pattern_O = ''
		atom_pattern_dist = []
		current_res_O_pos = O(res[o]).position
		x = N(res[o+1]).position - current_res_O_pos
		res_o_N_pos = N(res[ o ]).position
		res_o_m1_O_pos = O(res[o-1]).position
		res_o_p1_O_pos = O(res[o+1]).position
		res_o_p2_N_pos = N(res[o+2]).position
		
		# First round of processing from enhanced atom selection
		for j in tree.query_ball_point(current_res_O_pos, 3.9):
			atom = notH.atoms[j]
			atom_pos = atom.position
			ref_name = atom_to_pattern(atom)

			# Check to ensure backbone atoms within (+/-) 2 residues are not included
			if atom.name in bb_atoms and abs(res[o].resid - atom.resid) < 3:
				continue
			else:
				if ref_name:
					# Eliminate all Carbon atoms
					if "C" == ref_name:
						continue
					else:
						# Variables for calculating the vector angles of each atom
						v = atom_pos - current_res_O_pos

						# Eliminate any atoms that have vector angles less than 90.0
						if abs(vangle(x, v)) < 90.0:
							continue
						elif ref_name:
								cloud.append([atom, distance(atom_pos, current_res_O_pos)])

		# Second round of processing
		if cloud is None:
			atom_pattern_O = "Z:"
			atom_pattern_dist += [0.0,0.0,0.0,0.0,0.0]
		else:
			processed_cloud = find_closest_atom(cloud, O(res[o]))
	
			if processed_cloud == None or len(processed_cloud) < 1:
				atom_pattern_O = "Z:"
				atom_pattern_dist += [0.0,0.0,0.0,0.0,0.0]
			else:
				# Rank atoms based on the atom type preferences noted above
				processed_cloud.sort(key = lambda atom: (atom_ranking[atom_to_pattern(atom[0])], atom[1]))
				
				# Generate the atom pattern with their corresponding distances
				for atom in processed_cloud:
					atom_pattern_O += atom_to_pattern(atom[0]) + ":"
					apos = atom[0].position
					atom_pattern_dist += [distance(current_res_O_pos,   apos),distance(res_o_N_pos,   apos),
						distance(res_o_m1_O_pos, apos),distance(res_o_p1_O_pos, apos),
						distance(res_o_p2_N_pos, apos)]
		
		# Write the output to the open file
		ovals = [atom_pattern_O,u.trajectory.frame,
			distance(res_o_N_pos, current_res_O_pos),distance(res_o_m1_O_pos, current_res_O_pos),
			distance(N(res[o+1]).position, current_res_O_pos),distance(res_o_p1_O_pos, current_res_O_pos),
			distance(res_o_p2_N_pos, current_res_O_pos)] + atom_pattern_dist
		
		ret[i] = (nvals,ovals)
		
	return ret	


def process_residues(prot):
	'''process residues of an mdanalysis protein selection.
	This must be called to initialize stage before processing frame.
	TODO: refactor into class object to eliminate global(!) state'''
	global isNT
	global isCT
	isNT = set()
	isCT = set()
	startch = 0
	res = prot.atoms.residues
	residues = []
	for (i,r) in enumerate(res):
		if i == startch:
			addindex(isNT, r, 'N')
		if i >= startch+2:
			if ('OXT' in r.atoms.names or 'OT2' in r.atoms.names) and len(residues) > 0: #end of chain
				addindex(isCT, r, 'OXT')
				addindex(isCT, r, 'OT1')
				addindex(isCT, r, 'OT2')
				addindex(isCT, r, 'O')
				
				startch = i+1 #update start of chain			
			elif r.resname != 'PRO':
				residues.append(r)
	return residues
		
if __name__ == '__main__':
	start_time = time.time()
	
	# User-friendly argument handling
	parser = argparse.ArgumentParser()
	parser.add_argument("--topo", required=True,  help="input file with (.prmtop) extension")
	parser.add_argument("--traj", required=True, help="input file with (.dcd) extension")
	parser.add_argument("-r","--residues", help="space deliminated residues to compute (all if unspecified)",nargs='*')
	parser.add_argument("--start",help="frame to start on",type=int,default=0)
	parser.add_argument("--end",help="frame to end at",type=int,default=-1)
	parser.add_argument("--flush",help="flush I/O as it is written",action='store_true')
	parser.add_argument("--outputdir",help="directory in which to put DUMP files",default="output")
	args = parser.parse_args()
	
	# Set up universe, selections and data structures
	u = MDAnalysis.Universe(args.topo, args.traj)
	prot = u.select_atoms("protein")
	notH = u.select_atoms("not (name H*)")
	bb_atoms = ["N", "O"]
	res = prot.atoms.residues
	
	#identify residues of interest
	startch = 0
	residues = []
	
	residues = process_residues(prot)	
	if args.residues:
		residues = [] #override with user specified
		for r in args.residues:
			residues.append(res[int(r)-1])
	######################################################################################################################
	## Main code - opens all files, iterates through each timestep once while iterating through residues.               ##
	######################################################################################################################
	
	# Holds all open files
	open_O_files = []
	open_N_files = []
	
	if not os.path.exists(args.outputdir):
		os.makedirs(args.outputdir)
	
	# Open all files at the beginning
	for r in residues:
		i = r.id
		file_name_O = '%s/%s.O.%d.gz' % (args.outputdir,r.name,i)
		file_name_N = '%s/%s.N.%d.gz' % (args.outputdir,r.name,i)
		
		f_O = gzip.open(file_name_O,'w')
		f_N = gzip.open(file_name_N,'w')
		
		open_O_files.append(f_O)
		open_N_files.append(f_N)
	
	end = args.end
	if end < 0: end = u.trajectory.n_frames
	
	# Iterate through each desired frame
	prot = u.select_atoms("protein")
	notH = u.select_atoms("not (name H*)")
	
	for ts in u.trajectory[args.start:end]:
		
		print("Processing frame #: " + str(ts.frame))
		resdata = process_frame(u, residues, prot, notH)

		for r, Nfile, Ofile in zip(residues, open_N_files, open_O_files):
			i = r.id-1;
			assert res[i] == r
			assert i in resdata
			(ndata,odata) = resdata[i]
			# Write the output to the open file
			Nfile.write("{0}|{1}|{2:.5f}|{3:.5f}|{4:.5f}|{5:.5f}|{6:.5f}".format(*ndata[0:7]))
			for x in ndata[7:]:
				Nfile.write('|{0:.5f}'.format(x))
			Nfile.write('\n')
			if args.flush: Nfile.flush()
				
			# Write the output to the open file
			Ofile.write("{0}|{1}|{2:.5f}|{3:.5f}|{4:.5f}|{5:.5f}|{6:.5f}".format(*odata[0:7]))
			for x in odata[7:]:
				Ofile.write('|{0:.5f}'.format(x))
			Ofile.write('\n')
			if args.flush: Ofile.flush()

	
	print("Program time: {0:.3f} seconds.".format(time.time() - start_time))
