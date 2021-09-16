#!/usr/bin/env python3

'''Run shiftx on a trajectory/topology and output in the same .shift
   format that we use.'''
   
import numpy as np
import glob, re 
import argparse
import sys, os, gzip, tempfile, shutil, subprocess
from collections import defaultdict
import MDAnalysis

parser = argparse.ArgumentParser()
parser.add_argument('--topo',help="topology file of simulation (e.g. prmtop)")
parser.add_argument('--traj',help="trajectory file of simulation containing explicit wrapped waters (e.g. dcd)")
parser.add_argument('--shiftx',help='location of shiftx distribution',default='$HOME/build/shiftx2-v109A-linux-20150821/')
parser.add_argument('--out',help='output file (default stdout)',type=argparse.FileType('w'),default=sys.stdout)

args = parser.parse_args()
    
# create tmp directory
tmpdir = tempfile.mkdtemp()

# dump trajectory to pdb files, cpptraj is fastest at this
p = subprocess.Popen('cpptraj > /dev/null',stdin=subprocess.PIPE,shell=True)
cmds = '''parm {0}
trajin {1}
strip :WAT,Cl-,Na+
trajout {2}/out pdb multi nobox
run
quit'''.format(args.topo,args.traj,tmpdir)
p.communicate(cmds)

# need to rename some residues
for fname in glob.glob('%s/out*'%tmpdir):
    subprocess.call("sed 's/HIE/HIS/g' -i %s"%fname,shell=True)
    subprocess.call("sed 's/CYX/CYS/g' -i %s"%fname,shell=True)


# run shiftx in batch mode
cmd = "java  -Xmx1900m  -cp {0}/bin:{0}/lib/weka.jar ShiftXp -b '{1}/out*' -atoms BACKBONE -ph 5.0 -temp 298.0 -dir {0} > /dev/null".format(args.shiftx,tmpdir)
subprocess.call(cmd,shell=True)

# read and collate shiftx output

data = {'H': defaultdict(list), 'N': defaultdict(list), 'CA': defaultdict(list)}
resmap = {}
for fname in glob.glob('%s/*.sxp'%tmpdir):
    f = open(fname)
    f.readline() #header
    for line in f:
        (resid,resname,rname,atom,shift) = line.split(',')
        if atom in data:
            data[atom][int(resid)].append(float(shift))
            resmap[int(resid)] = resname
            
# output mean and stddev
o = args.out
o.write('#ID\tName\tN\tH\tC\tNstd\tHstd\tCstd\tCoverage\n')
for resi in sorted(resmap.keys()):
    if data['N'][resi]: #e.g. PRO doesn't have values
        o.write('%d\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t0.0\n'%(resi,resmap[resi],
            np.mean(data['N'][resi]),np.mean(data['H'][resi]),np.mean(data['CA'][resi]),
            np.std(data['N'][resi]),np.std(data['H'][resi]),np.std(data['CA'][resi])))
        
shutil.rmtree(tmpdir)
