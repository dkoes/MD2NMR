#!/usr/bin/env python

'''Calculate chemical shifts of a directory of residue traces'''

import numpy as np
import glob, re 
import argparse
import sys, os, gzip
import shiftres, dump
import multiprocessing
import MDAnalysis
from multiprocessing import Pool


def shifts_from_file(db,fname,resid,resname,args):
    '''return summary statistics for a file'''
    ndata = shiftres.read_resdump(fname)
    oname = fname.replace('.N.','.O.')
    odata = shiftres.read_resdump(oname)
    shifts = np.array(shiftres.compute_shifts(db[resname],ndata,odata,args,args.verbose))
    means = shifts[:,1:].mean(axis=0)
    stds = shifts[:,1:].std(axis=0)
    return (resid,resname,means,stds,shifts)

def thread_call(a):
    return shifts_from_file(db,a[0],a[1],a[2],args)

def doframes(framerange):
    (start,end) = framerange
    u = MDAnalysis.Universe(args.topo, args.traj)
    prot = u.select_atoms("protein")
    notH = u.select_atoms("not (name H*)")
    residues = dump.process_residues(prot)    

    ret = []
    for _ in u.trajectory[start:end]:
        ret.append(dump.process_frame(u, residues, prot, notH))
        
    return ret
    
def doshifts(rinfo):
    (r,rname) = rinfo
    return shifts_from_dicts(db, resdata, r, rname, args)
            
def shifts_from_dicts(db,resd,r, resname,args):
    '''return summary statistics for residue r from a list of dump dicts'''
    ndata = dict()
    odata = dict()
    for d in resd:
        #convert format
        nvals = d[r][0]
        pattern = nvals[0]
        frame = nvals[1]
        values = np.array(nvals[2:])
        assert frame not in ndata
        ndata[frame] = (pattern,values)
        
        ovals = d[r][1]
        pattern = ovals[0]
        frame = ovals[1]
        values = np.array(ovals[2:])
        assert frame not in odata
        odata[frame] = (pattern,values)
        
    shiftinfo = shiftres.compute_shifts(db[resname],ndata,odata,args,args.verbose)
    shifts = np.array(shiftinfo)
    means = shifts[:,1:].mean(axis=0)
    stds = shifts[:,1:].std(axis=0)
    return (r+1,resname,means,stds,shifts)        

        

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d',"--database",help="database file or directory containing .tri files",required=True)
    parser.add_argument('-i','--inputdir',help="directory containing dump files")
    parser.add_argument('--resdir',help='directory to output per-residue, per-frame shift calculations')
    parser.add_argument('--topo',help="topology file of simulation (e.g. prmtop)")
    parser.add_argument('--traj',help="trajectory file of simulation containing explicit wrapped waters (e.g. dcd)")
    parser.add_argument('--out',help='output file (default stdout)',type=argparse.FileType('w'),default=sys.stdout)
    parser.add_argument('--resmap',help='file containing mapping from resid to resname',default='')
    parser.add_argument('-v','--verbose',help='output informative messages',action='store_true')
    parser.add_argument('--cpus',help='number of cores to use',type=int,default=multiprocessing.cpu_count())
    parser.add_argument("-r","--residues", help="space deliminated residues to compute (all if unspecified)",nargs='*')    
    shiftres.add_shift_args(parser)
    args = parser.parse_args()
    
    db = shiftres.read_db(args.database)
    
    #read in residue mapping if present
    resmap = dict()
    if args.resmap:
        for line in open(args.resmap):
            (i,resn) = line.split()
            resmap[int(i)] = resn
        
    if args.inputdir:
        files = glob.glob('%s/*.N.*[0-9]'%args.inputdir)
        if not files:
            files = glob.glob('%s/*.N.*[0-9].gz'%args.inputdir)
        if not files:
            sys.stderr.write("No compatible files found in %s\n"%args.inputdir)
            sys.exit(-1)
        #setup inputs for parallel processing
        inputs = []
        for fname in files:
            #are files names with resname?
            m = re.search(r'([A-Za-z]{3})\.N\.(\d+)',fname)
            if not m:
                continue
            resname = None
            i = int(m.group(2))
            if m.group(1) == 'UMP': # a DUMP file
                if i in resmap:
                    resname = resmap[i]
            else:
                resname = m.group(1)
    
            if not resname:
                sys.stderr.write("Could not identify residue name for %s, need to use --resmap?\n" % fname)
                sys.exit(1)
                
            inputs.append((fname,i,resname))
            
        #count lines
        if fname.endswith('.gz'):
            for n, l in enumerate(gzip.open(fname),1): pass
        else:
            for n, l in enumerate(open(fname),1): pass
        
        #sort by resid
        inputs.sort(key=lambda (fname,i,resname):i)    
        
        #setup pool
        pool = Pool(args.cpus)
        
        results = pool.map(thread_call,inputs,1)    
        
    elif args.topo and args.traj:
        #read entirety of distance descriptors into memory
        u = MDAnalysis.Universe(args.topo, args.traj)
        prot = u.select_atoms("protein")
        notH = u.select_atoms("not (name H*)")
        res = prot.atoms.residues
        residues = dump.process_residues(prot)
        if args.residues:
            residues = [] #override with user specified
            for r in args.residues:
                residues.append(res[int(r)-1])
        pool = Pool(args.cpus)

        #parallelize over frames
        n = len(u.trajectory)
        step = int(n / args.cpus)+1
        franges = [[i,i+step] for i in xrange(0,n,step)]
        franges[-1][1]=n

        resdata = pool.map(doframes, franges)

        #flatten partitions
        resdata = [item for sublist in resdata for item in sublist]
        pool.close()

        #create new pool that knows about resdata    
        pool = Pool(args.cpus)

        #compute shifts            
        results = pool.map(doshifts, [(r.id-1,r.name) for r in residues])
        
    else:
        sys.stderr("Require --inputdir of dump files or topology/trajectory information\n")
        sys.exit(-1)

    if args.resdir and not os.path.exists(args.resdir):
        os.makedirs(args.resdir)
        
    args.out.write('#ID\tName\tN\tH\tC\tNstd\tHstd\tCstd\tCoverage\n')    
    for (resid,resname,means,stds,shifts) in results:
        args.out.write('%d\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n' % (resid,resname,means[0],means[1],means[2],stds[0],stds[1],stds[2],float(len(shifts))/n))
        if args.resdir: #output individual residue data
            out = gzip.open('%s/%s.%d.shifts.gz' % (args.resdir,resname,resid),'w')
            out.write('#Frame\tN\tH\tC\tNdist\tOdist\n')
            for sh in shifts:
                out.write('%d\t'%sh[0])
                out.write('\t'.join(map(lambda x: '%.5f'%x,sh[1:]))+'\n')
            
        
