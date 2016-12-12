#!/usr/bin/env python

'''Calculate chemical shifts of a directory of residue traces'''

import numpy as np
import glob, re 
import argparse
import sys
import shiftres
import multiprocessing
from multiprocessing import Pool


def shifts_from_file(db,fname,resname,thresholds,verbose):
    '''return summary statistics for a file'''
    ndata = shiftres.read_resdump(fname)
    oname = fname.replace('.N.','.O.')
    odata = shiftres.read_resdump(oname)
    shifts = np.array(shiftres.compute_shifts(db[resname],ndata,odata,thresholds,verbose))[:,1:]
    means = shifts.mean(axis=0)
    stds = shifts.std(axis=0)
    return (means,stds)

def thread_call(a):
    return shifts_from_file(db,a[0],a[2],thresholds,args.verbose)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d',"--database",help="database file or directory containing .tri files",required=True)
    parser.add_argument('-i','--inputdir',help="directory containing dump files",required=True)
    parser.add_argument('--out',help='output file (default stdout)',type=argparse.FileType('w'),default=sys.stdout)
    parser.add_argument('--resmap',help='file containing mapping from resid to resname',default='')
    parser.add_argument('-v','--verbose',help='output informative messages',action='store_true')
    parser.add_argument('--cpus',help='number of cores to use',default=multiprocessing.cpu_count())
    shiftres.set_default_threshold_args(parser)
        
    args = parser.parse_args()

    thresholds = {'N': (args.Nmin,args.Nmax),
                  'H': (args.Hmin,args.Hmax),
                  'C': (args.Cmin,args.Cmax)}
    
    db = shiftres.read_db(args.database)
    
    #read in residue mapping if present
    resmap = dict()
    if args.resmap:
        for line in open(args.resmap):
            (i,resn) = line.split()
            resmap[int(i)] = resn
        
    files = glob.glob('%s/*.N.*'%args.inputdir)
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
        
    #sort by resid
    inputs.sort(key=lambda (fname,i,resname):i)    
    
    #setup pool
    pool = Pool(args.cpus)
    
    results = pool.map(thread_call,inputs,1)

    args.out.write('ID\tName\tN\tH\tC\tNstd\tHstd\tCstd\n')    
    for ((fname,i,resname),(means,stds)) in zip(inputs,results):
        args.out.write('%d\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n' % (i,resname,means[0],means[1],means[2],stds[0],stds[1],stds[2]))
