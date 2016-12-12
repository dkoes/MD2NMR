#!/usr/bin/env python

'''Calculate chemical shifts of a single residue trace'''

import numpy as np
import re, gzip
import argparse, cPickle
import sys, os
import makedb

def read_resdump(fname):
    '''Read a dump file and return (pattern,values) indexed by frame'''
    if fname.endswith('.gz'):
        f = gzip.open(fname)
    else:
        f = open(fname)
    
    ret = dict()
    for line in f:
        vals = line.split('|')
        pattern = vals[0]
        values = np.array(vals[2:],np.float)
        i = int(vals[1])
        ret[i] = (pattern,values)
    return ret

def read_db(dbname):
    '''Read database and return, or create if dbname is directory'''
    if os.path.isdir(dbname):
        return makedb.make(dbname)
    else:
        return cPickle.load(open(dbname))

def compute_shifts(resdb,ndata,odata,thresholds=None,verbose=False):
    '''Use residue database to find closes match and return frame indexed N,H,C shifts along with distances to match'''
    # database is indexed by atom->pattern which returns (tree,shifts)
    ret = []
    for i in sorted(ndata.keys()):
        
        #N shifts
        (pattern,values) = ndata[i]                
        if pattern in resdb['N']:
            (tree,shifts) = resdb['N'][pattern]
            (dist,index) = tree.query(values,k=1,p=1)
            (C,H,N) = shifts[index]
        else:
            if verbose:
                sys.stderr.write("Missing pattern %s at resid %d\n" % (pattern,i))
            continue #skip O
        
        if thresholds:
            if N < thresholds['N'][0]:
                continue
            if N > thresholds['N'][1]:
                continue
            if H < thresholds['H'][0]:
                continue
            if H > thresholds['H'][1]:
                continue
            if C < thresholds['C'][0]:
                continue
            if C > thresholds['C'][1]:
                continue                     
        #O shifts
        (pattern,values) = odata[i]                
        if pattern in resdb['O']:
            (tree,shifts) = resdb['O'][pattern]
            (disto,index) = tree.query(values,k=1,p=1)
            (Co,Ho,No) = shifts[index]
        else:
            (Co,Ho,No) = resdb['defaultO']
            disto = -1.0

        #combine N and O
        ret.append((i,C+Co,H+Ho,N+No,dist,disto))
    return ret

def set_default_threshold_args(parser):
    '''Add default threshold values to the passed argparse parser'''
    parser.add_argument("--Nmin",help="N min for outlier removal",default=75.636)
    parser.add_argument("--Nmax",help="N max for outlier removal",default=161.275)
    parser.add_argument("--Hmin",help="H min for outlier removal",default=19.2456)
    parser.add_argument("--Hmax",help="H max for outlier removal",default=31.3854)    
    parser.add_argument("--Cmin",help="C min for outlier removal",default=102.247)
    parser.add_argument("--Cmax",help="C max for outlier removal",default=146.769)
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d',"--database",help="database file or directory containing .tri files",required=True)
    parser.add_argument('-O',"--oxygen",help="O dump file",required=True)
    parser.add_argument('-N',"--nitrogen",help="N dump file",required=True)
    parser.add_argument('--out',help='output file (default stdout)',type=argparse.FileType('w'),default=sys.stdout)
    parser.add_argument('--resname',help='set residue name',default='')
    parser.add_argument('-v','--verbose',help='output informative messages',action='store_true')
    set_default_threshold_args(parser)
        
    args = parser.parse_args()
    
    #figure out resname - earlier versions of the code didn't include this in the dump file name 
    resname = args.resname
    if not resname:
        m = re.search(r'([A-Za-z]{3})\.[NO]\.',args.nitrogen)
        if not m or m.group(1) == 'UMP':
            sys.stderr.write("Could not figure out residue name, please specify.\n")
            sys.exit(1)
        else:
            resname = m.group(1).upper()
            
    thresholds = {'N': (args.Nmin,args.Nmax),
                  'H': (args.Hmin,args.Hmax),
                  'C': (args.Cmin,args.Cmax)}
    odata = read_resdump(args.oxygen)
    ndata = read_resdump(args.nitrogen)
    db = read_db(args.database)
    shifts = compute_shifts(db[resname],ndata,odata,thresholds,args.verbose)
    args.out.write('Frame#\tN\tH\tC\tNdist\tOdist\n')
    for s in shifts:
        args.out.write('%d\t%.3f\t%.3f\t%.3f\t%.2f\t%.2f\n' % s)
    