#!/usr/bin/env python

'''Calculate chemical shifts of a single residue trace'''

import numpy as np
import re, gzip
import argparse, cPickle, gzip
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
    elif dbname.endswith('.gz'):
        return cPickle.load(gzip.open(dbname))
    else:
        return cPickle.load(open(dbname))

def compute_shifts(resdb,ndata,odata,args,verbose=False):
    '''Use residue database to find closes match and return frame indexed N,H,C shifts along with distances to match'''
    # database is indexed by atom->pattern which returns (tree,shifts)
    ret = []
    refs = np.array([args.Nref,args.Href,args.Cref])

    for i in sorted(ndata.keys()):
        #N shifts
        (pattern,values) = ndata[i]                
        if pattern in resdb['N']:
            (tree,shifts) = resdb['N'][pattern]
            (dists,indices) = tree.query(values,k=args.Nnum_matches,p=1,distance_upper_bound=args.Nmax_match*len(values))
            #annoyingly unpackes singletons
            if type(dists) is float:
                dists = [dists]
                indices = [indices]
            #this is ridiculous, but query will return values even if there isn't anything
            dists = filter(lambda x: np.isfinite(x), dists)
            indices = filter(lambda i: i < len(shifts), indices)
            
            if dists:
                dist = np.mean(dists)/len(values)
                (C,H,N) = np.mean(shifts.reshape(-1,3)[indices,:],axis=0)
            else:
                continue #skip
        else:
            if verbose:
                sys.stderr.write("Missing pattern %s at frame %d\n" % (pattern,i))
            continue #skip O
                         
        #O shifts
        (pattern,values) = odata[i]                
        if pattern in resdb['O']:
            (tree,shifts) = resdb['O'][pattern]
            
            (dists,indices) = tree.query(values,k=args.Onum_matches,p=1,distance_upper_bound=args.Omax_match*len(values))
            #annoyingly unpackes singletons
            if type(dists) is float:
                dists = [dists]
                indices = [indices]
            #this is ridiculous, but query will return values even if there isn't anything
            dists = filter(lambda x: np.isfinite(x), dists)
            indices = filter(lambda i: i < len(shifts), indices)
            
            if dists:
                disto = np.mean(dists)/len(values)
                (Co,Ho,No) = np.mean(shifts.reshape(-1,3)[indices,:],axis=0)
            else:
                (Co,Ho,No) = resdb['defaultO'] # use default
                disto = -1.0          
        else:
            (Co,Ho,No) = resdb['defaultO']
            disto = -1.0

        #combine N and O
        mrst = [C+Co,H+Ho,N+No]
        if args.mrst:
            shift = mrst
        else: #correct
            shift = refs-mrst
        ret.append((i,shift[0],shift[1],shift[2],dist,disto))
    return ret

def add_shift_args(parser):
    '''add arguments to parser that effect shift calculation'''
    parser.add_argument('--mrst',help='report MRST values instead of shifts',action='store_true')    
    parser.add_argument('--Nref',help='N reference value',default=229.853)
    parser.add_argument('--Href',help='H reference value',default=31.792)
    parser.add_argument('--Cref',help='C reference value',default=178.947)
    parser.add_argument('--Nmax-match',help='Maximum value of a match for N',type=float,default=float('inf'))
    parser.add_argument('--Nnum-matches',help='Maximum number of templates to match for N',type=int,default=1)
    parser.add_argument('--Omax-match',help='Maximum value of a match for O (default used for non-matching)',type=float,default=float('inf'))
    parser.add_argument('--Onum-matches',help='Maximum number of templates to match for O (default used for not-matching)',type=int,default=1)    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d',"--database",help="database file or directory containing .tri files",required=True)
    parser.add_argument('-O',"--oxygen",help="O dump file",required=True)
    parser.add_argument('-N',"--nitrogen",help="N dump file",required=True)
    parser.add_argument('--out',help='output file (default stdout)',type=argparse.FileType('w'),default=sys.stdout)
    parser.add_argument('--resname',help='set residue name',default='')
    parser.add_argument('-v','--verbose',help='output informative messages',action='store_true')
    add_shift_args(parser)
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

    odata = read_resdump(args.oxygen)
    ndata = read_resdump(args.nitrogen)
    db = read_db(args.database)
    shifts = compute_shifts(db[resname],ndata,odata,args,args.verbose)
    args.out.write('Frame#\tN\tH\tC\tNdist\tOdist\n')
    for s in shifts:
        args.out.write('%d\t%.3f\t%.3f\t%.3f\t%.2f\t%.2f\n' % s)
    
