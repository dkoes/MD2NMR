#!/usr/bin/env python

'''Create python database from tri files'''

import numpy as np
import glob, re, collections
import argparse, cPickle
from scipy import spatial

def make(tridir):
    '''Create an optimized python database for looking up calculated chemical shifts by pattern'''
    db = dict()
    for tri in glob.glob(tridir+'/*.tri'):
        m = re.search(r'%s/*(.*)\.([NO])\.tri'%tridir.rstrip('/'),tri)
        defaultO = np.zeros(3)
        if m:
            res = m.group(1)
            atom = m.group(2)
            if res not in db: db[res] = dict()
            if atom not in db[res]: db[res][atom] = dict()
            
            cnt = 0
            valuesbypattern = collections.defaultdict(list)
            shiftsbypattern = collections.defaultdict(list)
            for line in open(tri):
                vals = line.split('|')
                pattern = vals[0]
                N = float(vals[1])
                H = float(vals[2])
                C = float(vals[3])
                coords = np.array(vals[4:],np.float)                
                valuesbypattern[pattern].append(coords)
                shiftsbypattern[pattern].append((N,H,C))
                if atom == 'O':
                    defaultO += np.array([N,H,C])
                    cnt += 1
                    
            #convert data for each pattern into a kdtree and the corresponding shifts
            for (pattern,vals) in valuesbypattern.iteritems():                
                tree = spatial.cKDTree(np.array(vals),copy_data=True)
                shifts = shiftsbypattern[pattern]
                db[res][atom][pattern] = (tree,shifts)
            if atom == 'O':
                db[res]['defaultO'] = defaultO/cnt
                
    return db
                

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d',"--tridir",help="directory containing .tri files")
    parser.add_argument('-o',"--output",help="database output file",default="nmr.db")
    args = parser.parse_args()
    out = open(args.output,'w')
    db = make(args.tridir)
    cPickle.dump(db,out,cPickle.HIGHEST_PROTOCOL)