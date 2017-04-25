#!/usr/bin/env python

'''Create python database from tri files'''

import numpy as np
import glob, re, collections
import argparse, cPickle
from scipy import spatial

def make(tridir,keepoutliers=False,dynamicfilter=False):
    '''Create an optimized python database for looking up calculated chemical shifts by pattern'''
    db = dict()
    
    #scan through all files to compute overall mean - this could be more efficient
    ntotal = []
    ototal = []
    alltotal = []
    for tri in glob.glob(tridir+'/*.tri'):
        m = re.search(r'%s/*(.*)\.([NO])\.tri'%tridir.rstrip('/'),tri)
        if m:
            atom = m.group(2)            
            for line in open(tri):
                vals = line.split('|')
                nhc = np.array(vals[1:4],np.float)
                alltotal.append(nhc)
                if atom == 'N':
                    ntotal.append(nhc)
                elif atom == 'O':
                    ototal.append(nhc)
                
    nmean = np.mean(ntotal,axis=0)
    nstd = np.std(ntotal,axis=0)
    omean = np.mean(ototal,axis=0)
    ostd = np.std(ototal,axis=0)
    allmean = np.mean(alltotal, axis=0)
    allstd = np.std(alltotal, axis=0)
                
    db['defaultO'] = omean
    for tri in glob.glob(tridir+'/*.tri'):
        m = re.search(r'%s/*(.*)\.([NO])\.tri'%tridir.rstrip('/'),tri)
        if m:
            res = m.group(1)
            atom = m.group(2)
            if res not in db: db[res] = dict()
            if atom not in db[res]: db[res][atom] = dict()
            
            allres = []
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
                allres.append(np.array([N,H,C]))           
                    
            #convert data for each pattern into a kdtree and the corresponding shifts
            for (pattern,vals) in valuesbypattern.iteritems():    
                shifts = shiftsbypattern[pattern]         
                if not keepoutliers:
                    if dynamicfilter:
                        newshifts = []
                        newvals = []

                        for (nhc,v) in zip(shifts,vals):
                            if atom == 'N':
                                mean = nmean
                                std = nstd
                            else:
                                mean = omean
                                std = ostd
                            for i in xrange(3):
                                if nhc[i] < mean[i]-std[i]*3 or nhc[i] > mean[i]+std[i]*3:
                                    break
                            else: #executed if we _didn't_ break out of the leep 
                                newvals.append(v)
                                newshifts.append(nhc)
                        vals = newvals
                        shifts = newshifts
                    else: #not dynamic                       
                        NMIN = 77.068
                        NMAX = 160.337
                        HMIN = 19.566
                        HMAX = 31.142
                        CMIN = 102.217
                        CMAX = 146.892
      
                        newshifts = []
                        newvals = []                        
                        for (nhc,v) in zip(shifts,vals):
                            if nhc[0] >= NMIN and nhc[0] <= NMAX and \
                               nhc[1] >= HMIN and nhc[1] <= HMAX and \
                               nhc[2] >= CMIN and nhc[2] <= CMAX:
                                newvals.append(v)
                                newshifts.append(nhc)
                        vals = newvals
                        shifts = newshifts
                
                if len(vals):
                    tree = spatial.cKDTree(np.array(vals),copy_data=True)
                    db[res][atom][pattern] = (tree,np.array(shifts))
            if atom == 'O':
                db[res]['defaultO'] = np.mean(allres,axis=0)
                
    return db
                

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d',"--tridir",help="directory containing .tri files",required=True)
    parser.add_argument('-o',"--output",help="database output file",default="nmr.db")
    parser.add_argument('--disable-filter',help='do not remove outliers (3 sigma)',action='store_true')
    parser.add_argument('--dynamic-filter',help='compute per-residue std and mean for filter, rather than using predefined constants',action='store_true')
    args = parser.parse_args()
    out = open(args.output,'w')
    db = make(args.tridir,args.disable_filter,args.dynamic_filter)
    #histidine variants
    db['HIS'] = db['HSE']
    db['HIE'] = db['HSE']
    db['HSD'] = db['HSE']
    #cysteine variants
    db['CYX'] = db['CYS']
    cPickle.dump(db,out,cPickle.HIGHEST_PROTOCOL)
