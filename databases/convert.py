#!/usr/bin/env python2

import sys, gzip, cPickle, scipy

'''Convert python2 dictionaries that have ckdtrees that aren't compatible
with python3/newer scipy to dictionary with raw values'''

db = cPickle.load(gzip.open(sys.argv[1]))
out = open(sys.argv[2],'w')
 
untree = {}
for r in db.iterkeys():
    r = r.decode('utf-8')
    if type(db[r]) != dict:
        untree[r] = db[r]
    else:
        untree[r] = {}
        for a in db[r].iterkeys():
            a = a.decode('utf-8')
            if type(db[r][a]) != dict:
                untree[r][a] = db[r][a]
            else:
                untree[r][a] = {}
                for pat in db[r][a].iterkeys():
                    pat = pat.decode('utf-8')
                    (tree,shifts) = db[r][a][pat]
                    untree[r][a][pat] = (tree.data, shifts)


cPickle.dump(untree,out,cPickle.HIGHEST_PROTOCOL)
