# Databases of ROI templates with pre-calculated chemical shifts.

This directory contains databases for use with MD2NMR scripts.  Each database includes the distance descriptors and calculated values for residue regions of interest.

The newest database should always be used.  Older databases are provided for posterity. 

1. **db_12_2016.db**  The database used in [our first methods paper](http://www.sciencedirect.com/science/article/pii/S2210271X16304790).  This contains sufficient templates to provide high coverage for six reference structures (1ENH, 1IGD, 1HIK, 3OBL, 1QZM, 1UBQ) for a significant amount of simulation.
2. **db_2_2017.db** The database used to generate figures for a poster for the 2017 Biophysical Society Meeting.  This addes templates to provide high coverage for 2A3D and 1L2Y.
3. **db_2_2017_py3.db** Same as above, but for usage with Python3 code.  No longer pickle scipy.spatial.cKDTree.
