# MD2NMR

This project provides an *ab initio*, template-based approach for directly calculating chemical shifts from molecular dynamics simulations.  The calculated shifts are in no way empirically fit and are meant to represent the actual chemical environment induced in the simulation.  They can be used to compare simulations to experimental data and as indicator values for the local chemical environment of backbone atoms.


## Usage

```
$ shifts.py --topo 1enh_10ns_0.prmtop --traj 1enh_10ns_0.dcd  -d ../nmr.db
#ID	Name	N	H	C	Nstd	Hstd	Cstd	Coverage
3	ARG	114.219	6.916	53.369	10.642	1.058	4.698	1.000
4	THR	110.062	6.463	59.614	12.136	1.083	5.120	0.999
5	ALA	121.457	6.503	52.190	9.067	0.948	4.606	1.000
6	PHE	117.500	6.479	55.289	8.599	1.085	5.014	0.999
7	SER	113.636	7.643	55.922	8.497	0.958	4.498	1.000
...
```

A topology and trajectory file are provided as input (any format supported by [MDAnalysis](http://www.mdanalysis.org/)) along with a template database of pre-computed shifts (included in this repository).

The output includes the residue number, the residue name, the N, H, and C chemical shifts for the residue with their standard deviations across the trajectory, and the coverage.  The coverage is the percent of frames where a suitable matching template was found for the N region of the residue (the O region has a much smaller effect on the environment and a default average value is substituted when a match cannot be found).

## Options

```
$ ~/git/MD2NMR/shifts.py --help
usage: shifts.py [-h] -d DATABASE [-i INPUTDIR] [--resdir RESDIR]
                 [--topo TOPO] [--traj TRAJ] [--out OUT] [--resmap RESMAP]
                 [-v] [--cpus CPUS] [--mrst] [--Nref NREF] [--Href HREF]
                 [--Cref CREF] [--Nmax-match NMAX_MATCH]
                 [--Nnum-matches NNUM_MATCHES] [--Omax-match OMAX_MATCH]
                 [--Onum-matches ONUM_MATCHES]

optional arguments:
  -h, --help            show this help message and exit
  -d DATABASE, --database DATABASE
                        database file or directory containing .tri files
  -i INPUTDIR, --inputdir INPUTDIR
                        directory containing dump files
  --resdir RESDIR       directory to output per-residue, per-frame shift
                        calculations
  --topo TOPO           topology file of simulation (e.g. prmtop)
  --traj TRAJ           trajectory file of simulation containing explicit
                        wrapped waters (e.g. dcd)
  --out OUT             output file (default stdout)
  --resmap RESMAP       file containing mapping from resid to resname
  -v, --verbose         output informative messages
  --cpus CPUS           number of cores to use
  --mrst                report MRST values instead of shifts
  --Nref NREF           N reference value
  --Href HREF           H reference value
  --Cref CREF           C reference value
  --Nmax-match NMAX_MATCH
                        Maximum value of a match for N
  --Nnum-matches NNUM_MATCHES
                        Maximum number of templates to match for N
  --Omax-match OMAX_MATCH
                        Maximum value of a match for O (default used for non-
                        matching)
  --Onum-matches ONUM_MATCHES
                        Maximum number of templates to match for O (default
                        used for not-matching)
```

## Advanced usage

If you plan to perform multiple analyses or computing shifts directly from the trajectory is too cpu/memory intensive, you can calculate the intermediate distance representation of the trajectory using `dump.py`.

```
$ dump.py --help
usage: dump.py [-h] [-r [RESIDUES [RESIDUES ...]]] [--start START] [--end END]
               [--flush] [--outputdir OUTPUTDIR]
               topology_filename trajectory_filename

positional arguments:
  topology_filename     input file with (.prmtop) extension
  trajectory_filename   input file with (.dcd) extension

optional arguments:
  -h, --help            show this help message and exit
  -r [RESIDUES [RESIDUES ...]], --residues [RESIDUES [RESIDUES ...]]
                        space deliminated residues to compute (all if
                        unspecified)
  --start START         frame to start on
  --end END             frame to end at
  --flush               flush I/O as it is written
  --outputdir OUTPUTDIR
                        directory in which to put DUMP files
```

```
dump.py --topo 1enh_10ns_0.prmtop --traj 1enh_10ns_0.dcd --outputdir dump
shifts.py -d db_12_2016.db -i dump --resdir res1 
shifts.py -d db_12_2016.db -i dump --resdir res2 --Nmax-match .25 --Omax-match .25
```


