# MD2NMR

This project provides an *ab initio*, template-based approach for directly calculating chemical shifts from molecular dynamics simulations.  The calculated shifts are in no way empirically fit and are meant to represent the actual chemical environment induced in the simulation.  They can be used to compare simulations to experimental data and as indicator values for the local chemical environment of backbone atoms.

### Update (9/16/2021)

The code has been updated to use Python3 and MDAnalysis 2.0.  The Python3 code uses a different database format (`databases/db_2_2017_py3.db.gz`).
Unfortunately, due to MDAnalysis changes, the code is now approximately 2X slower.  The old Python2 implementation (0.2.0) can be used instead by reverting to MDAnalysis 0.14.0 (e.g. `python2 -m pip install MDAnalysis=0.14.0`).

## Usage

```
$ ./shifts.py --topo example/1enh_10ns_0.prmtop --traj example/small.dcd  -d databases/db_2_2017_py3.db.gz
#ID	Name	N	H	C	Nstd	Hstd	Cstd	Coverage
3	ARG	114.903	6.960	53.272	10.982	1.059	5.197	1.000
4	THR	112.799	6.824	60.177	13.145	1.230	5.310	1.000
5	ALA	120.528	6.752	51.084	8.670	1.007	4.975	1.000
6	PHE	117.034	6.515	55.837	9.036	1.009	5.482	1.000
7	SER	112.478	7.459	57.736	9.023	1.129	4.316	1.000
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


