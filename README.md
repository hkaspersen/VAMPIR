# VAMPIR - Virulence, Amr, Mlst and Plasmid analysis In R

This script takes reports from the software ARIBA 
(https://github.com/sanger-pathogens/ariba) and creates summary reports, 
basic statistics reports and visualizations based on user-defined 
settings.

Author: Håkon Kaspersen, Norwegian Veterinary Institute

# Usage
Make sure that your .Rprofile file have the correct library path before 
beginning. The .Rprofile file is located in your home folder, and you 
get there by typing cd. Type in nano .Rprofile, and paste in the 
following:

```
path <- "/work/projects/nn9305k/lib/R/"

.libPaths(c(path, .libPaths()))
```

Usage:

```
module load R/3.5.0

Rscript VAMPIR.R [options] -o output_folder
```

For help:

```
Rscript VAMPIR.R -h
```
Help screen:

```
Usage: VAMPIR.R [options] -o output_folder


Options:
        -h, --help
                Show this help message and exit

        -u MUT, --mut=MUT
                Directory of megaRes reports.

        -a ACQ, --acq=ACQ
                Directory of resFinder reports.

        -i INTRINSIC, --intrinsic=INTRINSIC
                List of intrinsic genes of interest, used with -u.
                Type 'all' for including all reported genes.
                Can partially match gene names, f. ex. 'gyr' will match all gyr genes identified.
                Example: -i gyr,par,mar

        -c ACQUIRED, --acquired=ACQUIRED
                List of acquired genes of interest, used with -a.
                Type 'all' for including all reported genes.
                Can partially match gene names, f. ex. 'qnr' will match all qnr genes identified.
                Example: -c blaTEM,oqxAB,qnr

        -v VIR, --vir=VIR
                Directory of ARIBA virulence reports.

        -r VIRGENES, --virgenes=VIRGENES
                Virulence genes of interest, use with -v.
                Type 'all' for including all reported genes.

        -d DATABASE, --database=DATABASE
                Virulence database used: virfinder, vfdb or vfdb_core

        -m MLST, --mlst=MLST
                Directory of ARIBA MLST reports.
        
        -p PLASMID, --plasmid_mob=PLASMID
                Directory of Mob suite plasmid reports.

        -q PLASMID, --plasmid_ariba=PLASMID
                Directory of ARIBA plasmid reports.

        -o OUTPUT, --output=OUTPUT
                Output directory.
                One folder for each analysis will be created
                at given location.

        --version
                Print version info.
```

## Tracks

- **Intrinsic AMR gene analysis** (-u, genes: -i)
	+ This track analyses reports from the MEGAres database, and 
gives reports based on which genes are specified by the user in -i.

- **Acquired AMR gene analysis** (-a, genes: -c)
	+ This track analyses reports from the ResFinder database, and 
gives reports based on which genes are specified by the user in -c.

- **Virulence gene analysis** (-v, database: -d, genes: -r)
	+ This track analyses virulence reports from ARIBA and gives a 
summary report and a detailed report on which virulence genes were 
found. Note that only detailed reports will be generated when
specifying the vfdb or vfdb_core databases.

- **Multilocus sequence typing analysis** (-m)
	+ This track takes summary reports on MLST from ARIBA and gives 
a summary report on sequence types and alleles, as well as a neighbor 
joining tree based on allele distances.

- **Plasmid typing analysis** (-p, -q)
	+ This track takes mob-suite plasmid reports (-p) or plasmidFinder
reports from ARIBA (-q) and gives summary reports on which plasmid
types were identified.

# Output files

- ***_report.txt**: A tab separated text file containing columns with 
genes, and 1/0 for present/absent for each isolate (row).

- ***_flags.txt**: A tab separated text file containing the quality 
control values (flags) for each gene/mutation found in the respective 
report. The column "flag_result" determines in the respective gene/mutation 
passed quality control (1) or not (0). Note that all reported genes and 
mutations are presented here.

- ***_stats.txt**: A tab separated text file containing the summary 
statistics of the respective analysis type. It presents the percentage 
of isolates where the given gene is present, as well as a 95 % 
confidence interval.

- **virulence_detailed_report.txt**: A tab separated text file 
containing the specific sub types of virulence genes, rather than the grouped 
virulence genes in virulence_summary_report.txt.

- **virulence_summary_report.txt**: A tab separated text file containing 
the grouped virulence gene results.

- **intrinsic_mut_report.txt**: A tab separated text file containing the 
detailed mutations found in specified intrinsic genes.

- **mlst_tree.svg**: A figure presenting the phylogenetic tree based on 
allele distances.

- **mlst_tree.newick**: A newick-format phylogenetic tree based on 
allele distances.
