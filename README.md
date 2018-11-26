# Analysis of AMR- and virulence genes, plasmid typing and multilocus sequence typing

This script takes reports from the software ARIBA 
(https://github.com/sanger-pathogens/ariba) and creates summary reports, 
basic statistics reports and visualizations based on user-defined 
settings.

# Usage
```
Rscript amr_vir_sequence_typing.R [options] -o output_folder
```

For help:

```
Rscript amr_vir_sequence_typing.R -h
```
Help screen:

```
Usage: amr_vir_sequence_typing.R [options] -o output_folder


Options:
        -h, --help
                Show this help message and exit

        -u MUT, --mut=MUT
                Location of intrinsic gene reports.

        -a ACQ, --acq=ACQ
                Location of acquired gene reports.

        -i INTRINSIC, --intrinsic=INTRINSIC
                List of intrinsic genes of interest.
                     Type 'all' for including all reported genes.

        -c ACQUIRED, --acquired=ACQUIRED
                List of acquired genes of interest.
                     Type 'all' for including all reported genes.

        -v VIR, --vir=VIR
                Location of ARIBA virulence reports.

        -m MLST, --mlst=MLST
                Location of ARIBA MLST reports.

        -p PLASMID, --plasmid=PLASMID
                Location of ARIBA plasmid reports.

        -o OUTPUT, --output=OUTPUT
                Output directory location.
                     One folder for each analysis will be created
                     at given location

```

## Tracks
The user can specify five different tracks, depending on what they may 
need. The following tracks are available:

- Intrinsic AMR gene analysis (-u, genes: -i)
	+ This track analyses reports from the MEGAres database, and 
gives reports based on which genes are specified by the user in -i.

- Acquired AMR gene analysis (-a, genes: -c)
	+ This track analyses reports from the ResFinder database, and 
gives reports based on which genes are specified by the user in -c.

- Virulence gene analysis (-v)
	+ This track analyses virulence reports from ARIBA and gives a 
summary report and a detailed report on which virulence genes were 
found.

- Multilocus sequence typing analysis (-m)
	+ This track takes summary reports on MLST from ARIBA and gives 
a summary report on sequence types and alleles, as well as a neighbor 
joining tree based on allele distances.

- Plasmid typing analysis (-p)
	+ This track takes plasmidFinder reports from ARIBA and gives 
summary reports on which plasmid types were identified.
