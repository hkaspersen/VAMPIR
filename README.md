# Analysis of AMR- and virulence genes, plasmid typing and multilocus sequence typing

This script takes reports from the software ARIBA 
(https://github.com/sanger-pathogens/ariba) and creates summary reports, 
basic statistics reports and visualizations based on user-defined 
settings.

Author: Håkon Kaspersen, Norwegian Veterinary Institute

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

# Output file descriptions

- *_report.txt: A tab separated text file containing columns with genes, 
and 1/0 for present/absent for each isolate (row).

- *_flags.txt: A tab separated text file containing the quality control 
values (flags) for each gene/mutation found in the respective report. 
The column "flag_result" determines in the respective gene/mutation 
passed quality control (1) or not (0). Note that all reported genes and 
mutations are presented here.

- *_stats.txt: A tab separated text file containing the summary 
statistics of the respective analysis type. It presents the percentage 
of isolates where the given gene is present, as well as a 95 % 
confidence interval.

- virulence_detailed_report.txt: A tab separated text file containing 
the specific subtypes of virulence genes, rather than the grouped 
virulence genes in virulence_summary_report.txt.

- virulence_summary_report.txt: A tab separated text file containing the 
grouped virulence gene results.

- intrinsic_mut_report.txt: A tab separated text file containing the 
detailed mutations found in specified intrinsic genes.

- mlst_tree.svg: A figure presenting the phylogenetic tree based on 
allele distances.

- mlst_tree.newick: A newick-format phylogenetic tree based on allele 
distances.
