# cragr

cragr is an R package for CRAG (**C**ell f**R**ee dn**A** fra**G**mentation) analysis. Specifically, cragr processes cfDNA whole genome sequencing (WGS) data by first obtaining integrated fragmentation scores (IFS) and then identifying hotspots.

## Citation 

Cite our paper:

Zhou, X., & Liu, Y. (2020). De novo characterization of cell-free DNA fragmentation hotspots boosts the power for early detection and localization of multi-cancer. _bioRxiv_.

## Installation

System requirements:

* R 4.x
* tabix and bgzip (included in htslib: [http://www.htslib.org/download/](http://www.htslib.org/download/))
* bedtools ([https://bedtools.readthedocs.io/en/latest/](https://bedtools.readthedocs.io/en/latest/))

Recommended installation steps:

1. Install R devtools: `install.packages(“devtools”)`
2. Install cragr: `devtools::install_github("epifluidlab/cragr")`

## Getting started

cragr provides both R scripts and R functions for the analysis. The quickest way to start is using the R script for your dataset.

The R script can be found at:

	system.file("extdata/scripts", "cragr.R", package = "cragr")

Alternatively, it can be found at the source code directory `inst/ext/scripts/`.

### Input

The workflow requires bgzipped and indexed *fragment data files* as input. A fragment data file is essentially a 0-based BED file ([https://genome.ucsc.edu/FAQ/FAQformat.html](https://genome.ucsc.edu/FAQ/FAQformat.html)), with each row representing a cfDNA fragment. For example:

	14      19000035        19000198        .       27      +
	14      19000044        19000202        .       42      +
	14      19000045        19000202        .       20      -
	14      19000049        19000202        .       12      +

The columns are: chromosome name, fragment start, fragment end, feature name (which can be a dot `.` because cragr doesn’t use this information), mapping quality score (for cfDNA fragments derived from paired-end sequencing, the score is the smaller MAPQ of the two paired reads), and strand (can be omitted because CRAG analysis isn’t strand-aware).

You can prepare the fragment data file from a query-sorted BAM file, or refer to our FinaleDB database for more details:

Zheng, H., Zhu, M. S., & Liu, Y. (2021). FinaleDB: a browser and database of cell-free DNA fragmentation patterns. _Bioinformatics_, 37(16), 2502-2503.

### Start analysis

The workflow are divided into two stages:

1.  Read fragment data and calculate raw IFS scores (without GC
	correction).
2.  Perform GC-correction if needed, and then call hotspots.

The rationale behind this 2-stage workflow is this: for deeply-sequenced samples, stage 1 may require much more memory to complete than stage 2. It is often advisable to process the sequencing data for one chromosome at a time, and then merge the raw IFS scores together for stage 2 analysis. If you are working in an HPC environment, this 2-stage style may be beneficial in terms of cost-effectiveness.

#### Stage 1

The following example perform stage 1 analysis for chromosome 14:

	Rscript cragr.R ifs \
	-i frag.bed.gz \
	-o output.raw_ifs.bed.gz \
	--gc-correct \
	--genome GRCh37 \
	--exclude-region encode.blacklist.hs37-1kg \
	--chrom 14

#### Stage 2

The following example perform stage 2 analysis for chromosome 14:

	Rscript cragr.R peak \
	-i output.raw_ifs.bed.gz \
	-o ifs.bedGraph.gz \
	--gc-correct \
	--genome GRCh37 \
	--exclude-region encode.blacklist.hs37-1kg \
	--chrom 14

The output file `ifs.bedGraph.gz` is the track of IFS scores.

Finally, the hotspot regions can be called from IFS scores:

	zcat ifs.bedGraph.gz |
	bioawk -t 'substr($1,1,1)!="#" && $16<=0.2 && $16!="."' |
	bedtools slop -g {input.chrom_sizes} -i - -b 90 -header |
	bedtools merge -header -i - -d 200 |
	bgzip >> hotspot.bed.gz

### Run with snakemake

For large-scale applications, we suggest use workflow management tools such as snakemake, nextflow, etc. We have included a snakemake script in cragr, located at `system.file("extdata/scripts", "cragr.smk", package = "cragr")`. We encourage the users to check it out for more details.

## License

For academic research, please refer to MIT license. For commerical usage, please contact the authors.