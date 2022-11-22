# cragr

cragr is an R package for CRAG (**C**ell f**R**ee dn**A** fra**G**mentation) analysis. Specifically, cragr processes cfDNA whole genome sequencing (WGS) data by first obtaining integrated fragmentation scores (IFS) and then identifying hotspots.

## Citation

Cite our paper:

Zhou X, Zheng H, Fu H, McKillip KL, Pinney SM, Liu Y. (2022) CRAG: De novo characterization of cell-free DNA fragmentation hotspots in plasma whole-genome sequencing. Genome Medicine in press; preprint doi: https://doi.org/10.1101/2020.07.16.201350


## Installation

System requirements:

- R 4.1.x or higher
- tabix and bgzip (included in [htslib](http://www.htslib.org/download/))
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- For Mac OSX, install xcode first

Recommended installation steps:

1. Install R devtools: `install.packages("devtools")`
2. Install cragr: `devtools::install_github("epifluidlab/cragr")`

cragr also requires BSgenome packages for the reference genome of interest, in order to perform GC correction.

* For GRCh37: install [BSgenome.Hsapiens.1000genomes.hs37d5](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.1000genomes.hs37d5.html)
* For GRCh38: install [BSgenome.Hsapiens.NCBI.GRCh38](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.NCBI.GRCh38.html)
* Currently cragr only supports GRCh37 and GRCh38.

## Getting started

cragr provides both R scripts and R functions for the analysis. The quickest way to start is using the R script for your dataset.

The R script can be found at:

    system.file("extdata/scripts", "cragr.R", package = "cragr")

Alternatively, it can be found at the source code directory `inst/ext/scripts/`.

### Input

The workflow requires bgzipped and indexed _fragment data files_ as input. A fragment data file is essentially a 0-based BED file ([https://genome.ucsc.edu/FAQ/FAQformat.html](https://genome.ucsc.edu/FAQ/FAQformat.html)), with each row representing a cfDNA fragment. For example:

    14      19000035        19000198        .       27      +
    14      19000044        19000202        .       42      +
    14      19000045        19000202        .       20      -
    14      19000049        19000202        .       12      +

The columns are: chromosome name, fragment start, fragment end, feature name (which can be a dot `.` because cragr doesn’t use this information), mapping quality score (for cfDNA fragments derived from paired-end sequencing, the score is the smaller MAPQ of the two paired reads), and strand (can be omitted because CRAG analysis isn’t strand-aware).

You can prepare the fragment data file from a query-sorted BAM file, or refer to our FinaleDB database for more details:

> Zheng, H., Zhu, M. S., & Liu, Y. (2021). FinaleDB: a browser and database of cell-free DNA fragmentation patterns. _Bioinformatics_, 37(16), 2502-2503.

#### Example datasets

We provide two example datasets for testing the package.

The example shipped with the package (`inst/extdata/example_frag.bed.gz`) is based on shallow whole-genome sequencing result for a cfDNA plasma sample, that only fragments from chr14 and chr21 are included. The dataset is quite small so that users may quickly test the workflow without consuming much computing resources. However, because of the low coverage, while IFS scores can be calculated without problems, in the end hotspots can't be identified.

We also provide another example which is based on the BH01 dataset (Snyder 2016), a deeply-sequenced cfDNA plasma sample. The fragment file can be downloaded at https://zenodo.org/record/5227520.

For you information, below are steps from the BAM file to the fragment file:

```
> samtools sort -n -o qsorted.bam -@ 4 sample.bam
> samtools view -h -f 3 -F 3852 -G 48 --incl-flags 48 \
  qsorted.bam |
  bamToBed -bedpe -mate1 -i stdin |
  awk -F'\t' -v OFS="\t" '{if ($1!=$4) next; if ($9=="+") {s=$2;e=$6} else {s=$5;e=$3} if (e>s) print $1,s,e,$7,$8,$9}' |
  sort -k1,1V -k2,2n |
  bgzip > sample.frag.bed.gz
> tabix -p bed sample.frag.bed.gz
```

#### Additional data files

In addition to the input fragment files, cragr needs several additional data files as well.

* Blacklist regions: it is recommended to exclude certain problematic regions from the analysis, due to various reasons such as low mappability. Please refer to [ENCODE Blacklist](https://www.nature.com/articles/s41598-019-45839-z) for more details. The region file should be in BED format. See this [example](https://github.com/epifluidlab/cragr/blob/3d419a49/inst/extdata/wgEncodeDacMapabilityConsensusExcludable.hs37-1kg.bed).
* Chromosome sizes: can be obtained from the FASTA file of the reference genome: `samtools faidx input.fa`. See this [example](https://github.com/epifluidlab/cragr/blob/3d419a49/inst/extdata/human_g1k_v37.chrom.sizes).
* High-mappability file: a BED file containing regions with high mappability scores.

### Start analysis

The analysis consists of several stages:

1. Read fragment data and calculate raw IFS scores (those without GC correction).
2. Perform GC-correction if needed, and calcualte p-values and FDR values for each genomic interval based on negative binomial model.
3. Call hotspots.
4. (Optional) given a cfDNA fragment dataset, calculate IFS scores for certain genomic intervals (for example, the hotspots obtained from Stage 3).

#### Stage 1

The following example perform stage 1 analysis for chromosome 14:

    Rscript cragr.R ifs \
    -i frag.bed.gz \
    -o output.raw_ifs.chr14.bed.gz \
    --gc-correct \
    --genome GRCh37 \
    --exclude-region encode.blacklist.hs37-1kg.bed \
    --chrom 14

It is recommended to calculate IFS scores chromosome by chromosome because this enables the most potential for parallelization. In this case, after all tasks are finished, concatenate the raw IFS tracks together.

#### Stage 2

The following example perform stage 2 analysis for chromosome 14:

    Rscript cragr.R peak \
    -i output.raw_ifs.bed.gz \
    -o ifs.bedGraph.gz \
    --gc-correct \
    --genome GRCh37

The output file `ifs.bedGraph.gz` is the track of IFS scores. It can be visualized using standard genome browser such as IGV.

#### Stage 3

Finally, the hotspot regions can be called from IFS scores:

    zcat ifs.bedGraph.gz |
    awk -F'\t' -v OFS="\t" 'substr($1,1,1)!="#" && $17<=FDR_CUTOFF && $17!="."' |
    bedtools slop -g CHROM_SIZES -i - -b FLANK -header |
    bedtools merge -header -i - -d MERGE_GAP -c 17 -o min > hotspot.bed

Pay attention to these parameters:

* `FDR_CUTOFF`: in the IFS track, if the FDR value of an interval is smaller than the cutoff, the interval will be considered belonging to a hotspot.
* `CHROM_SIZES`: path to the chromosome sizes file.
* `FLANK`: in Stage 1 analysis, there are two important parameters: `--window-size` (200 by default) and `--step-size` (20 by default). `FLANK` must be `(WINDOW_SIZE - STEP_SIZE) / 2`. Therefore, by default, `FLANK = (200 - 20) / 2 = 90`.
* `MERGE_GAP`: if the distance between two hotspot intervals are less than `MERGE_GAP`, they will be merged into one larger hotspot.


The output hotspot file contains fourth columns as shown below:

```
22      16077180        16077500        0.04534177218
22      16078380        16078920        0.05546186323
22      16079260        16079740        0.05591242971
22      16080180        16080380        0.07216067306
```

The fourth column is the FDR (BH-corrected) associated with the hotspot. Lower values indicate more significant hotspots.

### Stage 4

After obtaining hotspots, one can calculate IFS scores over the hotspots for any cfDNA samples:

```
Rscript cragr.R signal \
-i frag.bed.gz \
--hotspot hotspot.bed \
-o ifs_hotspot.bed.gz \
--gc-correct \
--genome GRCh37 \
--exclude-region encode.blacklist.hs37-1kg.bed \
--chrom 22
```

## Script usage

### `ifs`

```
Options:
        -i INPUT, --input=INPUT
                Path to the input fragment file. The file should be in bgzip-compressed BED format, alongside with the .tbi index file.

        -o OUTPUT, --output=OUTPUT
                Path to the output file.

        --genome=GENOME
                Which reference genome the input fragment file is based on. Should be either GRCh37 or GRCh38.

        -g, --gc-correct
                Perform GC correction.

        --gc-correct-method=GC-CORRECT-METHOD
                GC correction method. Should be either standard or caret [standard].

        --gc-correct-n=GC-CORRECT-N
                Maximal number of data points for GC correction model training [1000000].

        -m HIGH-MAPPABILITY, --high-mappability=HIGH-MAPPABILITY
                Path to the mappability file, which should be in BED format [NULL].

        --chrom=CHROM
                Perform the analysis for the specified chromosome.

        --min-mapq=MIN-MAPQ
                Minimal MAPQ for fragments included in the analysis.

        --min-fraglen=MIN-FRAGLEN
                Minimal length for fragments included in the analysis.

        --max-fraglen=MAX-FRAGLEN
                Maximal length for fragments included in the analysis.

        --exclude-region=EXCLUDE-REGION
                BED files defining regions to be excluded from the analysis, separated by colon. Default is the ENCODE Blacklist: https://www.nature.com/articles/s41598-019-45839-z, which is included in this R package

        -w WINDOW-SIZE, --window-size=WINDOW-SIZE
                Size of the sliding window [200]

        -s STEP-SIZE, --step-size=STEP-SIZE
                Step size of the sliding window [20]

        -t THREAD, --thread=THREAD
```

### `peak`

```
Options:
        -i INPUT, --input=INPUT
                Path to the input file. If there are multiple input files, they should be separated by colons

        -o OUTPUT, --output=OUTPUT
                Path to output file

        --genome=GENOME
                Genome of the input

        -g, --gc-correct
                Whether to perform GC correction

        --gc-correct-method=GC-CORRECT-METHOD
                Methods used in GC correction. Should be either standard or caret [standard]

        --gc-correct-n=GC-CORRECT-N
                Maximal sample size for GC correction model training [1e6L]

        -w WINDOW-SIZE, --window-size=WINDOW-SIZE
                Size of the sliding window [200]

        -s STEP-SIZE, --step-size=STEP-SIZE
                Step size of the sliding window [20]

        -t THREAD, --thread=THREAD
```

### `signal`

```
Options:
        -i INPUT, --input=INPUT
                Path to the input file. If there are multiple input files, they should be separated by colons

        --hotspot=HOTSPOT
                Path to the hotspot file

        -o OUTPUT, --output=OUTPUT
                Path to output file

        --genome=GENOME
                Reference genome to use. Should be either GRCh37 or GRCh38.

        -g, --gc-correct
                Whether to perform GC correction.

        --gc-correct-method=GC-CORRECT-METHOD
                Methods used in GC correction. Should be either standard or caret [standard]

        --gc-correct-n=GC-CORRECT-N
                Maximal sample size for GC correction model training [1e6L]

        -m HIGH-MAPPABILITY, --high-mappability=HIGH-MAPPABILITY
                Path to the mappability file. Default is NULL, i.e. do NOT exclude fragments from low-mappability regions

        --chrom=CHROM
                Perform the analysis only for a selected group of chromosomes. Separated by colons, such as 12:16:X. If not provided, all chromosomes found in the input file will be used

        --min-mapq=MIN-MAPQ
                Minimal MAPQ for fragments included in the analysis

        --min-fraglen=MIN-FRAGLEN
                Minimal length for fragments included in the analysis

        --max-fraglen=MAX-FRAGLEN
                Maximal length for fragments included in the analysis

        --exclude-region=EXCLUDE-REGION
                BED files defining regions to be excluded from the analysis, separated by colon. Default is the ENCODE Blacklist: https://www.nature.com/articles/s41598-019-45839-z, which is included in this R package

        -w WINDOW-SIZE, --window-size=WINDOW-SIZE
                Size of the sliding window [200]

        -t THREAD, --thread=THREAD
```

## Run with snakemake

For large-scale applications, we suggest use workflow management tools such as snakemake, nextflow, etc. We have included a snakemake script in cragr, located at `system.file("extdata/scripts", "cragr.smk", package = "cragr")`. We encourage the users to check it out for more details.

## License

For academic research, please refer to MIT license. For commerical usage, please contact the authors.
