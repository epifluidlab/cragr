
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cragr

<!-- badges: start -->
<!-- badges: end -->

cragr is an R package for CRAG analysis (**C**ell f**R**ee dn**A**
fra**G**mentation (CRAG): De novo characterization of cell-free DNA
fragmentation hotspots).

For more details of CRAG, refer to our bioRxiv preprint (Zhou and Liu
2020).

## Installation

You can install the development version from
[GitHub](https://github.com/epifluidlab/cragr) with:

``` r
# install.packages("devtools")
devtools::install_github("haizi-zh/cragr")
```

In addition, usually we need to perform GC-correction for IFS scores.
Therefore, you need to install corresponding `BSGenome` packages.
Currently, we only support analysis for
[hs37-1kg](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use).
As a result, you need to install it from
[Bioconductor](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.1000genomes.hs37d5.html):

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")
```

## Getting started

cragr provides both R APIs and an R script for the analysis. The
quickest way to start is using the R script for your dataset.

You can find the script in R package installation:

``` r
system.file("extdata/scripts", "cragr.R", package = "cragr")
```

Or simply clone the source code from GitHub and find the script in the
directory `inst/ext/scripts/`.

### Preparation

The pipeline requires *fragment data files* as input. A fragment data
file is essentially a BED file, with each row representing a cfDNA
fragment. Below is an example:

    14      19000035        19000198        .       27      +
    14      19000044        19000202        .       42      +
    14      19000045        19000202        .       20      -
    14      19000049        19000202        .       12      +

The six columns are: chrom, start, end, feature name (which we do not
use), MAPQ (which is the smallest MAPQ of the two paired reads), and
strand (which we do not use).

You can easily prepare the fragment data file from a query-sorted BAM
file, or refer to our FinaleDB paper for more details (Zheng, Zhu, and
Liu 2020).

### Start analysis

The pipeline are divided into two stages:

1.  Read fragment data and calculate raw IFS scores (without GC
    correction).
2.  Perform GC-correction if needed, and then call hotspots.

The rationale behind this 2-stage workflow is this: for deeply-sequenced
samples, stage 1 may require much more memory to complete than stage 2.
If you are working in an HPC environment, this 2-stage design may be
beneficial in terms of cost-effectiveness.

Stage 1

The following example perform stage 1 analysis for chromosome 14:

    Rscript cragr.R stage1 \
    -i frag.bed.gz \
    --output-ifs ifs_raw.bed.gz \
    -m high_mappability.bed.gz \
    --chrom 14

Here, `-m` is a pre-calculated BED file defining high-mappability
regions. Each interval must **match** window size and step size used in
the analysis. For example, the CRAG preprint suggests choose sliding
windows of 200 bp, at a step size of 20 bp. In this case, the
high-mappability file for chr14 should be like:

    14      19002020        19002220        0.90624999
    14      19002040        19002240        0.90624999
    14      19002060        19002260        0.90624999
    14      19003780        19003980        0.9062843997

Stage 2

The following example perform stage 2 analysis for chromosome 14:

    Rscript cragr.R stage2 \
    -i ifs_raw.bed.gz \
    --output-ifs ifs.bedGraph.gz \
    --output-hotspot hotspot.bed.gz \
    --chrom 14

Then you will get two output file:

`ifs.bedGraph.gz` is the track of IFS scores:

    #chrom  start   end     score   cov     gc      score0  z_score pval    pval_adjust     pval_local
    14      20195550        20195570        155.70130800196 82      0.35    163.646238889208        0.272700302902589       0.69838482119691        1       0.222145318665377
    14      20195570        20195590        164.705321220779        86      0.37    171.8828729757  0.713220028696098       0.892403670170209       1       0.477082038689824

`hotspot.bed.gz` contains all called hotspots. Adjacent and overlapping
hotspots are merged togeter:

    #chrom  start   end     z_score score   pval    pval_adjust     pval_local      cov     score0
    14      20254300        20254540        -3.16060640969945       85.526155302845 1.86548322442011e-07    9.65582497336755e-06    1.44516815909016e-08    44      89.5529678676775
    14      20255200        20255520        -3.15469665138396       85.6469479458914        4.97775337339482e-06    0.000164618460431793    3.43963679886297e-08    46.5    91.165760047382
    14      20330280        20330640        -2.97445785387944       89.3309430549676        1.29560333360055e-06    5.19151009461698e-05    1.40943389645918e-08    45.25   90.3545474112741

### Run with snakemake

For large-scale real applications, we suggest use workflow management
tools such as snakemake, nextflow, etc. We provide a snakemake file.
Similar to `cragr.R`, you can find it:

In R package installation:

``` r
system.file("extdata/scripts", "cragr.smk", package = "cragr")
```

Or simply clone the source code from GitHub and find the script in the
directory `inst/ext/scripts/`.

To run the snakemake pipeline, you need to prepare the workspace
directory as below:

1.  Fragment data files are placed in directory `frag`, and named as
    \`{sample\_id}.frag.bed.gz"
2.  High-mappability data file is placed in directory `data` and named
    as `mappability.hs37-1kg.w200.s20.0_9.bed.gz`

Then run start the snakemake pipeline:

    snakemake --cores {# of cpu cores} -s cragr.smk \
    --config SCRIPT_HOME={directory containing your cragr.R} \
    -n {sample_id}.ifs.bedGraph.gz

Note: usually snakemake workflow is more related to your specific
computing environment. As a result, the example above is only for your
reference. Have a look in both the snakemake file and the command above,
and modify it accordingly.

## Citation

*De novo* characterization of cell-free DNA fragmentation hotspots
boosts the power for early detection and localization of multi-cancer  
Xionghui Zhou, Yaping Liu  
bioRxiv 2020.07.16.201350; doi:
<https://doi.org/10.1101/2020.07.16.201350>

## License

For academic research, please refer to MIT license. For commerical
usage, please contact the authors.

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-zheng2020" class="csl-entry">

Zheng, Haizi, Michelle S Zhu, and Yaping Liu. 2020. ???FinaleDB: A Browser
and Database of Cell-Free DNA Fragmentation Patterns.??? Edited by
Robinson Peter. *Bioinformatics*, December, btaa999.
<https://doi.org/10.1093/bioinformatics/btaa999>.

</div>

<div id="ref-zhou2020" class="csl-entry">

Zhou, Xionghui, and Yaping Liu. 2020. ???De Novo Characterization of
Cell-Free DNA Fragmentation Hotspots Boosts the Power for Early
Detection and Localization of Multi-Cancer.??? *bioRxiv*, July.
<https://doi.org/10.1101/2020.07.16.201350>.

</div>

</div>
