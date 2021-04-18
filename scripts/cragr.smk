# cragr pipeline

localrules: all, merge

# >>> Configuration >>>
FULL_CORES = config.get("FULL_CORES", 16)
PART_CORES = config.get("PART_CORES", 4)
MEM_PER_CORE = config.get("MEM_PER_CORE", 1800)
WALL_TIME_MAX = config.get("WALL_TIME_MAX", 2880)
CORE_FACTOR = float(config.get("CORE_FACTOR", 1))

SCRIPT_HOME = config.get("SCRIPT_HOME")
WALL_TIME_MAX = 2880
MEM_PER_CORE = 1800


def ifs_raw_cores(chrom, file_path):
    import math
    import os.path

    if chrom == "X":
        chrom_idx = 7
    elif chrom == "Y":
        chrom_idx = 20
    else:
        chrom_idx = int(chrom)

    # Factor: if the file size is 4.4G
    # chr1: 8 cores, chr22: 2 cores
    chrom_factor = 8 - 0.286 * (chrom_idx - 1)

    file_size = os.path.getsize(file_path) / pow(1024, 3)
    size_factor = file_size / 4

    return int(max([math.ceil(chrom_factor * size_factor), 2]))


rule raw_ifs:
    input:
        frag="frag/{sid}.frag.bed.gz",
        frag_idx="frag/{sid}.frag.bed.gz.tbi",
        mappability="data/mappability.hs37-1kg.w200.s20.0_9.bed.gz",
        mappability_idx="data/mappability.hs37-1kg.w200.s20.0_9.bed.gz.tbi",
    output: 
        ifs=temp("temp/{sid}.ifs.raw.chr{chrom}.bed.gz")
    log: "temp/{sid}.ifs.raw.chr{chrom}.log"
    params:
        label=lambda wildcards: f"ifs.raw.{wildcards.sid}.chr{wildcards.chrom}",
        script_home=lambda wildcards: SCRIPT_HOME
    threads: lambda wildcards, input, attempt: int(ifs_raw_cores(wildcards.chrom, input.frag) * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=300,
        attempt=lambda wildcards, threads, attempt: attempt
    shell:
        """
        tmpdir=$(mktemp -d)

        Rscript {params.script_home}/cragr.R ifs \
        -i {input.frag} \
        --output-ifs "$tmpdir"/output.bed.gz \
        -m {input.mappability} \
        --chrom {wildcards.chrom} \
        --verbose

        mv "$tmpdir"/output.bed.gz {output.ifs}.tmp
        mv {output.ifs}.tmp {output.ifs}
        """


rule hotspot:
    input: 
        ifs=temp("temp/{sid}.ifs.raw.chr{chrom}.bed.gz"),
        mappability="data/mappability.hs37-1kg.w200.s20.0_9.bed.gz",
    output: 
        ifs="result/{sid}.ifs.chr{chrom}.bedGraph.gz",
        hotspot="result/{sid}.hotspot.chr{chrom}.bed.gz",
    log: "temp/{sid}.ifs.raw.chr{chrom}.log"
    params:
        label=lambda wildcards: f"hotspot.{wildcards.sid}.chr{wildcards.chrom}",
        script_home=lambda wildcards: SCRIPT_HOME
    threads: lambda wildcards, input, attempt: int(2 * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=300,
        attempt=lambda wildcards, threads, attempt: attempt
    shell:
        """
        tmpdir=$(mktemp -d)

        Rscript {params.script_home}/cragr.R hotspot \
        -i {input.ifs} \
        --output-ifs "$tmpdir"/ifs.bedGraph.gz \
        --output-hotspot "$tmpdir"/hotspot.bed.gz \
        -m {input.mappability} \
        --chrom {wildcards.chrom} \
        --verbose

        mv "$tmpdir"/ifs.bedGraph.gz {output.ifs}.tmp
        mv {output.ifs}.tmp {output.ifs}
        mv "$tmpdir"/hotspot.bed.gz {output.hotspot}.tmp
        mv {output.hotspot}.tmp {output.hotspot}
        """