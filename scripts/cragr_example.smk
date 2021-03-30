# cragr pipeline

localrules: all, merge

# >>> Configuration >>>
FULL_CORES = config.get("FULL_CORES", 16)
PART_CORES = config.get("PART_CORES", 4)
MEM_PER_CORE = config.get("MEM_PER_CORE", 1800)
WALL_TIME_MAX = config.get("WALL_TIME_MAX", 2880)
CORE_FACTOR = float(config.get("CORE_FACTOR", 1))

# Default base directory for data files. Default: ./data
DATA_DIR = config.get("DATA_DIR", os.path.abspath("data"))


def determine_ifs_cores(chrom, factor=1):
    if chrom == "X":
        cores = 6
    elif chrom == "Y":
        cores = 3
    else:
        chrom = int(chrom)
        if chrom >= 1 and chrom <= 7:
            cores = 6
        elif chrom >= 8 and chrom <= 15:
            cores = 4
        elif chrom >= 16 and chrom <= 22:
            cores = 3
        else:
            raise ValueError

    return round(factor * cores)


rule ifs:
    input: 
        frag="frag/{sample}.hg19.frag.bed.gz",
        gc="data/human_g1k_v37.gc20bp.bed.gz",
        mappability="data/wgEncodeDukeMapabilityUniqueness35bp.hs37-1kg.20bp.bedGraph.gz"
    output: 
        ifs=temp("results/{sample}.chr{chrom}.ifs.raw.bedGraph.gz"),
        ifs_idx=temp("results/{sample}.chr{chrom}.ifs.raw.bedGraph.gz.tbi")
    log: "results/{sample}.chr{chrom}.ifs.raw.log"
    params:
        label=lambda wildcards: f"ifs.{wildcards.sample}.chr{wildcards.chrom}",
    threads: lambda wildcards, attempt: round(determine_ifs_cores(wildcards.chrom, CORE_FACTOR) * (1 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=300,
        attempt=lambda wildcards, threads, attempt: attempt
    shell:
        """
        tmpdir=$(mktemp -d)
        Rscript scripts/cragr.R ifs \
        --input {input.frag} --prefix "$tmpdir"/{wildcards.sample}.chr{wildcards.chrom} \
        --gc {input.gc} --mappability {input.mappability} \
        --chrom {wildcards.chrom} \
        --min-fraglen 50 --max-fraglen 500 2>&1 | tee {log}

        output_name={wildcards.sample}.chr{wildcards.chrom}.ifs.raw.bedGraph.gz
        mv "$tmpdir"/"$output_name" {output.ifs}.tmp
        mv "$tmpdir"/"$output_name".tbi {output.ifs_idx}.tmp
        mv {output.ifs}.tmp {output.ifs}
        mv {output.ifs_idx}.tmp {output.ifs_idx}
        """

rule hotspot:
    input: 
        ifs="results/{sample}.chr{chrom}.ifs.raw.bedGraph.gz",
        ifs_idx="results/{sample}.chr{chrom}.ifs.raw.bedGraph.gz.tbi",
        gc="data/human_g1k_v37.gc20bp.bed.gz",
        mappability="data/wgEncodeDukeMapabilityUniqueness35bp.hs37-1kg.20bp.bedGraph.gz"
    output: 
        ifs=temp("results/{sample}.chr{chrom}.ifs.bedGraph.gz"),
        ifs_idx=temp("results/{sample}.chr{chrom}.ifs.bedGraph.gz.tbi"),
        hotspot=temp("results/{sample}.chr{chrom}.hotspot.bed.gz"),
        hotspot_idx=temp("results/{sample}.chr{chrom}.hotspot.bed.gz.tbi"),
    log: "results/{sample}.chr{chrom}.hotspot.log"
    params:
        label=lambda wildcards: f"hotspot.{wildcards.sample}.chr{wildcards.chrom}",
    threads: lambda wildcards, attempt: round(determine_ifs_cores(wildcards.chrom, CORE_FACTOR) * (1 + 0.5 * attempt) * 0.667)
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=300,
        attempt=lambda wildcards, threads, attempt: attempt
    shell:
        """
        tmpdir=$(mktemp -d)
        Rscript scripts/cragr.R hotspot \
        --input {input.ifs} --prefix "$tmpdir"/{wildcards.sample}.chr{wildcards.chrom} \
        --gc {input.gc} --mappability {input.mappability} \
        --chrom {wildcards.chrom} \
        --min-fraglen 50 --max-fraglen 500 2>&1 | tee {log}

        output_name={wildcards.sample}.chr{wildcards.chrom}.ifs.bedGraph.gz
        mv "$tmpdir"/"$output_name" {output.ifs}.tmp
        mv "$tmpdir"/"$output_name".tbi {output.ifs_idx}.tmp
        mv {output.ifs}.tmp {output.ifs}
        mv {output.ifs_idx}.tmp {output.ifs_idx}

        output_name={wildcards.sample}.chr{wildcards.chrom}.hotspot.bed.gz
        mv "$tmpdir"/"$output_name" {output.hotspot}.tmp
        mv "$tmpdir"/"$output_name".tbi {output.hotspot_idx}.tmp
        mv {output.hotspot}.tmp {output.hotspot}
        mv {output.hotspot_idx}.tmp {output.hotspot_idx}
        """


rule merge:
    input: 
        ifs=expand("results/{{sample}}.chr{chrom}.ifs.bedGraph.gz", chrom=list(range(1, 23))),
        ifs_idx=expand("results/{{sample}}.chr{chrom}.ifs.bedGraph.gz.tbi", chrom=list(range(1, 23))),
        hotspot=expand("results/{{sample}}.chr{chrom}.hotspot.bed.gz", chrom=list(range(1, 23))),
        hotspot_idx=expand("results/{{sample}}.chr{chrom}.hotspot.bed.gz.tbi", chrom=list(range(1, 23))),
    output:
        ifs="results/{sample}.ifs.bedGraph.gz",
        ifs_idx="results/{sample}.ifs.bedGraph.gz.tbi",
        hotspot="results/{sample}.hotspot.bed.gz",
        hotspot_idx="results/{sample}.hotspot.bed.gz.tbi",
    shell:
        """
        tmpdir=$(mktemp -d)

        set +o pipefail
        zcat {input.ifs[0]} | head -n 1 | bgzip > $tmpdir/ifs.bedGraph.gz
        for input_ifs in {input.ifs}
        do
            echo Writing $input_ifs ...
            zcat $input_ifs | tail -n +2 | bgzip >> $tmpdir/ifs.bedGraph.gz
        done
        set -o pipefail
        tabix -p bed $tmpdir/ifs.bedGraph.gz

        set +o pipefail
        zcat {input.hotspot[0]} | head -n 1 | bgzip > $tmpdir/hotspot.bed.gz
        for input_hotspot in {input.hotspot}
        do
            echo Writing $input_hotspot ...
            zcat $input_hotspot | tail -n +2 | bgzip >> $tmpdir/hotspot.bed.gz
        done
        set -o pipefail
        tabix -p bed $tmpdir/hotspot.bed.gz

        mv "$tmpdir"/ifs.bedGraph.gz {output.ifs}.tmp
        mv "$tmpdir"/ifs.bedGraph.gz.tbi {output.ifs_idx}.tmp
        mv {output.ifs}.tmp {output.ifs}
        mv {output.ifs_idx}.tmp {output.ifs_idx}

        mv "$tmpdir"/hotspot.bed.gz {output.hotspot}.tmp
        mv "$tmpdir"/hotspot.bed.gz.tbi {output.hotspot_idx}.tmp
        mv {output.hotspot}.tmp {output.hotspot}
        mv {output.hotspot_idx}.tmp {output.hotspot_idx}
        """

sample_ids, = glob_wildcards("frag/{sample,Pilot2.+}.hg19.frag.bed.gz")
# Head & neck cancer samples
responder_ids = [12, 25, 36, 39, 68, 78, 80, 82]
non_responder_ids = [1, 3, 6, 9, 14, 16, 22, 28, 30 ,34, 42, 45, 63, 65, 71, 75, 86, 89]
sample_ids = responder_ids + non_responder_ids

rule all:
    input: expand("results/Pilot2-{sample}.ifs.bedGraph.gz", sample = sample_ids)