# cragr pipeline

localrules: all, merge_ifs, call_hotspot_pois, call_hotspot_nb, signal_level

# >>> Configuration >>>
FULL_CORES = config.get("FULL_CORES", 16)
PART_CORES = config.get("PART_CORES", 4)
MEM_PER_CORE = config.get("MEM_PER_CORE", 1800)
WALL_TIME_MAX = config.get("WALL_TIME_MAX", 2880)
CORE_FACTOR = float(config.get("CORE_FACTOR", 1))

WALL_TIME_MAX = 2880
MEM_PER_CORE = 1800


FDR_NB = config.get("FDR_NB", 0.25)
FDR_POIS = config.get("FDR_POIS", 0.01)
GC_CORRECT = config.get("GC_CORRECT", "TRUE")


def find_main_script():
    import subprocess

    return subprocess.run(
        [
            "Rscript",
            "-e",
            'cat(paste0(system.file("extdata/scripts/cragr.R", package = "cragr")))',
        ],
        capture_output=True,
        shell=False,
        text=True,
        check=True,
    ).stdout


def stage_ifs_cores(chrom, file_path, single_chrom=False):
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

    if single_chrom:
        file_size = file_size * (300 / (23 - 0.8 * chrom_idx))

    # 4.4G ~ 1, 1G ~ 0.375
    size_factor = 0.35 + 0.15 * file_size

    return int(max([math.ceil(chrom_factor * size_factor), 2]))


def stage_peak_cores(chrom):
    if chrom == "X":
        chrom_idx = 7
    elif chrom == "Y":
        chrom_idx = 20
    else:
        chrom_idx = int(chrom)

    if chrom_idx <= 6:
        return 4
    elif chrom_idx <= 14:
        return 3
    else:
        return 2


rule stage_ifs:
    input:
        # frag="frag/{sid}.frag.chr{chrom}.bed.gz",
        # frag_idx="frag/{sid}.frag.chr{chrom}.bed.gz.tbi",
        frag="frag/{sid}.frag.bed.gz",
        frag_idx="frag/{sid}.frag.bed.gz.tbi",
        mappability="data/mappability.hs37-1kg.w200.s20.0_9.bed.gz",
        mappability_idx="data/mappability.hs37-1kg.w200.s20.0_9.bed.gz.tbi",
    output:
        ifs=temp("temp/{sid}.ifs.raw.chr{chrom}.bed.gz")
    log: "log/{sid}.chr{chrom}.stage_ifs.log"
    params:
        label=lambda wildcards: f"cragr.stage_ifs.{wildcards.sid}.chr{wildcards.chrom}",
        main_script=lambda wildcards: find_main_script(),
        gc_correct=lambda wildcards: GC_CORRECT
    threads: lambda wildcards, input, attempt: int(stage_ifs_cores(wildcards.chrom, input.frag) * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=300,
        attempt=lambda wildcards, threads, attempt: attempt
    shell:
        """
        tmpdir=$(mktemp -d)

        Rscript {params.main_script} ifs \
        -i {input.frag} \
        -o "$tmpdir"/output.bed.gz \
        --gc-correct {params.gc_correct} \
        -m {input.mappability} \
        --chrom {wildcards.chrom} \
        --verbose 2>&1 | tee {log}

        mv "$tmpdir"/output.bed.gz {output.ifs}.tmp
        mv {output.ifs}.tmp {output.ifs}
        """


rule stage_peak:
    input:
        ifs="temp/{sid}.ifs.raw.chr{chrom}.bed.gz",
        mappability="data/mappability.hs37-1kg.w200.s20.0_9.bed.gz",
    output:
        ifs=temp("temp/{sid}.ifs.chr{chrom}.bedGraph.gz"),
        # hotspot="temp/{sid}.hotspot.chr{chrom}.bed.gz",
    log: "log/{sid}.chr{chrom}.stage_peak.log"
    params:
        label=lambda wildcards: f"cragr.stage_peak.{wildcards.sid}.chr{wildcards.chrom}",
        main_script=lambda wildcards: find_main_script(),
        gc_correct=lambda wildcards: GC_CORRECT
    threads: lambda wildcards, input, attempt: int(stage_peak_cores(wildcards.chrom) * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=300,
        attempt=lambda wildcards, threads, attempt: attempt
    shell:
        """
        tmpdir=$(mktemp -d)

        Rscript {params.main_script} peak \
        -i {input.ifs} \
        -o "$tmpdir"/ifs.bedGraph.gz \
        --gc-correct {params.gc_correct} \
        -m {input.mappability} \
        --chrom {wildcards.chrom} \
        --verbose 2>&1 | tee {log}

        mv "$tmpdir"/ifs.bedGraph.gz {output.ifs}.tmp
        mv {output.ifs}.tmp {output.ifs}
        """


rule merge_ifs:
    input:
        ifs=expand("temp/{{sid}}.ifs.chr{chrom}.bedGraph.gz", chrom=range(1, 23)),
        # hotspot=expand("temp/{{sid}}.hotspot.chr{chrom}.bed.gz", chrom=range(1, 23))
    output:
        ifs="result/{sid}.ifs.bedGraph.gz",
        ifs_idx="result/{sid}.ifs.bedGraph.gz.tbi",
        # hotspot="result/{sid}.hotspot.bed.gz",
        # hotspot_idx="result/{sid}.hotspot.bed.gz.tbi",
    shell:
        """
        tmpdir=$(mktemp -d)
        output_ifs="$tmpdir"/ifs.bedGraph.gz

        echo "Merging IFS data"
        idx=0
        for input_file in {input.ifs}
        do
            echo "Processing $input_file"
            idx=$((idx + 1))
            if [ $idx == 1 ]; then
                zcat < "$input_file" | bgzip > "$output_ifs"
            else
                zcat < "$input_file" | awk 'substr($0,1,1)!="#"' | bgzip >> "$output_ifs"
            fi
        done

        tabix -p bed "$output_ifs"

        mv "$output_ifs" {output.ifs}.tmp
        mv {output.ifs}.tmp {output.ifs}
        mv "$output_ifs".tbi {output.ifs_idx}
        """


# Convert IFS bedGraph files to bigWig files
rule ifs_bigwig:
    input:
        ifs="result/{sid}.ifs.bedGraph.gz",
        ifs_idx="result/{sid}.ifs.bedGraph.gz.tbi",
        chrom_sizes="data/human_g1k_v37.chrom.sizes",
    output: "result/{sid}.ifs.bw"
    shell:
        """
        tmpdir=$(mktemp -d)
        echo "Writing to temporary directory: $tmpdir"

        for chrom in 1 {{10..19}} 2 {{20..22}} {{3..9}}
        do
            echo "Extracting chr$chrom ..."
            tabix {input.ifs} $chrom | bioawk -t '{{print $1,$2,$3,$4}}' >> $tmpdir/ifs.bed
        done
        echo "Converting to bigWig ..."
        bedGraphToBigWig "$tmpdir"/ifs.bed {input.chrom_sizes} $tmpdir/ifs.bw
        unlink "$tmpdir"/ifs.bed

        mv "$tmpdir"/ifs.bw {output}.tmp
        mv {output}.tmp {output}
        """


rule call_hotspot_nb:
    input:
        ifs="result/{sid}.ifs.bedGraph.gz",
        ifs_idx="result/{sid}.ifs.bedGraph.gz.tbi",
    output:
        hotspot="result/{sid}.hotspot.nb.bed.gz",
    log: "log/{sid}.hotspot.nb.log",
    params:
        label=lambda wildcards: f"cragr.hotspot_nb.{wildcards.sid}",
        main_script=lambda wildcards: find_main_script(),
        fdr_nb=lambda wildcards: FDR_NB,
    shell:
        """
        tmpdir=$(mktemp -d)

        Rscript {params.main_script} hotspot \
        -i {input.ifs} \
        -o "$tmpdir"/hotspot.bed.gz \
        --fdr {params.fdr_nb} \
        --hotspot-method nb \
        --verbose 2>&1 | tee {log}

        mv "$tmpdir"/hotspot.bed.gz {output.hotspot}.tmp
        mv {output.hotspot}.tmp {output.hotspot}
        """

rule call_hotspot_pois:
    input:
        ifs="result/{sid}.ifs.bedGraph.gz",
        ifs_idx="result/{sid}.ifs.bedGraph.gz.tbi",
    output:
        hotspot="result/{sid}.hotspot.pois.bed.gz",
    log: "log/{sid}.hotspot.pois.log",
    params:
        label=lambda wildcards: f"cragr.hotspot_pois.{wildcards.sid}",
        main_script=lambda wildcards: find_main_script(),
        fdr_pois=lambda wildcards: FDR_POIS,
    shell:
        """
        tmpdir=$(mktemp -d)

        Rscript {params.main_script} hotspot \
        -i {input.ifs} \
        -o "$tmpdir"/hotspot.bed.gz \
        --fdr {params.fdr_pois} \
        --pval 0.00001 \
        --hotspot-method pois \
        --verbose 2>&1 | tee {log}

        mv "$tmpdir"/hotspot.bed.gz {output.hotspot}.tmp
        mv {output.hotspot}.tmp {output.hotspot}
        """


rule call_hotspot_nb_fdr:
    input:
        ifs="result/{sid}.ifs.bedGraph.gz",
        ifs_idx="result/{sid}.ifs.bedGraph.gz.tbi",
    output:
        hotspot="result/{sid}.hotspot.nb.fdr{fdr}.bed.gz",
    log: "log/{sid}.hotspot.nb.fdr{fdr}.log",
    params:
        label=lambda wildcards: f"cragr.hotspot_nb.{wildcards.sid}",
        main_script=lambda wildcards: find_main_script(),
        fdr_nb=lambda wildcards: FDR_NB,
    shell:
        """
        tmpdir=$(mktemp -d)

        Rscript {params.main_script} hotspot \
        -i {input.ifs} \
        -o "$tmpdir"/hotspot.bed.gz \
        --fdr {wildcards.fdr} \
        --hotspot-method nb \
        --verbose 2>&1 | tee {log}

        mv "$tmpdir"/hotspot.bed.gz {output.hotspot}.tmp
        mv {output.hotspot}.tmp {output.hotspot}
        """


rule call_hotspot_pois_threshold:
    input:
        ifs="result/{sid}.ifs.bedGraph.gz",
        ifs_idx="result/{sid}.ifs.bedGraph.gz.tbi",
    output:
        hotspot="result/{sid}.hotspot.pois.pval{threshold}.bed.gz",
    log: "log/{sid}.hotspot.pois.pval{threshold}.log",
    params:
        label=lambda wildcards: f"cragr.hotspot_nb.{wildcards.sid}",
        main_script=lambda wildcards: find_main_script(),
        fdr_nb=lambda wildcards: FDR_NB,
    shell:
        """
        tmpdir=$(mktemp -d)

        Rscript {params.main_script} hotspot \
        -i {input.ifs} \
        -o "$tmpdir"/hotspot.bed.gz \
        --pval {wildcards.threshold} \
        --hotspot-method pois \
        --verbose 2>&1 | tee {log}

        mv "$tmpdir"/hotspot.bed.gz {output.hotspot}.tmp
        mv {output.hotspot}.tmp {output.hotspot}
        """


rule signal_level:
    input:
        hotspot="result/{sid}.hotspot.{hotspot}.bed.gz",
        signal="data/roadmap/{signal}.pval.signal.bedGraph.gz",
        signal_tbi="data/roadmap/{signal}.pval.signal.bedGraph.gz.tbi",
    output:
        "result/{sid}.hotspot.{hotspot}.{signal}.tsv"
    log: "log/{sid}.hotspot.{hotspot}.{signal}.log"
    params:
        label=lambda wildcards: f"cragr.{wildcards.sid}.{wildcards.hotspot}.{wildcards.signal}",
        main_script=lambda wildcards: find_main_script(),
    shell:
        """
        Rscript {params.main_script} signal \
        -i {input.hotspot} \
        --signal {input.signal} \
        --output {output} \
        --verbose 2>&1 | tee {log}
        """
