# cragr pipeline

# localrules: merge_ifs, signal_level

# >>> Configuration >>>
FULL_CORES = config.get("FULL_CORES", 16)
PART_CORES = config.get("PART_CORES", 4)
MEM_PER_CORE = config.get("MEM_PER_CORE", 1800)
WALL_TIME_MAX = config.get("WALL_TIME_MAX", 2880)
CPU_FACTOR = float(config.get("CPU_FACTOR", 1))

WALL_TIME_MAX = 2880
MEM_PER_CORE = 1800

FDR_NB = config.get("FDR_NB", 0.25)
FDR_POIS = config.get("FDR_POIS", 0.01)
GC_CORRECT = config.get("GC_CORRECT", "TRUE")

SCRIPT_PATH = config.get("SCRIPT_PATH", ".") 
DISABLE_PARALLEL = bool(config.get("DISABLE_PARALLEL", 0))
R = config.get("R", "Rscript")
# <<< Configuration <<<


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
        slurm_job_label=lambda wildcards: f"cragr.stage_ifs.{wildcards.sid}.chr{wildcards.chrom}",
    threads: lambda wildcards, input, attempt: int(CPU_FACTOR * stage_ifs_cores(wildcards.chrom, input.frag) * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=300,
        attempt=lambda wildcards, threads, attempt: attempt
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        {R} -e 'sessionInfo()'

        {R} {SCRIPT_PATH}/cragr.R ifs \
        -i {input.frag} \
        -o "$tmpdir"/output.bed.gz \
        --gc-correct \
        --genome GRCh37 \
        --exclude-region encode.blacklist.hs37-1kg \
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
        ifs=temp("temp/{sid}.ifs.{gc_type}.chr{chrom}.bedGraph.gz"),
    log: "log/{sid}.{gc_type}.chr{chrom}.stage_peak.log"
    params:
        slurm_job_label=lambda wildcards: f"cragr.stage_peak.{wildcards.sid}.{wildcards.gc_type}.chr{wildcards.chrom}",
    threads: lambda wildcards, input, attempt: int(CPU_FACTOR * stage_peak_cores(wildcards.chrom) * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=300,
        attempt=lambda wildcards, threads, attempt: attempt
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        if [ "{wildcards.gc_type}" == "gc" ]
        then
            gc_correct=TRUE
        else
            gc_correct=FALSE
        fi

        {R} -e 'sessionInfo()'

        {R} {SCRIPT_PATH}/cragr.R peak \
        -i {input.ifs} \
        -o "$tmpdir"/ifs.bedGraph.gz \
        --gc-correct \
        --genome GRCh37 \
        --exclude-region encode.blacklist.hs37-1kg \
        -m {input.mappability} \
        --chrom {wildcards.chrom} \
        --verbose 2>&1 | tee {log}

        mv "$tmpdir"/ifs.bedGraph.gz {output.ifs}.tmp
        mv {output.ifs}.tmp {output.ifs}
        """


rule merge_ifs:
    input:
        ifs=expand("temp/{{sid}}.ifs.{{gc_type}}.chr{chrom}.bedGraph.gz", chrom=range(1, 23)),
    output:
        ifs="result/{sid}.ifs.{gc_type,(gc|nogc)}.bedGraph.gz",
        ifs_idx="result/{sid}.ifs.{gc_type}.bedGraph.gz.tbi",
    threads: lambda wildcards, input, attempt: int(2 * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=300,
        attempt=lambda wildcards, threads, attempt: attempt
    params:
        slurm_job_label=lambda wildcards: f"cragr.merge_ifs.{wildcards.sid}.{wildcards.gc_type}",
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

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
        ifs="result/{sid}.ifs.{gc_type}.bedGraph.gz",
        ifs_idx="result/{sid}.ifs.{gc_type}.bedGraph.gz.tbi",
        chrom_sizes="data/human_g1k_v37.chrom.sizes",
    output: "result/{sid}.ifs.{gc_type,(gc|nogc)}.bw"
    threads: lambda wildcards, input, attempt: int(2 * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=300,
        attempt=lambda wildcards, threads, attempt: attempt
    params:
        slurm_job_label=lambda wildcards: f"cragr.ifs_bigwig.{wildcards.sid}.{wildcards.gc_type}",
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

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


# Dried IFS files: only keep regions whether the adjusted p-values are less than 0.75
rule ifs_dried:
    input: "result/{sid}.ifs.{gc_type}.bedGraph.gz",
    output:
        ifs="result/{sid}.ifs.{gc_type}.dried.bedGraph.gz",
        ifs_tbi="result/{sid}.ifs.{gc_type}.dried.bedGraph.gz.tbi",
    threads: lambda wildcards, input, attempt: int(2 * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=300,
        attempt=lambda wildcards, threads, attempt: attempt
    params:
        slurm_job_label=lambda wildcards: f"cragr.ifs_dried.{wildcards.sid}.{wildcards.gc_type}",
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        zcat {input} | bioawk -t 'substr($0,1,1)=="#" || ($16!="." && $16<=0.75) || ($10!="." && $10<=0.75)' | bgzip > $tmpdir/output.bedGraph.gz
        tabix -p bed $tmpdir/output.bedGraph.gz

        mv $tmpdir/output.bedGraph.gz {output.ifs}.tmp
        mv $tmpdir/output.bedGraph.gz.tbi {output.ifs_tbi}.tmp
        mv {output.ifs}.tmp {output.ifs}
        mv {output.ifs_tbi}.tmp {output.ifs_tbi}

        rm -rf $tmpdir
        """



rule call_hotspot:
    input:
        ifs="result/{sid}.ifs.{gc_type}.bedGraph.gz",
        ifs_idx="result/{sid}.ifs.{gc_type}.bedGraph.gz.tbi",
        chrom_sizes="data/human_g1k_v37.chrom.sizes",
    output:
        hotspot="result/{sid}.hotspot.{gc_type,(gc|nogc)}.bed.gz",
    threads: lambda wildcards, input, attempt: int(2 * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=300,
        attempt=lambda wildcards, threads, attempt: attempt
    params:
        slurm_job_label=lambda wildcards: f"cragr.hotspot.{wildcards.sid}",
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        set +e
        zcat {input.ifs} | \
        head -n 1000 | bioawk -t 'substr($1,1,1)=="#" && $1!="#chrom"' > $tmpdir/output.bed
        set -e

        echo "#chrom\tstart\tend" >> $tmpdir/output.bed
        bgzip $tmpdir/output.bed

        zcat {input.ifs} |
        bioawk -t 'substr($1,1,1)!="#" && $16<=0.2 && $16!="."' |
        bedtools slop -g {input.chrom_sizes} -i - -b 90 -header |
        bedtools merge -header -i - -d 200 |
        bgzip >> $tmpdir/output.bed.gz

        mv $tmpdir/output.bed.gz {output.hotspot}.tmp
        mv {output.hotspot}.tmp {output.hotspot}

        rm -rf $tmpdir
        """


# rule signal_level:
#     input:
#         hotspot="result/{sid}.hotspot.{hotspot}.bed.gz",
#         signal="data/roadmap/{signal}.pval.signal.bedGraph.gz",
#         signal_tbi="data/roadmap/{signal}.pval.signal.bedGraph.gz.tbi",
#     output:
#         "result/{sid}.hotspot.{hotspot}.{signal}.tsv"
#     log: "log/{sid}.hotspot.{hotspot}.{signal}.log"
#     params:
#         slurm_job_label=lambda wildcards: f"cragr.{wildcards.sid}.{wildcards.hotspot}.{wildcards.signal}",
#     shell:
#         """
#         Rscript {params.main_script} signal \
#         -i {input.hotspot} \
#         --signal {input.signal} \
#         --output {output} \
#         --verbose 2>&1 | tee {log}
#         """
