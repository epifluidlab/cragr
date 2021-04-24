# cragr pipeline

localrules: all, merge_ifs, call_hotspot_pois, call_hotspot_nb

# >>> Configuration >>>
FULL_CORES = config.get("FULL_CORES", 16)
PART_CORES = config.get("PART_CORES", 4)
MEM_PER_CORE = config.get("MEM_PER_CORE", 1800)
WALL_TIME_MAX = config.get("WALL_TIME_MAX", 2880)
CORE_FACTOR = float(config.get("CORE_FACTOR", 1))

SCRIPT_HOME = config.get("SCRIPT_HOME")
WALL_TIME_MAX = 2880
MEM_PER_CORE = 1800


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


rule stage_ifs_chrom:
    input:
        frag="frag/{sid}.frag.chr{chrom}.bed.gz",
        frag_idx="frag/{sid}.frag.chr{chrom}.bed.gz.tbi",
        # frag="frag/{sid}.frag.bed.gz",
        # frag_idx="frag/{sid}.frag.bed.gz.tbi",
        mappability="data/mappability.hs37-1kg.w200.s20.0_9.bed.gz",
        mappability_idx="data/mappability.hs37-1kg.w200.s20.0_9.bed.gz.tbi",
    output: 
        ifs=temp("temp/{sid}.ifs.raw.chr{chrom}.bed.gz")
    log: "log/{sid}.chr{chrom}.stage_ifs.log"
    params:
        label=lambda wildcards: f"cragr.stage_ifs.{wildcards.sid}.chr{wildcards.chrom}",
        script_home=lambda wildcards: SCRIPT_HOME
    threads: lambda wildcards, input, attempt: int(stage_ifs_cores(wildcards.chrom, input.frag, single_chrom=True) * (0.5 + 0.5 * attempt))
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
        -o "$tmpdir"/output.bed.gz \
        -m {input.mappability} \
        --chrom {wildcards.chrom} \
        --verbose 2>&1 | tee {log}

        mv "$tmpdir"/output.bed.gz {output.ifs}.tmp
        mv {output.ifs}.tmp {output.ifs}
        """


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
        script_home=lambda wildcards: SCRIPT_HOME
    threads: lambda wildcards, input, attempt: int(stage_ifs_cores(wildcards.chrom, input.frag) * (0.5 + 0.5 * attempt))
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
        -o "$tmpdir"/output.bed.gz \
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
        script_home=lambda wildcards: SCRIPT_HOME
    threads: lambda wildcards, input, attempt: int(stage_peak_cores(wildcards.chrom) * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=300,
        attempt=lambda wildcards, threads, attempt: attempt
    shell:
        """
        tmpdir=$(mktemp -d)

        Rscript {params.script_home}/cragr.R peak \
        -i {input.ifs} \
        -o "$tmpdir"/ifs.bedGraph.gz \
        -m {input.mappability} \
        --chrom {wildcards.chrom} \
        --verbose 2>&1 | tee {log}

        mv "$tmpdir"/ifs.bedGraph.gz {output.ifs}.tmp
        mv {output.ifs}.tmp {output.ifs}
        """


rule merge_ifs:
    input: 
        ifs=expand("temp/{{sid}}.ifs.chr{chrom}.bedGraph.gz", chrom=range(21, 23)),
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


rule call_hotspot_nb:
    input: 
        ifs="result/{sid}.ifs.bedGraph.gz",
        ifs_idx="result/{sid}.ifs.bedGraph.gz.tbi",
    output:
        hotspot="result/{sid}.hotspot.nb.bed.gz",
    log: "log/{sid}.hotspot.nb.bed.gz",
    params:
        label=lambda wildcards: f"cragr.hotspot_nb.{wildcards.sid}",
        script_home=lambda wildcards: SCRIPT_HOME
    shell:
        """
        tmpdir=$(mktemp -d)

        Rscript {params.script_home}/cragr.R hotspot \
        -i {input.ifs} \
        -o "$tmpdir"/hotspot.bed.gz \
        --fdr 0.25 \
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
    log: "log/{sid}.hotspot.pois.bed.gz",
    params:
        label=lambda wildcards: f"cragr.hotspot_pois.{wildcards.sid}",
        script_home=lambda wildcards: SCRIPT_HOME
    shell:
        """
        tmpdir=$(mktemp -d)

        Rscript {params.script_home}/cragr.R hotspot \
        -i {input.ifs} \
        -o "$tmpdir"/hotspot.bed.gz \
        --fdr 0.01 \
        --pval 0.00001 \
        --hotspot-method pois \
        --verbose 2>&1 | tee {log}

        mv "$tmpdir"/hotspot.bed.gz {output.hotspot}.tmp
        mv {output.hotspot}.tmp {output.hotspot}
        """    
