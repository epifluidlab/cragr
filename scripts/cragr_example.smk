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



WALL_TIME_MAX = 2880
MEM_PER_CORE = 1800


rule ifs:
    input:
        frag="{prefix}/frag/{sid}.frag.gz",
        frag_idx="{prefix}/frag/{sid}.frag.gz.tbi",
        mappability="{prefix}/{DATA_DIR}/mappability.hs37-1kg.w200.s20.0_9.bed.gz",
        mappability_idx="{prefix}/{DATA_DIR}/mappability.hs37-1kg.w200.s20.0_9.bed.gz.tbi",
    output: 
        ifs=temp("{prefix}/temp/{sid}.ifs.chr{chrom}.bedGraph.gz"),
        ifs_idx=temp("{prefix}/temp/{sid}.ifs.chr{chrom}.bedGraph.gz.tbi")
    log: "{prefix}/temp/{sid}.ifs.chr{chrom}.log"
    params:
        label=lambda wildcards: f"ifs.{wildcards.sid}.chr{wildcards.chrom}",
        script_home=lambda wildcards: SCRIPT_HOME
    threads: lambda wildcards, attempt: int(2 * (0.5 + 0.5 * attempt))
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
        --output-ifs "$tmpdir"/output.bedGraph.gz \
        -m {input.mappability} \
        --chrom {wildcards.chrom} \
        --verbose

        mv "$tmpdir"/output.bedGraph.gz {output.ifs}.tmp
        mv {output.ifs}.tmp {output.ifs}
        mv "$tmpdir"/output.bedGraph.gz.tbi {output.ifs_idx}.tmp
        mv {output.ifs_idx}.tmp {output.ifs_idx}
        """

