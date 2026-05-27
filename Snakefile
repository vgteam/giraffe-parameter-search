"""
Workflow file for parameter testing. 

Snakemake works backwards to create requested files. Here, the requested files are the .png for each graph 
described in the config and the .gam file for each alignment.

Workflow process sketch:
    - Rule all requests graphs specified in config. If there are no graphs listed, still requests final statistic tsv.
    - Rule giraffe_real_reads aligns reads according to config specification.
    - Rule stats_tsv creates a tsv for each parameter set of all statistics available from mapping. 
    - Rule stats_tsv_aggregate combines these into a single comprehensive tsv sorted by hash (includes parameter values).
    - Additional included snakefile "make_plots.smk" creates requested graphs using aggregate tsv.
"""

import parameter_search #file used for converting hashes to parameter strings
import pandas as pd 
import subprocess
import yaml
import datetime

configfile: "snake_config.yaml"

supported_statistics = ["mapped", "reads_correct", "reads_incorrect", "accuracy_percent", "mapq60", "wrong_mapq60", "softclips", "clipped_or_unmapped", "speed", 
        "speed_from_log", "memory(GB)"]


# =============================================
# Define macros for stats collection and graphing 
# =============================================

# Parameters
PARAM_SEARCH = parameter_search.ParameterSearch() #function used for converting hashes to parameter strings
experiment = config["experiment"]
exp_config = config["experiments"][experiment]

params_df = pd.read_csv(config["params_file"], sep="\t", comment=None, index_col=0) #stores all information stored in params_file
params_df.index = params_df.index.astype(str) # we reference hashes to fill in values, but hashes in {param_set} are strings
PARAM_SETS = [str(hash) for hash in list(params_df.index)] #list of all parameter set hashes


# Mapping tsv statistics
STATS_INC = config["included_statistics"] # stats to include in mapping_stats.tsv
if "all" in STATS_INC:
    STATS_INC = supported_statistics
STATS_EXC = config["excluded_statistics"] # stats to exclude in mapping_stats.tsv

# check that the stats for mapping are something we can collect
for stat in STATS_INC:
    if stat not in supported_statistics:
        raise ValueError(f"Unsupported statistic '{a}' in included_statistics")

for stat in STATS_EXC:
    if stat not in supported_statistics:
        raise ValueError(f"Unsupported statistic '{b}' in excluded_statistics")


# Graphing statistics
X_VARS = params_df.columns.tolist() # names of parameters tested
Y_VARS = config["desired_2d_graphs"] # y-variables for 2d plots
YZ_VARS = config["desired_3d_graphs"] # lists of y and z variables for 3d plots

# check that the stats in the requested graphs are something we can collect
for stat in Y_VARS + [stat for pair in YZ_VARS for stat in pair]:
    if stat not in supported_statistics:
        raise ValueError(f"Unsupported statistic '{stat}' in desired_graphs")
    if stat in STATS_EXC:
        raise ValueError(f"A statistic ('{stat}') cannot be in excluded_statistics but also in desired_graphs")

# =============================================
# Define supported wildcard/config values
# =============================================

wildcard_constraints:
    sample = "HG002",
    realness = "real|sim",
    tech = "hifi|r10y2025",
    subset = "1k|10k|100k|1m",
    preset = "hifi|r10",
    root = re.escape(config["root"])

# =============================================
# Slurm and memory management utilities
# =============================================

MAPPER_THREADS = 64 #default threads for mapping
SLURM_PARTITIONS = [ 
    ("short", 60),
    ("medium", 12 * 60),
    ("long", 7 * 24 * 60)
]

def choose_partition(minutes):
    """
    Get a Slurm partition that can fit a job running for the given number of
    minutes, or raise an error.
    """
    for name, limit in SLURM_PARTITIONS:
        if minutes <= limit:
            return name
    raise ValueError(f"No Slurm partition accepts jobs that run for {minutes} minutes")

def numerical_subset(subset):
    """
    Convert string subset values to numerical values.
    Used for calculating auto_mapping_threads().
    """
    match subset:
        case "1k":
            return 1000
        case "10k":
            return 10000
        case "100k":
            return 100000
        case "1m":
            return 1000000
    raise ValueError(f"{subset} is an unsupported subset. There may not be a fasta reference for this subset.")

def auto_mapping_threads(wildcards):
    """
    Choose the number of threads to use map reads, based on the subset.
    """
    mapping_threads = 0
    number = numerical_subset(exp_config["subset"])

    if number > 100000:
        mapping_threads = MAPPER_THREADS
    elif number > 10000:
        mapping_threads = 16
    else:
        mapping_threads = 8

    return mapping_threads

def auto_mapping_memory(wildcards):
    """
    Determine the memory to use for Giraffe mapping, in MB, based on tech.
    """
    thread_count = auto_mapping_threads(wildcards)

    base_mb = 60000

    if wildcards["tech"] == "illumina" or wildcards["tech"] == "element":
        scale_mb = 200000
    elif wildcards["tech"] == "hifi":
        scale_mb = 240000
    elif wildcards["tech"] == "r10":
        scale_mb = 600000
    else:
        scale_mb = 210000

    # Scale down memory with threads
    return scale_mb / 64 * thread_count + base_mb

# =============================================
# File utilities
# =============================================

def find_fastq(tech):
    """
    Determine .fq reference file for mapping based on realness, tech, sample, and subset.
    Maybe change this to match realness after making experiments for everything. 
    """
    reference = config.get("reference", None)
    if reference:
        return reference
    else:
        directory = config["reads_dir"]
        realness = exp_config["realness"]
        tech = exp_config["tech"]
        sample = exp_config["sample"]
        subset = exp_config["subset"]
        if tech == "r10y2025":
            return (f"{directory}/{realness}/{tech}/{sample}/{sample}-{realness}-{tech}-{subset}.fq")
        if tech == "hifi":
            return (f"{directory}/{realness}/{tech}/{sample}/HG002Revio_hg002v1.0.1_hifi_revio_pbmay24.pri.unshuffled.{subset}.fq")
        else:
            raise ValueError(f"Unsupported tech: {tech}")


def graph_names():
    """
    Makes a list of names of all graphs requested in the config file.
    This is sufficiently complicated to get its own function, especially once more types of plots are added
    """
    filenames = []
    for x in X_VARS:
        for y in Y_VARS:
            filenames.append(f"{config['root']}/results/{config['experiment']}/graphs/{exp_config['tech']}.{exp_config['sample']}.{exp_config['subset']}.{x}_vs_{y}.png")
        for y, z in YZ_VARS:
            filenames.append(f"{config['root']}/results/{config['experiment']}/graphs/{exp_config['tech']}.{exp_config['sample']}.{exp_config['subset']}.{x}_vs_{y}_vs_{z}.png")
    return filenames

include: "make_plots.smk"


# =============================================
# Entry point
# =============================================

rule all:
    """
    Final list of files requested. Snakemake looks here first when running the snakefile. 
    """
    input:
        expand("{root}/results/{experiment}/stats/manifest.yaml", root=config["root"], experiment=config["experiment"]),
        expand("{root}/results/{experiment}/stats/{tech}.{sample}.{subset}.mapping_stats.tsv", root=config["root"], experiment=config["experiment"], tech=exp_config["tech"], sample=exp_config["sample"], subset=exp_config["subset"]),
        graph_names()

rule manifest:
    """
    Generates run/experiment manifest with wildcards and run information.
    """
    output:
        yaml = "{root}/results/{experiment}/stats/manifest.yaml"
    run:
        manifest = {
            "timestamp": datetime.datetime.now().isoformat(timespec='seconds'),
            "vg_version": subprocess.getoutput("vg version"),
            "experiment": wildcards.experiment,
            "experiment_config": exp_config,
            "reference": find_fastq(exp_config["tech"]),
            "notes": "This file was generated by running an experiment in giraffe-parameter-search."
        }
        with open(output[0], "w") as f:
            yaml.dump(manifest, f, default_flow_style=False, sort_keys=False)

# =============================================
# Alignment
# =============================================

rule giraffe_real_reads:
    """
    Maps real reads with vg giraffe.
    """
    input:
        gbz = exp_config["gbz_graph"],
        dist = exp_config["distance_index"],
        zipfile = exp_config["zipfile"],
        minindex = exp_config["minimizer_index"],
        fastq = find_fastq(exp_config["tech"])
    output:
        "{root}/results/{experiment}/{param_set}/{tech}.{sample}.{subset}.real.{param_set}.gam"
    log: 
        "{root}/results/{experiment}/{param_set}/{tech}.{sample}.{subset}.real.{param_set}.log"
    wildcard_constraints:
        realness="real"
    benchmark: 
        "{root}/results/{experiment}/{param_set}/{tech}.{sample}.{subset}.real.{param_set}.benchmark"
    threads: auto_mapping_threads
    resources:
        mem_mb=auto_mapping_memory,
        runtime=1200,
        slurm_partition=choose_partition(1200),
    params:
        preset = exp_config["preset"],
        flags = lambda wildcards: expand(PARAM_SEARCH.hash_to_parameter_string(wildcards.param_set))
    run:
        # try block and if error throw exception print bottom of log block
        try:
            shell("vg giraffe -t{threads} --progress --parameter-preset {params.preset} -Z {input.gbz} -d {input.dist} -m {input.minindex} -f {input.fastq} -z {input.zipfile} {params.flags} --output-format gam >{output} 2>{log}")
        except Exception as e:
            # print bottom 10 lines of log block
            with open(log[0], 'r') as f:
                lines = f.readlines()
                tail = "".join(lines[-10:])
            raise Exception(f"vg giraffe command failed for {params.flags}:\n{tail}") from e


rule giraffe_sim_reads:
    """
    Maps simulated reads with vg giraffe. Used to calculate reads_correct.
    """
    input:
        gbz = exp_config["gbz_graph"],
        dist = exp_config["distance_index"],
        zipfile = exp_config["zipfile"],
        minindex = exp_config["minimizer_index"],
        gam= config["reads_dir"] + "/sim/{tech}/{sample}/{sample}-sim-{tech}-{subset}.gam"
    output:
        "{root}/results/{experiment}/{param_set}/{tech}.{sample}.{subset}.sim.{param_set}.gam"
    log:
        "{root}/results/{experiment}/{param_set}/{tech}.{sample}.{subset}.sim.{param_set}.log"
    wildcard_constraints:
        realness="sim"
    threads: auto_mapping_threads
    resources:
        mem_mb=auto_mapping_memory,
        runtime=600,
        slurm_partition=choose_partition(600)
    params:
        preset = exp_config["preset"],
        flags = lambda wildcards: expand(PARAM_SEARCH.hash_to_parameter_string(wildcards.param_set))
    run:
        try:
            shell("vg giraffe -t{threads} --progress --parameter-preset {params.preset} --track-provenance --set-refpos -Z {input.gbz} -d {input.dist} -m {input.minindex} -G {input.gam} -z {input.zipfile} {params.flags} --output-format gam >{output} 2>{log}")
        except Exception as e:
            # print bottom 10 lines of log block
            with open(log[0], 'r') as f:
                lines = f.readlines()
                tail = "".join(lines[-10:])
            raise Exception(f"vg giraffe command failed for {params.flags}:\n{tail}") from e

# =============================================
# Statistics
# =============================================

rule compare_alignments:
    """
    Run vg gamcompare on sim alignment gam and reference gam.
    This generates:
        - a txt file with number of correct reads used in stats_tsv
        - an annotated gam file used in vg_filter
    """
    input:
        gam = "{root}/results/{experiment}/{param_set}/{tech}.{sample}.{subset}.sim.{param_set}.gam",
        truth_gam = lambda wildcards: expand(
            "{reads_dir}/sim/{tech}/{sample}/{sample}-sim-{tech}-{subset}.gam", 
            reads_dir=config["reads_dir"], tech=wildcards.tech, sample=exp_config["sample"], subset=wildcards.subset
            )
    output:
        gam = "{root}/results/{experiment}/compared/{tech}.{sample}.{subset}.sim.{param_set}.compared.gam",
        tsv = "{root}/results/{experiment}/compared/{tech}.{sample}.{subset}.sim.{param_set}.compared.tsv",
        compare = "{root}/results/{experiment}/compared/{tech}.{sample}.{subset}.sim.{param_set}.compare.txt"
    params:
        # In the v2.0 CHM13-based graphs we now use a non-HG002 Y, and that's where our
        # Y truth positions are. For other graphs (like v2 prerelease 3 or v1.1) we won't
        # have the right CHM13 Y.
        # So we just ignore Y in CHM13. 
        ignore_flag="--ignore 'CHM13#0#chrY'"
    threads: 16
    resources:
        mem_mb=200000,
        runtime=800,
        slurm_partition=choose_partition(800)
    run:
        # Note that vg gamcompare computes the eligible read count directly
        # from the truth, so if reads vanish altogether (like with
        # GraphAligner's unmapped reads), they will count as wrong if they were
        # eligible.
        try:
            shell("vg gamcompare --threads 16 --range 200 {params.ignore_flag} {input.gam} {input.truth_gam} --output-gam {output.gam} -T > {output.tsv} 2>{output.compare}")
        except subprocess.CalledProcessError as error:
            with open(log[0]) as f:
                log_tail = f.readlines()[-10:] #grab last 10 lines of log file and print
            raise RuntimeError(f"command 'vg gamcompare' failed with exit code {error.returncode}.") from error


rule vg_filter:
    """
    Run vg filter on compared alignment gam files.
    This generates:
        - a tsv with name, score, correctness, softclip_total, identity, mapping_quality, and length used in stats_tsv
    """
    input: 
        gam = "{root}/results/{experiment}/compared/{tech}.{sample}.{subset}.sim.{param_set}.compared.gam"
    output:
        tsv = "{root}/results/{experiment}/{param_set}/{tech}.{sample}.{subset}.{param_set}.vg_filter_stats.tsv"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=100,
        slurm_partition=choose_partition(100)
    shell:
        "vg filter -t {threads} -T \"name;score;correctness;softclip_total;identity;mapping_quality;length\" {input.gam} > {output.tsv}" 

rule vg_stats:
    """
    Run vg stats on real alignment gam files.
    This generates:
        - a txt file with alignment score, mapping quality, and speed used in stats_tsv
    """
    input: 
        gam = "{root}/results/{experiment}/{param_set}/{tech}.{sample}.{subset}.real.{param_set}.gam"
    output:
        txt = "{root}/results/{experiment}/{param_set}/{tech}.{sample}.{subset}.real.{param_set}.gamstats.txt"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=100,
        slurm_partition=choose_partition(100)
    shell:
        "vg stats -p {threads} -a {input.gam} >{output.txt}"


rule stats_tsv:
    """
    Create tsvs containing supported_statistics for each parameter set, including the paramters as columns at the front.
    The individual tsvs are aggregated into the complete tsv in stats_tsv_aggregate. This is because Snakemake
    requires the output file name to include all wildcards the input files contain. 
    """
    input:
        mapping_log = "{root}/results/{experiment}/{param_set}/{tech}.{sample}.{subset}.real.{param_set}.log",
        vg_stats = "{root}/results/{experiment}/{param_set}/{tech}.{sample}.{subset}.real.{param_set}.gamstats.txt",
        gamcompare = "{root}/results/{experiment}/compared/{tech}.{sample}.{subset}.sim.{param_set}.compare.txt",
        #clipped_or_unmapped = "{root}/results/{experiment}/{param_set}/{tech}.{sample}.{subset}.{param_set}.clipped_or_unmapped.tsv",
        vg_filter = "{root}/results/{experiment}/{param_set}/{tech}.{sample}.{subset}.{param_set}.vg_filter_stats.tsv"
    output:
        temp("{root}/results/{experiment}/stats/{tech}.{sample}.{subset}.{param_set}.mapping_stats.tsv")
    threads: 1
    resources:
        mem_mb=2000,
        runtime=100,
        slurm_partition=choose_partition(100)
    run:
        # copy from df with all paramter hashes and values, but only the row with our current parameter set
        df = params_df.loc[[wildcards.param_set]].copy()
        # add headers for included stats
        df[STATS_INC] = -1 
        # remove columns of excluded stats
        df = df.drop(columns=STATS_EXC)

        # get speed_from_log, memory from log file
        with open(input.mapping_log, "r") as f:
            text = f.read()
        if "speed_from_log" in STATS_INC:
            speed = float(re.search(r"Mapping speed:\s+([\d.]+)\s+reads per second per thread", text).group(1))
            df["speed_from_log"] = round(speed, 2)
        if "memory(GB)" in STATS_INC:
            memory= float(re.search(r"Memory footprint:\s+([\d.]+)\s+GB", text).group(1))
            df["memory(GB)"] = round(memory, 2)

        # get speed(r/p/s) from vg stats file
        with open(input.vg_stats, "r") as f:
            text = f.read()
        if "speed" in STATS_INC:
            speed = float(re.search(r"Speed:\s+([\d.]+)\s+reads/second", text).group(1))
            df["speed"] = round(speed, 2)

        # get reads_correct and accuracy_percent from compare file
        with open(input.gamcompare, "r") as f:
            text = f.read()
        if "reads_correct" in STATS_INC:
            reads_correct = int(re.search(r"([\d,]+)\s+reads correct", text).group(1).replace(",", ""))
            df["reads_correct"] = reads_correct
        if "accuracy_percent" in STATS_INC:
            accuracy_percent = float(re.search(r"([\d.]+)%\s+accuracy", text).group(1))
            df["accuracy_percent"] = accuracy_percent
        
        #get softclips, clipped_or_unmapped, mapped, mapq60, and wrong_mapq60 from vg filter file
        filter_df = pd.read_csv(input.vg_filter, sep='\t')

        if "reads_incorrect" in STATS_INC:
            reads_incorrect = (filter_df["correctness"] == "incorrect").sum()
            df["reads_incorrect"] = reads_incorrect

        # get from softclip_total
        if "softclips" in STATS_INC:
            softclips = filter_df["softclip_total"].sum()
            df["softclips"] = softclips

        # look for read names where read is unmapped (score = 0), add those lengths
        # add softclips
        if "clipped_or_unmapped" in STATS_INC:
            clipped_or_unmapped = (filter_df[filter_df["score"] == 0]["length"].sum()) + softclips
            df["clipped_or_unmapped"] = clipped_or_unmapped

        # unmapped reads have score=0
        if "mapped" in STATS_INC:
            mapped = (filter_df["score"] != 0).sum()
            df["mapped"] = mapped

        # count number of reads with mq=60
        if "mapq60" in STATS_INC:
            mapq60 = (filter_df["mapping_quality"] == 60).sum()
            df["mapq60"] = mapq60

        # count number of reads with mq=60 which have correctness = incorrect
        if "wrong_mapq60" in STATS_INC:    
            wrong_mapq60 = (filter_df[filter_df["mapping_quality"] == 60]["correctness"] == "incorrect").sum()
            df["wrong_mapq60"] = wrong_mapq60

        df.to_csv(output[0], sep="\t", index=True)

rule stats_tsv_aggregate:
    """
    Combine all per-param_set stats tsvs into one file.
    """
    input: 
        lambda wildcards: expand(
            "{root}/results/{experiment}/stats/{tech}.{sample}.{subset}.{param_set}.mapping_stats.tsv", 
            root=wildcards.root, experiment=wildcards.experiment, tech=wildcards.tech, sample=wildcards.sample, subset=wildcards.subset, param_set=PARAM_SETS
            )        
    output:
        "{root}/results/{experiment}/stats/{tech}.{sample}.{subset}.mapping_stats.tsv"
    run:
        dfs = [pd.read_csv(f, sep="\t", index_col=0) for f in input]
        pd.concat(dfs).to_csv(output[0], sep="\t", index=True)


"""
rule save_results:
    
    Saves graphs, mapping stats tsv, parameter search config
    
    input:
        path = config["results_path"]
"""