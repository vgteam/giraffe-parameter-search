"""
This file is used to create the graphs requested in config.yaml.
Requires running the main Snakefile to give it information on what to plot. You should not need to run this snakefile individually.

Uses file scatter.py for plotting, which should be in the same directory as this file.
"""

rule scatter_2d:
    """
    Create a scatter plot with parameter values and ONE dependent stat variable.
    """
    input:
        tsv = "{root}/results/{experiment}/stats/{tech}.{sample}.{subset}.mapping_stats.tsv"
    output:
        png = "{root}/results/{experiment}/graphs/{tech}.{sample}.{subset}.{x}_vs_{y}.png"
    threads: 1
    resources:
        mem_mb=512,
        runtime=10,
        slurm_partition=choose_partition(10)
    run:
        infile = open(input.tsv)
        header = infile.readline().split()
        parameter_col = str(header.index(wildcards.x)+1)
        stat_col = str(header.index(wildcards.y)+1) 
        infile.close()
        shell("cat {input.tsv} | grep -v '#' | awk '{{print $" + parameter_col + " \"\\t\" $" + stat_col + "}}' | ./scatter.py --title '" + wildcards.x + " vs. " + wildcards.y + "' --x_label " + wildcards.x + " --y_label '" + wildcards.y + "' --legend_overlay 'best' --save {output.png} /dev/stdin")

rule scatter_3d:
    """
    Create a scatter plot with parameter values and TWO dependent stat variables.
    """
    input:
        tsv = "{root}/results/{experiment}/stats/{tech}.{sample}.{subset}.mapping_stats.tsv"
    output:
        png = "{root}/results/{experiment}/graphs/{tech}.{sample}.{subset}.{x}_vs_{y}_vs_{z}.png"
    threads: 1
    resources:
        mem_mb=512,
        runtime=10,
        slurm_partition=choose_partition(10)
    run:
        infile = open(input.tsv)
        header = infile.readline().split()
        parameter_col = str(header.index(wildcards.x)+1)
        y_stat_col = str(header.index(wildcards.y)+1) 
        z_stat_col = str(header.index(wildcards.z)+1) 
        infile.close()
        #shell("cat <(cat {input.tsv} | grep -v '#' | awk '{{print $" + parameter_col + " \"\\t\" $" + y_stat_col + "}}' | sed 'Y PARAMETER') <(cat {input.tsv} | grep -v '#' | awk '{{print $" + parameter_col + " \"\\t\" $" + z_stat_col + "}}' | sed 'Z PARAMETER') | ./scatter.py --title '" + wildcards.x + " vs. " + wildcards.y + " vs. " + wildcards.z "' --x_label " + wildcards.x + " --y_label '" + wildcards.y + "' '" + wildcards.z + "'  --y_per_category --categories 'Y PARAMETER NAME' 'Z PARAMETER NAME' --legend_overlay 'best' --save {output.plot} /dev/stdin")
        shell("cat <(cat {input.tsv} | grep -v '#' | awk '{{print $" + parameter_col + " \"\\t\" $" + y_stat_col + "}}' | sed 's/^/" + wildcards.y + "\\t/') <(cat {input.tsv} | grep -v '#' | awk '{{print $" + parameter_col + " \"\\t\" $" + z_stat_col + "}}' | sed 's/^/" + wildcards.z + "\\t/') | ./scatter.py --title '" + wildcards.x + " vs. " + wildcards.y + " vs. " + wildcards.z + "' --x_label " + wildcards.x + " --y_label '" + wildcards.y + "/" + wildcards.z + "' --y_per_category --categories '" + wildcards.y + "' '" + wildcards.z + "' --legend_overlay 'best' --save {output.png} /dev/stdin")
    
ruleorder: scatter_2d > scatter_3d