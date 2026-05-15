# giraffe-parameter-search
Workflow for testing vg giraffe parameters

Experiment config file config.yaml  
Main snakefile 'Snakefile'  
Sub-snakefile for graphing 'make_graphs.smk'  
Graphing util 'scatter.py'  
Parameter search util 'parameter_search.py'
Parameter search util config file parameter_search_config.tsv

## How To Use
1. Edit `parameter_search_config.tsv` to contain the parameters you want to generate values for, and what way you want to generate them. See `parameter_search.py` for more information.
2. Run:
    ```
    python3 parameter_search.py --count x
    ```
    This will create the `hash_to_parameters` file.
3. Specify the graphs you want to generate in `config.yaml`.  
3. Run:

    ```
    (umask 002; snakemake --configfile config.yaml --config experiment=map_hifi_10k -j128 --rerun-incomplete --use-singularity --singularity-args "-B /private" --latency-wait 120 --executor slurm --keep-going)
    ```
    Replace map_hifi_10k with experiment of choice from config.

## TODO
- Make graphing work (have snakefile correctly pass independent and dependent variable wildcards)
- Add capability for more graphs
- Add capability for more realness, tech, subset settings
- Add variant calling sub snakefile
- Add more documentation (dedicated wiki)