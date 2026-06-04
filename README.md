# giraffe-parameter-search
Workflow for testing vg giraffe parameters

Contents:
* Experiment config file `config.yaml` 
* Main snakefile `Snakefile`  
* Sub-snakefile for graphing `make_plots.smk`  
* Graphing util `scatter.py`  
* Parameter search util `parameter_search.py`
* Parameter search util config file `parameter_search_config.tsv`


> [!NOTE]
> **for more information, visit the Wiki tab.** 

## How To Use
1. Edit `parameter_search_config.tsv` to contain the parameters you want to generate values for, and what distribution you want to sample from. See `parameter_search.py` and the wiki for more information.
2. Run:
    ```
    python3 parameter_search.py --count x --benchmark-default
    ```
    This will create the `hash_to_parameters` file. You do NOT need to delete the old version of this file before running the script again. It will delete the file for you unless using `--append`.
3. Specify the statistics you want to measure and graphs you want to generate in `config.yaml`.  
3. Run:

    ```
    snakemake --configfile config.yaml --config experiment=my_experiment --rerun-incomplete --executor slurm --keep-going
    ```
    Replace `my_experiment` with experiment of choice from config.
4. Your results will include a table of statistics, plots, and a manifest including the datetime, vg version, and config parameters.

## Dependencies 
Have the following:
* vg: [Install vg here](https://github.com/vgteam/vg/tree/master#installation)
* pandas
* yaml
* scipy
* bidict
* matplotlib

You can run:
```
conda create -n giraffe-param-search -c conda-forge -c bioconda snakemake=9.13.7 \
    snakemake-executor-plugin-slurm \
    bidict matplotlib pandas pyyaml scipy
```

## TODO
- Add mismapped statistic
- Add variant calling sub snakefile
