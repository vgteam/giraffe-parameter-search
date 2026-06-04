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
1. Edit `parameter_search_config.tsv` to contain the parameters you want to generate values for, and what way you want to generate them. See `parameter_search.py` for more information.
2. Run:
    ```
    python3 parameter_search.py --count x
    ```
    This will create the `hash_to_parameters` file.
3. Specify the statistics you want to measure and graphs you want to generate in `config.yaml`.  
3. Run:

    ```
    snakemake --configfile config.yaml --config experiment=my_experiment --rerun-incomplete --executor slurm --keep-going
    ```
    Replace my_experiment with experiment of choice from config.

## Dependencies 
Have the following:
* vg
* pandas
* yaml
* scipy
* bidict
* matplotlib

## TODO
- Add mismapped statistic
- Add variant calling sub snakefile
