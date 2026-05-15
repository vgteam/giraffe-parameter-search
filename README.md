# giraffe-parameter-search
Workflow for testing vg giraffe parameters


Experiment config file config.yaml  
Main snakefile 'Snakefile'  
Sub-snakefile for graphing 'make_graphs.smk'  
Graphing util 'scatter.py'  
Parameter search util 'parameter_search.py'  

## TODO
- Make graphing work (have snakefile correctly pass independent and dependent variable wildcards)
- Add capability for more graphs
- Add capability for more realness, tech, subset settings
- Add variant calling sub snakefile
- Add more documentation (dedicated wiki)