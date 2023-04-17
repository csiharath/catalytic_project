# Catalytic project

`model_correction.py` is a script adding informations into an sbml model such as KEGG ID and genes names, and preprocessing data of the model to create lists talen as argument for km predictions

## Requirements

- 

## Usage 

`
python3 code/model_correction.py -m <sbmlmodel> -i -g -d <dataset> -k -c
` 
### Required argument

* `-m, --model <sbmlmodel>` : Metabolic model in sbml format.

### Optional arguments

* `-i, --keggid` : Option to add kegg id in the model. Requires to add a dataset with option `--data`.
* `-g, --genename` : Option to add gene's name in the model. Requires to add a dataset with option `--data`.
* `-d, --data <dataset>` : Csv file containing colums of:
  - reaction's name in the model
  - the corresponding KEGG ID
  - a list og associated genes

*  `-k, --km` : Option to find parametes to calculates km.
*  `-c, --kcat` : Option to find parametes to calculates kcat.

