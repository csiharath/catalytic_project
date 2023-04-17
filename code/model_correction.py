"""
Preprocessing data for km prediction
Usage:
    model_correction.py -m <sbmlmodel> -k -g -d <dataset>

Positional arguments:
    sbmlmodel            Metabolic model in sbml format.
    dataset              Csv file containing colums of reaction's name in the model, the corresponding KEGG ID, and a list og associated genes
"""

import argparse
import textwrap
import time

import cobra
import pandas as pd
import numpy as np
import pickle
import model_correction_utils as mcu


############################# Command Line Interface #############################

parser = argparse.ArgumentParser(description='Preprocessing of dataset for km and kcat prediction',\
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-m', '--model', dest='sbmlmodel', action='store',\
    required=True,\
    help=textwrap.dedent('''\
        Metabolic model in sbml format.
        \n'''))
parser.add_argument('-k', '--keggid', dest='kegg_id', action='store_true',\
    help=textwrap.dedent('''\
        Option to add kegg id in the model. Requires to add a dataset with option --data.
        \n'''))
parser.add_argument('-g', '--genename', dest='gene_name', action='store_true',\
    help=textwrap.dedent('''\
        Option to add gene's name in the model. Requires to add a dataset with option --data.
        \n'''))
parser.add_argument('-d', '--data', dest='dataset', action='store', default='',\
    help=textwrap.dedent('''\
        Csv file containing colums of:
            - reaction's name in the model
            - the corresponding KEGG ID
            - a list og associated genes
        \n'''))
    
############################# Raising parsing errors #############################

args = parser.parse_args()

if (args.kegg_id or args.gene_name) and args.dataset == '':
    parser.error("--keggid and --genename require to add dataset with --data")

############################# Variables Declaration ##############################

start_time = time.time()

if args.kegg_id or args.gene_name:
    dataset = args.dataset

model_file = args.sbmlmodel
add_keggid = args.kegg_id
add_genes = args.gene_name

model, errors = cobra.io.validate_sbml_model(model_file)

####################### Adding informations into the model #######################

if add_keggid or add_genes:
    df_react = pd.read_csv(dataset, sep=";")
    df_react = df_react.set_index('Reaction name')

if add_keggid:
    print(f"Adding kegg id to {model_file} model.")
    mcu.add_kegg_id_to_model(model, df_react)

if add_genes:
    print(f"Adding genes names to {model_file} model.")
    mcu.add_gene_reaction_rule(model, df_react)

if add_keggid or add_genes:
    cobra.io.write_sbml_model(model, filename=model_file.split(".xml")[0]+"_updated.xml")


######################## Defining km prediction parameters #######################

enzyme_list = [(enzyme.id, enzyme.notes['kegg_id']) for enzyme in model.reactions if 'kegg_id' in enzyme.notes]

print("Looking for km predection parameters in database.")
dict_km_parameters = mcu.find_compounds_AAseq(model, enzyme_list)
dict_aaseq = dict_km_parameters[1]
dict_km_parameters = dict_km_parameters[0]

km_arguments = mcu.create_km_arguments(dict_km_parameters)

substrates = km_arguments[0]
enzymes = km_arguments[1]

##################################################################################

unique_compounds = np.unique(substrates).tolist()

dict_compounds = mcu.build_dict_compounds(unique_compounds)

##################################################################################

print(len(enzymes))
print(len(substrates))
print(len(dict_aaseq.keys()))
print(len(dict_compounds.keys()))

# with open("data/enzyme.p", "wb") as e:
#     pickle.dump(enzymes, e)

# with open("data/substrat.p", "wb") as s:
#     pickle.dump(substrates, s)

# with open("data/compound.p", "wb") as c:
#     pickle.dump(dict_compounds, c)

# with open("data/aaseq.p", "wb") as aa:
#     pickle.dump(dict_aaseq, aa)

# Windows :
with open("data\enzyme.p", "wb") as e:
    pickle.dump(enzymes, e)

with open("data\substrat.p", "wb") as s:
    pickle.dump(substrates, s)

with open("data\compound.p", "wb") as c:
    pickle.dump(dict_compounds, c)

with open("data\seq.p", "wb") as aa:
    pickle.dump(dict_aaseq, aa)

print("--- %s seconds ---" % (time.time() - start_time))