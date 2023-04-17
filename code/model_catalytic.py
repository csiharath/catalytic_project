"""
Preprocessing data for km prediction
Usage:
    model_correction.py -m <sbmlmodel> -i -g -d <dataset> -k -c

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
import model_catalytic_utils as mcu


############################# Command Line Interface #############################

parser = argparse.ArgumentParser(description='Preprocessing of dataset for km and kcat prediction',\
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-m', '--model', dest='sbmlmodel', action='store',\
    required=True,\
    help=textwrap.dedent('''\
        Metabolic model in sbml format.
        \n'''))
parser.add_argument('-i', '--keggid', dest='kegg_id', action='store_true',\
    help=textwrap.dedent('''\
        Option to add kegg id in the model. Requires to add a dataset with option --data.
        \n'''))
parser.add_argument('-g', '--genename', dest='gene_name', action='store_true',\
    help=textwrap.dedent('''\
        Option to add gene's name in the model. Requires to add a dataset with option --data.
        \n'''))
parser.add_argument('-k', '--km', dest='km', action='store_true',\
    help=textwrap.dedent('''\
        Option to find parametes to calculates km.
        \n'''))
parser.add_argument('-c', '--kcat', dest='kcat', action='store_true',\
    help=textwrap.dedent('''\
        Option to find parametes to calculates kcat.
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
    parser.error("--keggid and --genename require to add dataset with --data.")

if not args.km and not args.kcat:
    parser.error("At least one of the options --km or -kcat is required.")

############################# Variables Declaration ##############################

start_time = time.time()

if args.kegg_id or args.gene_name:
    dataset = args.dataset

model_file = args.sbmlmodel
add_keggid = args.kegg_id
add_genes = args.gene_name
km = args.km
kcat = args.kcat

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

if kcat and km:
    print("Looking for km and kcat prediction parameters in database.")
elif km:
    print("Looking for km prediction parameters in database.")
elif kcat:
    print("Looking for kcat prediction parameters in database.")

response_tuple = mcu.find_compounds_AAseq(model, enzyme_list)

dict_aaseq = response_tuple[1]
dict_km_parameters = response_tuple[0]
dict_genes_id = response_tuple[2]
arguments = mcu.create_km_kcat_arguments(dict_km_parameters, km, kcat)

km_sub = arguments[0]
km_enz = arguments[1]
kcat_prod = arguments[2]
kcat_sub = arguments[3]
kcat_enz = arguments[4]


##################################################################################


unique_compounds = np.unique([c for enzyme in dict_km_parameters for c in dict_km_parameters[enzyme]['substrates'] + dict_km_parameters[enzyme]['products']]).tolist()

dict_compounds = mcu.build_dict_compounds(unique_compounds)

##################################################################################

if km:
    with open("data/km_enzyme.p", "wb") as kme:
        pickle.dump(km_enz, kme)

    with open("data/km_substrat.p", "wb") as kms:
        pickle.dump(km_sub, kms)

if kcat:
    with open("data/kcat_enzyme.p", "wb") as kcate:
        pickle.dump(kcat_enz, kcate)

    with open("data/kcat_substrat.p", "wb") as kcats:
        pickle.dump(kcat_sub, kcats)

    with open("data/kcat_product.p", "wb") as kcatp:
        pickle.dump(kcat_prod, kcatp)

with open("data/compound.p", "wb") as c:
    pickle.dump(dict_compounds, c)

with open("data/seq.p", "wb") as aa:
    pickle.dump(dict_aaseq, aa)
# Windows :
# if km:
#     with open("data\km_enzyme.p", "wb") as kme:
#         pickle.dump(km_enz, kme)

#     with open("data\km_substrat.p", "wb") as kms:
#         pickle.dump(km_sub, kms)

# if kcat:
#     with open("data\kcat_enzyme.p", "wb") as kcate:
#         pickle.dump(kcat_enz, kcate)

#     with open("data\kcat_substrat.p", "wb") as kcats:
#         pickle.dump(kcat_sub, kcats)

#     with open("data\kcat_product.p", "wb") as kcatp:
#         pickle.dump(kcat_prod, kcatp)


# with open("data\compound.p", "wb") as c:
#     pickle.dump(dict_compounds, c)

# with open("data\seq.p", "wb") as aa:
#     pickle.dump(dict_aaseq, aa)

print("--- %s seconds ---" % (time.time() - start_time))
