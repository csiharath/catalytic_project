# import cobra
import requests
import pandas as pd
# import numpy
from collections import defaultdict
from tqdm import tqdm


def add_kegg_id_to_model(model, dataframe_react):
    """
    """
    for reaction in model.reactions:
        if 'kegg_id' not in reaction.notes:
            try :
                reaction.notes = {"kegg_id": dataframe_react.loc[dataframe_react['Reaction name'] == reaction.id].iloc[0]['KEGG reaction ID']}
                # if dataframe_react.loc[reaction.id]['KEGG reaction ID'] != "NaN":
                #     reaction.notes = {"kegg_id": dataframe_react.loc[reaction.id]['KEGG reaction ID']}
                # else:
                #     print(dataframe_react.loc[reaction.id]['KEGG reaction ID'])
            except :
                print(f"Warning: No KEGG ID found for reaction {reaction.id}")


def add_gene_reaction_rule(model, data:pd.DataFrame):
    """
    """
    reaction2genes = defaultdict(list)

    for i in range(len(data.index)): 
        reaction_name = data.loc[i]['Reaction name']

        if not isinstance(reaction_name,str) :
            pass
        else:
            reaction_key = reaction_name
            reaction2genes[reaction_key].append(data.loc[i]["Gene name"])

    for reac_id, genes in reaction2genes.items():
        try:
            reac_obj = model.reactions.get_by_id(reac_id)
            gene_reaction_rule = " and ".join(genes)
            reac_obj.gene_reaction_rule = gene_reaction_rule

        except KeyError:
            continue
    return reaction2genes


def find_compounds_AAseq(model, list_of_enzymes):
    """
    model:cobra.core.model.Model object
    list_of_enzymes:list(tuple) -> List of tuples containing model reaciton IDs and corresponding KEGG IDs.
    """
    Errors = []
    Success = []
    dict_seq = {}
    aaseq = {}
    gene_id_dict = defaultdict(dict)
    seq = ""

    for enzyme in tqdm(list_of_enzymes):
        aa = False
        ortho = False
        ko = False
        url = "http://rest.kegg.jp/get/"

        reac_obj = model.reactions.get_by_id(enzyme[0])
        genes = [str(g.id).lower() for g in reac_obj.genes]
        genes_found = []

        # First request to get list of compounds and the KEGG Orthology ID from the KEGG ID
        r1 = requests.get(url=url+enzyme[1])

        if r1.status_code==200:
            # The content is one long string: splitting with \n to iter on every line
            info = str(r1.content).split('\'')[1].split('\\n')

            for line in info:

                # The line with all compounds:
                if 'EQUATION' in line:
                    eq = line.split("=")
                    sub = [c for c in eq[0].split() if 'C' in c]
                    prod = [c for c in eq[1].split() if 'C' in c]

                # Line(s) with KEGG Orthology ID to test
                elif 'ORTHOLOGY' in line or ortho:
                    seq = ""
                    if ortho:
                        id = line.split()[0]
                        if id == 'DBLINKS' or id =='///':
                            ko = True
                            Errors.append(f"\n[!]KEGG ID {enzyme[1]} has no corresponding KEGG Orthology ID for organism eco\n\n===\n")
                            break
                    else:
                        id = line.split()[1]
                        ortho = True

                    # Second request to get CDS ID from KEGG ORthology ID
                    r2 = requests.get(url=url+id)

                    if r2.status_code==200:
                        info = str(r2.content).split('\'')[1].split('\\n')

                        for line in info:

                            # Line corresponding to the organism target (e.coli)
                            if 'ECO:' in line:
                                
                                #For each line, we select the ones corresponding to the organism's ID
                                #And store in a list of tuples the KEGG gene IDs and corresponding gene names.
                                #We can then use the gene names to check if the KEGG IDs are the right ones.

                                line_list = line.split(" ")
                                start_index = line_list.index("ECO:")
                                gID_list = line_list[start_index+1:]
                                gID_gname = [(gID, gname) for gID, gname in zip([elem.split('(')[0] for elem in gID_list],\
                                                                           [elem.split('(')[1].strip(')') for elem in gID_list])]
                                
                                #Gene ID comparison between the KEGG response and the model's gene names.
                                for gID, gname in gID_gname:
                                    if gname.lower() in genes or (len(genes) == 1 and genes_found == []):
                                        Success.append(f"\n[>]Found {gname} gene id in reaction {reac_obj.name}. Adding its aa sequence to the dictionary.\n\n===\n")
                                    
                                        # Third request to get AA sequence from CDS ID
                                        r3 = requests.get(url=url+"eco:"+gID)

                                        if r3.status_code==200:
                                            info = str(r3.content).split('\'')[1].split('\\n')

                                            #Iterating over the lines of the api's response to select
                                            #the ones containing aa sequence.
                                            for line in info:

                                                #First, we look for a gene ID:
                                                if "NCBI-GeneID: " in line:
                                                    gene_id_dict[enzyme[0]]["ncbigene"] = line.split(": ")[1]
                                                elif "UniProt: " in line:
                                                    gene_id_dict[enzyme[0]]["uniprot"] =  line.split(": ")[1]
                                                #Add a check, outside of the loop, for when no gene is found.

                                                if 'AASEQ' in line:
                                                    aa = True
                                                    
                                                elif 'NTSEQ' in line:
                                                    aaseq[seq] = (enzyme[0], gname)
                                                    try:
                                                        genes.remove(gname.lower())
                                                        gene = gname.lower()
                                                    except:
                                                        gene = genes[0]
                                                        Success.append(f"Gene {gname} was associated to {genes[0]}")
                                                    aa = False
                                                    break

                                                elif aa:
                                                    seq += line.split()[0]
                                    else:
                                        genes_found.append(gID)
                                        continue
                                break
                
                

                if seq != "":
                    if enzyme[1] not in dict_seq:
                        dict_seq[enzyme[1]] = [{"substrates": sub, "products": prod, "enzyme":seq, "gene": gene}]
                    else:
                        dict_seq[enzyme[1]].append({"substrates": sub, "products": prod, "enzyme":seq, "gene": gene})
                    # break
                elif ko or genes == []:
                    break
                elif ortho:
                    aa = False
                    #aa=False --> The function will try to find another KEGG orthology ID that contains ECO organism Gene IDS.        

        else: 
            Errors.append(f"[!]KEGG ID {enzyme[1]} is not found in kegg.\n===")

    print(f"\nSUCCESSES:")
    for success_message in Success:
        print(f"{success_message}")

    print(f"\nERRORS:")
    for error_message in Errors:
        print(f"{error_message}")

    return dict_seq, aaseq, gene_id_dict


def create_km_kcat_arguments(dict_param, km = True, kcat = True):
    """
    """
    km_enz = []
    km_sub = []
    kcat_enz = []
    kcat_sub = []
    kcat_prod = []

    if km:
        for enzyme in dict_param:
            for gene in dict_param[enzyme]:
                compounds = gene['substrates'] + gene['products']
                # compounds = dict_param[enzyme]['substrates'] + dict_param[enzyme]['products']
                for c in compounds:
                    km_sub.append(c)
                    km_enz.append(gene['enzyme'])

    if kcat: 
        for enzyme in dict_param:
            for gene in dict_param[enzyme]:
                kcat_sub.append(";".join(gene['substrates']))
                kcat_prod.append(";".join(gene['products']))
                kcat_enz.append(gene['enzyme'])
    

    return km_sub, km_enz, kcat_prod, kcat_sub, kcat_enz


def build_dict_compounds(list_compounds, dict_compounds = {}):
    """
    """
    for compound in list_compounds:
        url = "http://rest.kegg.jp/get/"
        r = requests.get(url=url+compound)
        car = str(r.content)[1]
        info = str(r.content).split(car)[1].split('\\n')
        dict_compounds[compound] = info[1].split(None, 1)[1].split(';')[0]

    return dict_compounds

def add_gene_alt_ids(gene_ids_dict, model):
    for reaction_id, annotations in gene_ids_dict.items():
        reaction_obj = model.reactions.get_by_id(reaction_id)
        reaction_gene_obj = [g for g in reaction_obj.genes]

        if len(reaction_gene_obj) > 0:
            reaction_gene_obj[0].annotation = annotations
