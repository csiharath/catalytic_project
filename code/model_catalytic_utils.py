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
        try :
            if dataframe_react.loc[reaction.id]['KEGG reaction ID'] != "NaN":
                reaction.notes = {"kegg_id": dataframe_react.loc[reaction.id]['KEGG reaction ID']}
            else:
                print(dataframe_react.loc[reaction.id]['KEGG reaction ID'])
        except :
            print(f"Warning: No KEGG ID found for reaction {reaction.id}")


def add_gene_reaction_rule(model, data:pd.DataFrame):
    """
    """
    reaction2genes = defaultdict(list)

    for i in range(len(data.index)): 
        reaction_name = data.index[i]

        if not isinstance(reaction_name,str) :
            pass

        else:
            reaction_key = reaction_name
        for split_gene_name in data.loc[reaction_name]["Gene name"].split(","):
            if isinstance(split_gene_name, str):

                if "?" in split_gene_name:
                    continue

                elif "/" in split_gene_name:
                    genes = split_gene_name.split("/")
                    # print(genes, type(genes))
                    
                    for gene_name in genes :
                        reaction2genes[reaction_key].append(gene_name)

                elif "(" in split_gene_name :
                    gene = split_gene_name.split("(")[0].strip(" ")
                    reaction2genes[reaction_key].append(gene)
                else:
                    reaction2genes[reaction_key].append(split_gene_name)
            else : continue
    
    for reac_id, genes in reaction2genes.items():
        try:
            reac_obj = model.reactions.get_by_id(reac_id)
            gene_reaction_rule = " or ".join(genes)
            reac_obj.gene_reaction_rule = gene_reaction_rule

        except KeyError:
            continue
    return reaction2genes


def find_compounds_AAseq(model, list_of_enzymes, dict_seq = {}):
    Errors = []
    Success = []

    aaseq = {}

    for enzyme in tqdm(list_of_enzymes):
        aa = False
        ortho = False
        ko = False
        seq = ""
        url = "http://rest.kegg.jp/get/"

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
                                ############# modifier pour verifier les sous unites #############
                                
                                #For each line, we select the ones corresponding to the organism's ID
                                #And store in a list of tuples the KEGG gene IDs and corresponding gene names.
                                #We can then use the gene names to check if the KEGG IDs are the right ones.

                                line_list = line.split(" ")
                                start_index = line_list.index("ECO:")
                                idZ_list = line_list[start_index+1:]
                                idZ_gID = [(idZ, gID) for idZ, gID in zip([elem.split('(')[0] for elem in idZ_list],\
                                                                           [elem.split('(')[1].strip(')') for elem in idZ_list])]
                                
                                #Gene ID comparison between the KEGG response and the model's gene names.
                                for idZ, gID in idZ_gID:
                                    reac_obj = model.reactions.get_by_id(enzyme[0])

                                    if gID.lower() in [str(g.id).lower() for g in reac_obj.genes]:
                                        Success.append(f"\n[>]Found {gID} gene id in reaction {reac_obj.name}. Adding its aa sequence to the dictionary.\n\n===\n")
                                    
                                        # Third request to get AA sequence from CDS ID
                                        r3 = requests.get(url=url+"eco:"+idZ)
                                    

                                        if r3.status_code==200:
                                            info = str(r3.content).split('\'')[1].split('\\n')

                                            for line in info:
                                                if 'AASEQ' in line:
                                                    aa = True
                                                    
                                                elif 'NTSEQ' in line:
                                                    aaseq[seq] = enzyme[0]
                                                    break

                                                elif aa:
                                                    seq += line.split()[0]
                                    else:
                                        continue
                                break

                                    

                if seq != "":
                    dict_seq[enzyme[1]] = {"substrates": sub, "products": prod, "enzyme":seq}
                    break
                elif ko:
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

    return dict_seq, aaseq


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
            compounds = dict_param[enzyme]['substrates'] + dict_param[enzyme]['products']
            for c in compounds:
                km_sub.append(c)
                km_enz.append(dict_param[enzyme]['enzyme'])

    if kcat: 
        for enzyme in dict_param:
            kcat_sub.append(";".join(dict_param[enzyme]['substrates']))
            kcat_prod.append(";".join(dict_param[enzyme]['products']))
            kcat_enz.append(dict_param[enzyme]['enzyme'])

    return km_sub, km_enz, kcat_prod, kcat_sub, kcat_enz


def build_dict_compounds(list_compunds, dict_compounds = {}):
    """
    """
    for compound in list_compunds:
        url = "http://rest.kegg.jp/get/"
        r = requests.get(url=url+compound)
        car = str(r.content)[1]
        info = str(r.content).split(car)[1].split('\\n')
        dict_compounds[compound] = info[1].split()[1].split(';')[0]

    return dict_compounds
