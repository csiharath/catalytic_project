import pandas as pd
from contextlib import redirect_stdout
from tqdm import tqdm
import sys
import plotly.express as px
from itertools import cycle
import matplotlib.pyplot as plt
import cobra

from cobra.core.model import Model
from cobra.core.reaction import Reaction


def parcours(reaction, flux_dict, metabolites_to_exclude, max_iterations=10000, i=0,v=True, m_l=[], cofactors=set()) :
    """
    if not v :
        orig_stdout = sys.stdout
        f = open('out.txt', 'w')
        sys.stdout = f"""
    for m in reaction.metabolites:
        #print(f"\nChecking out metabolite {m.id}")

        if not m in m_l and not i >= max_iterations: # If the metabolite has not already been visited, and if we haven't reached the maximum number of iterations.
            m_l.append(m) # Adding the metabolite to the list of already visited metabolites.
            if not m.name in metabolites_to_exclude:
                for r in m.reactions:
                    #print(f"\nChecking out reaction {r.id}")

                    if r.id not in flux_dict.keys():

                        # This part checks if there is a flux for the given reaction, and if it is an exchange reaction.
                        # If there is a flux and it is not an exchange reaction, it is added to the flux dict.
                        # If there is a flux and it is an exchange reaction, it is only added to the flux dict and the recursion starts back at the next reaction.
                        # If it was an exchange reaction involving C_x or C_s, it behaves normally.
                        if r.flux != 0.0:
                            #print(f"\nNew reaction, adding {r.id} to flux dict.")
                            flux_dict[r.id] = r.flux
                            flux_dict = parcours(r, flux_dict, metabolites_to_exclude, max_iterations, i, v, m_l, cofactors)
                            i +=1
                        else :
                            #print(f"\nERROR -- flux == 0 for {r.id}")
                            pass

                    else :
                        #print(f"\nERROR -- id in dict for {r.id}")
                        continue
            else :
                continue
        else :
            #print(f"\nERROR -- metabolite {m.id} already visited")
            continue

    """if not v :
        f.close()
        sys.stdout = orig_stdout"""
    return flux_dict

def run_parcours(reaction:Reaction, model:Model, metabolites_to_exclude=["PPi", "CoA", "O2", "H+", "ATP", "CO2", "NADP+", "FADH2", "ADP", "Na+", "NADH", "FAD", "NADPH", "NAD+", "Pi", "H2O"], max_iterations=10000):
    flux_dict = {}
    print(f"\n[{reaction.id}] : Getting all related reactions and fluxes...")
    f = parcours(reaction, flux_dict, metabolites_to_exclude, max_iterations)
    for reaction, flux in f.items() :
        print_reactions(model.reactions.get_by_id(reaction), flux)
        print(f"FLUX : {flux} --- ID : {model.reactions.get_by_id(reaction).id} --- COMPARTMENT : {model.reactions.get_by_id(reaction).compartments}")
        print(f"\nSUBSYSTEM : {model.reactions.get_by_id(reaction).subsystem} --- GENE NAME : {model.reactions.get_by_id(reaction).name}\n\n---\n\n")

def build_reaction_df(optimized_model:Model, by_compartment = True) :
    # Builds dataframes of fluxes associated to compartments or subsystems
    # The dataframes are separated by compartment or subsystem.
    reactions_list = [r for r in optimized_model.reactions]
    if by_compartment :
        compartments_reactions_dict = {}



        for compartment in optimized_model.compartments :
            compartments_reactions_dict[str(compartment)] = {
                "flux" : [abs(r.flux) for r in reactions_list if str(compartment) in r.compartments],\
                "subSystem" : [r.subsystem for r in reactions_list if str(compartment) in r.compartments],\
                "id" : [r.id for r in reactions_list if str(compartment) in r.compartments],\
                "name" : [r.name for r in reactions_list if str(compartment) in r.compartments],\
                "compartment" : [compartment for r in reactions_list if str(compartment) in r.compartments],
                "direction" : [str(r.flux)[0] for r in reactions_list if str(compartment) in r.compartments],\
                "reactants" : [" + ".join([m.name for m in r.reactants]) for r in reactions_list if str(compartment) in r.compartments],\
                "products" : [" + ".join([m.name for m in r.products]) for r in reactions_list if str(compartment) in r.compartments]
            }

        return (pd.DataFrame(compartments_reactions_dict["C_c"]), \
                pd.DataFrame(compartments_reactions_dict["C_r"]), \
                pd.DataFrame(compartments_reactions_dict["C_s"]), \
                pd.DataFrame(compartments_reactions_dict["C_m"]), \
                pd.DataFrame(compartments_reactions_dict["C_p"]), \
                pd.DataFrame(compartments_reactions_dict["C_x"]), \
                pd.DataFrame(compartments_reactions_dict["C_l"]), \
                pd.DataFrame(compartments_reactions_dict["C_g"]), \
                pd.DataFrame(compartments_reactions_dict["C_n"]))
    else :

        subsystem_reactions_dict = {}
        
        for subsystem in optimized_model.groups :
            # Exception for beta-oxydation and Carnitine shuttle.
            # As there is a lot of different beta oxydation and carnitine shuttle subsystems,
            # they are grouped as one.
            if "beta oxidation of" in str(subsystem).lower() or "carnitine shuttle" in str(subsystem).lower() :
                subsystem_key = " ".join(str(subsystem).lower().split(" ")[:2]).lower()
            else :
                subsystem_key = str(subsystem).lower()

            subsystem_reactions_dict[subsystem_key] = {
                "flux" : [abs(r.flux) for r in reactions_list if str(subsystem) in r.subsystem],\
                "subSystem" : [str(r.subsystem).lower() for r in reactions_list if str(subsystem) in r.subsystem],\
                "id" : [r.id for r in reactions_list if str(subsystem) in r.subsystem],\
                "name" : [r.name for r in reactions_list if str(subsystem) in r.subsystem],\
                "compartment" : [str([comp for comp in r.compartments][0]) for r in reactions_list if str(subsystem) in r.subsystem],\
                "direction" : [str(r.flux)[0] for r in reactions_list if str(subsystem) in r.subsystem],\
                "reactants" : [" + ".join([m.name for m in r.reactants]) for r in reactions_list if str(subsystem) in r.subsystem],\
                "products" : [" + ".join([m.name for m in r.products]) for r in reactions_list if str(subsystem) in r.subsystem]
            }

        return subsystem_reactions_dict

def merge_dict_keys(dict_to_merge, key1, key2):
    data_to_add = dict_to_merge.pop(key2)

    if isinstance(dict_to_merge, dict):
        for key, value in data_to_add.items():
            if isinstance(value,list):
                for v in value:
                    dict_to_merge[key1][key].append(v)
    return dict_to_merge


def get_subsystem_fluxes(dfs, multiple=True):
    subsystem_fluxes = {}

    if multiple :
        for df in dfs :
            for line in df.iterrows():
                #print(line[1]["flux"])
                if "Transport" not in line[1]["subSystem"] and "Exchange" not in line[1]["subSystem"] and len(line[1]["subSystem"]) > 1:
                    try :
                        subsystem_fluxes[line[1]["subSystem"]] += abs(line[1]["flux"])
                    except KeyError :
                        subsystem_fluxes[line[1]["subSystem"]] = abs(line[1]["flux"])
            df_bar = pd.DataFrame(subsystem_fluxes,index=["Fluxes"]).T
        return df_bar
    else :
        for line in dfs.iterrows():
            #print(line[1]["flux"])
            if "Transport" not in line[1]["subSystem"] and "Exchange" not in line[1]["subSystem"] and len(line[1]["subSystem"]) > 1:
                try :
                    subsystem_fluxes[line[1]["subSystem"]] += abs(line[1]["flux"])
                except KeyError :
                    subsystem_fluxes[line[1]["subSystem"]] = abs(line[1]["flux"])
        df_bar = pd.DataFrame(subsystem_fluxes,index=["Fluxes"]).T
        return df_bar


def reformat_dict(origin_dict, by, data, subsystem):
    """
    Reformats a dictionary to make it contain "data" relative to "by"
    """
    new_dict = {}

    for prot_names, flux in zip(origin_dict[subsystem][by], origin_dict[subsystem][data]):
        for name in prot_names.split('/'):
            new_dict[name] = new_dict.get(name, 0.0) + flux  
    
    return new_dict


def combine_fluxes_by_gene(reaction_dict_1, reaction_dict_2):
    """
    Function taking as input two dictionaries containing fluxes associated with gene names.
    Returns a dictionary of lists used to plot a comparison of flux attribution to proteins between two models.
    The input dictionaries should follow the format :
    {protein_id : flux}
    """
    common_index = [] #Building a list containing every gene mentionned at least once in either dictionaries.

    for name_1 in reaction_dict_1.keys(): #Initializing the common_index list with all gene names of reaction_dict_1.
        if name_1 in common_index:
            continue
        else:
            common_index.append(name_1) #Adding the name to the index if does not already contain it.

    for name_2 in reaction_dict_2.keys():
            if name_2 in common_index:
                continue
            else:
                common_index.append(name_2)
    #At this point, the common_index list should be filled.

    fluxes_list_1 = []
    fluxes_list_2 = []
    for protein_name in common_index:
        try:
            fluxes_list_1.append(float(reaction_dict_1[protein_name]))
        except KeyError:
            #print(f"\nProtein {protein_name} not found in dict_1")
            fluxes_list_1.append(0.0)
        
        try:
            fluxes_list_2.append(float(reaction_dict_2[protein_name]))
        except KeyError:
            #print(f"\nProtein {protein_name} not found in dict_1")
            fluxes_list_2.append(0.0)
    
    #At this point, there should be 3 lists, one of protein ids and two of fluxes.
    #This next part checks it :
    #print(f"\nlen(names) = {len(common_index)} \t--\t len(fluxes_1) = {len(fluxes_list_1)} \t--\t len(fluxes_2) = {len(fluxes_list_2)}.")
    if len(common_index) == len(fluxes_list_1) and len(common_index) == len(fluxes_list_2) :
        final_dict = {"names" : common_index, "fluxes_1" : fluxes_list_1, "fluxes_2" : fluxes_list_2}

        return pd.DataFrame(final_dict)
    else :
        raise ValueError
### DATA VISUALISATION ###

def plot_treemap(df, model, title, path=['subSystem', 'id'], flux_filter=0.0, color_by = "subsystem") :
    ### Building colormap :


    if color_by == "subsystem" :
        # cols = ['#E48F72','#FC6955','#7E7DCD','#BC7196','#86CE00','#E3EE9E','#22FFA7','#FF0092','#C9FBE5','#B68E00','#00B5F7','#6E899C',
        # '#D626FF','#AF0038','#0D2A63','#6C4516','#DA60CA','#1616A7','#620042','#A777F1','#862A16','#778AAE','#6C7C32','#B2828D',
        # '#FC0080','#00A08B','#511CFB','#EB663B','#750D86','#B68100','#222A2A','#FB0D0D','#1CA71C','#E15F99',
        # '#2E91E5','#DC587D','#EEA6FB','#479B55','#FF9616','#F6F926','#0DF9FF','#FE00CE','#FED4C4','#6A76FC','#00FE35','#FD3216',
        # '#2E91E5','#E15F99','#1CA71C','#FB0D0D','#222A2A','#B68100','#750D86','#EB663B','#511CFB','#00A08B','#FB00D1',
        # '#FC0080','#B2828D','#6C7C32','#778AAE','#862A16','#A777F1','#620042','#1616A7','#DA60CA','#6C4516','#0D2A63','#AF0038']
        cols = ["#4D455D", "#E96479", "#F5E9CF", "#7DB9B6", "#539165", "#820000", "#FFB100"]
        # Used to plot fluxes compartment by compartment --> coloring by the next highest hierarchical category = subsystems.


        subsystems = set()
        for r in model.reactions :
            if len(r.subsystem) >0 and not "Transport" in r.subsystem and not "Exchange" in r.subsystem :
                subsystems.add(r.subsystem)

        cmap = {}
        for subsystem, color in zip(subsystems, cycle(set(cols))) :
            cmap[subsystem] = color
        ###

        df = df.loc[(df["subSystem"] != "Transport, mitochondrial")& (df["subSystem"] != "Transport, extracellular" ) & (df["subSystem"] != "Exchange reactions")& (df["subSystem"] != "Transport, peroxisomal" )]
        fig = px.treemap(df.loc[(df["flux"] >= flux_filter) & (df["name"] != "Null")].to_dict() , path=path,
                    values='flux', hover_name= "name" ,color='subSystem',hover_data = ["direction", "reactants", "products"], color_discrete_map=cmap)

        fig.update_layout(title_text=title, font_size=12)

    elif color_by == "compartment" :
        cols = ["#4D455D", "#E96479", "#F5E9CF", "#7DB9B6", "#539165", "#820000", "#FFB100"]
        # Used to plot fluxes subsystem by subsystem --> coloring by the next highest hierarchical category = compartments.

        compartments = set()
        for r in model.reactions :
            if len(r.compartments) > 0 and not "C_r" in r.compartments and not "C_l" in r.compartments and not "C_g" in r.compartments :
                for comp in r.compartments :
                    compartments.add(comp)
        cmap =  {}
        for compartment, color in zip(compartments, cycle(set(cols))) :
            cmap[compartment] = color
        ###

        df = df.loc[(df["subSystem"] != "Transport, mitochondrial")& (df["subSystem"] != "Transport, extracellular" ) & (df["subSystem"] != "Exchange reactions")& (df["subSystem"] != "Transport, peroxisomal" )]
        fig = px.treemap(df.loc[(df["flux"] >= flux_filter) & (df["name"] != "Null")].to_dict() , path=path,
                    values='flux', hover_name= "name" ,color='compartment', hover_data = ["direction", "reactants", "products"], color_discrete_map=cmap)
    return fig


def compartment_fluxes_barplots(model_1, model_2) :
    barplots = {}
    for compartments_iHep, compartments_G2 in zip(build_reaction_df(model_1), build_reaction_df(model_2)): # -> Iterates over the different compartments fluxes dataframes.

        color_HepG2 = 'salmon'
        color_iHep = 'lightblue'
        subS_HepG2 = get_subsystem_fluxes(compartments_G2, multiple=False)
        subS_iHep = get_subsystem_fluxes(compartments_iHep, multiple=False)
        df_both = pd.concat([subS_iHep, subS_HepG2],axis=1)
        sum_HepG2 = subS_HepG2.sum(0)
        sum_iHep = subS_iHep.sum(0)
        df_both.columns = ["iHep", "HepG2"]

        df_both = df_both.loc[(df_both["iHep"] != 0.0) & (df_both["HepG2"] != 0.0)]

        compartment_G2 = compartments_G2["compartment"][0]
        compartment_iHep = compartments_iHep["compartment"][0]



        if not compartment_iHep == compartment_G2 :
            print("Erreur compartements")
            break
        if not "C_r" in compartment_iHep and not "C_s" in compartment_iHep and not "C_x" in compartment_iHep and not "C_g" in compartment_iHep and not "C_l" in compartment_iHep :

            labels = list(df_both.index)
            normalized_fluxes_iHep = [float(val/sum_iHep) for val in df_both.loc[:,"iHep"].values]
            normalized_fluxes_G2 = [float(val/sum_HepG2) for val in df_both.loc[:,"HepG2"].values]

            xmax = max(max(normalized_fluxes_G2), max(normalized_fluxes_iHep))
            xmin = min(min(normalized_fluxes_G2), min(normalized_fluxes_iHep))

            fig, axes = plt.subplots(figsize=(5,7), ncols=2, sharey=True)
            fig.tight_layout()


            axes[0].barh(labels, normalized_fluxes_iHep, align='center', color=color_iHep, zorder=10)
            axes[0].set_title(f"iHep : {compartment_iHep}")
            axes[1].barh(labels, normalized_fluxes_G2 , align='center', color=color_HepG2, zorder=10)
            axes[1].set_title(f"Hep_G2 : {compartment_G2}")
            axes[0].set_xlabel("Fraction des comptages totaux")
            axes[1].set_xlabel("Fraction des comptages totaux")
            axes[0].invert_yaxis() # labels read top-to-bottom
            axes[0].invert_xaxis() # mirror data for both duildings
            plt.xlim = (xmin, xmax)
            plt.close()
            barplots[compartment_G2] = (fig, axes)
    return barplots


def subsystem_barplots(model_1:Model, model_2:Model, model_1_name:str, model_2_name:str, subsystems_dict:dict) :
    # Main model, used for labels is model_1 !
    barplots = {}
    # Take all reaction-associated fluxes, and returns them as a dictionary. 
    # by_compartment being False means that the dictionary's keys are subSystems.
    fluxes_by_subsystem_dict_1 = build_reaction_df(model_1, by_compartment=False)
    fluxes_by_subsystem_dict_2 = build_reaction_df(model_2, by_compartment=False)
    color_1 = "wheat"
    color_2 = "powderblue"

    for subsystem in subsystems_dict.keys() :
        # Take the fluxes_by_subsystems dictionaries and shrinks them, selecting only the fluxes and name 
        # columns and only one subsystem.
        flux_by_prot_dict_1 = reformat_dict(fluxes_by_subsystem_dict_1, "name", "flux", subsystem)
        flux_by_prot_dict_2 = reformat_dict(fluxes_by_subsystem_dict_2, "name", "flux", subsystem)

        # Combines the two dictionaries to make a final dataframe containing, for each gene,
        # its associated flux, in each model.
        df_flux = combine_fluxes_by_gene(flux_by_prot_dict_1, flux_by_prot_dict_2)
        df_flux_nonull = df_flux.loc[(df_flux["fluxes_1"] > 0.0) | (df_flux["fluxes_2"] > 0.0)]
        labels = list(df_flux_nonull["names"])
        try :

            xmax = max(max(df_flux_nonull["fluxes_1"]),max(df_flux_nonull["fluxes_2"]))
            xmin = min(min(df_flux_nonull["fluxes_1"]), min(df_flux_nonull["fluxes_2"]))
        except ValueError:
            print(f"\nNo fluxes associated with subsystem {subsystem} in either models.")
            continue
        edgecolors_1 = []
        edgecolors_2 = []
        for f1, f2 in zip(df_flux_nonull["fluxes_1"].tolist(), df_flux_nonull["fluxes_2"].tolist()):
            if float(f1) > float(f2):
                edgecolors_1.append("firebrick")
                edgecolors_2.append("powderblue")
            elif float(f1) < float(f2):
                edgecolors_2.append("firebrick")
                edgecolors_1.append("wheat")
            else :
                edgecolors_1.append("wheat")
                edgecolors_2.append("powderblue")
        fig, axes = plt.subplots(figsize=subsystems_dict[subsystem], ncols=2, sharey=True, facecolor='w')
        fig.tight_layout()
        fig.suptitle(f"{subsystem} protein fluxes.")
        axes[0].barh(labels, df_flux_nonull["fluxes_1"], align='center', color=color_1, edgecolor = edgecolors_1, linewidth=3,zorder=10)
        axes[0].set_title(f"{model_1_name}")
        axes[1].barh(labels, df_flux_nonull["fluxes_2"] , align='center', color=color_2, edgecolor = edgecolors_2, linewidth=3,zorder=10)
        axes[1].set_title(f"{model_2_name}")
        axes[0].set_xlabel("Fraction des comptages totaux")
        axes[1].set_xlabel("Fraction des comptages totaux")
        axes[0].invert_yaxis() # labels read top-to-bottom
        axes[0].invert_xaxis() # mirror data for both duildings


        ticks_0 = axes[0].get_xticks()
        ticks_1 = axes[1].get_xticks()
        max_0 = max(ticks_0)
        max_1 = max(ticks_1)
        if max_0 > max_1 :
            axes[1].set_xticks(ticks_0)
        else :
            axes[0].set_xticks(ticks_1)
            
        plt.xlim = (0.0, max(max_0, max_1))
        plt.close()
        barplots[subsystem] = (fig, axes)


    return barplots




def print_exchanges(optimized_model, filter ) :
    intakes = []
    secretions = []
    neutrals = []
    if filter == "non_null" :
        flux_comparison = 0.0
    elif filter == "all" :
        flux_comparison = 200000.0
    for reaction in optimized_model.boundary :
        if reaction.flux != flux_comparison :


            # Jolification
            spaces = 12
            spaces_str = ""
            for i in str(round(reaction.flux)) :
                spaces -= 1
            if "EX_t" in reaction.id :
                spaces -= 1
            for j in range(spaces) :
                spaces_str += " "
            # Fin de la jolification

            ml = [metab for metab in reaction.metabolites]
            m = [metab for metab in reaction.metabolites][0]
            if reaction.flux > 0.0 :
                secretions.append(f"{reaction.id} : {round(reaction.flux)}{spaces_str}ub : {reaction.upper_bound}\t---\t\tmetabolites : \t id : {m.id} --- metabolite name : {m.name} ; id : {m.id}")
            elif reaction.flux < 0.0 :
                intakes.append(f"{reaction.id} : {round(reaction.flux)}{spaces_str}ub : {reaction.upper_bound}\t---\t\tmetabolites : \t id : {m.id} --- metabolite name : {m.name} ; id : {m.id}")
            else :
                neutrals.append(f"{reaction.id} : {round(reaction.flux)}{spaces_str}ub : {reaction.upper_bound}\t---\t\tmetabolites : \t id : {m.id} --- metabolite name : {m.name} ; id : {m.id}")

    intakes.sort(key=lambda f : float(f.split(": ")[1].split("ub")[0]))
    secretions.sort(key=lambda f : float(f.split(": ")[1].split("ub")[0]), reverse=True)
    print("\n##########\nINTAKES :\n")
    for i in intakes :
        print(i)
    print("\n##########\nSECRETIONS :\n")
    for s in secretions :
        print(s)
    print("\n##########\nNEUTRALS : \n")
    for n in neutrals :
        print(n)


def print_reactions(reaction, flux = 0.0, v=True):
    '''
    Affiche les reactions d'un modele de façon plus lisible
    '''
    list_met = []
    liste_reactif = [i for i in reaction.reactants]
    liste_produit = [i for i in reaction.products]
    string = ""
    string2 = ""
    if not "EX_" in reaction.id :
        if flux > 0.0 :
                fleche = "-->"
        elif flux < 0.0 :
                fleche = "<--"
        else:
                fleche = "<=>"
    else :
        if flux > 0.0 :
            fleche = "<--"
        elif flux < 0.0 :
             fleche = "-->"
        else :
             fleche = "<=>"

    for i in liste_reactif:
        if i != liste_reactif[-1]:
            list_met.append(i.name)
            string += str(float(abs(reaction.metabolites[i])))
            string += " "
            string += str(i.name)
            string += " + "
            list_met.append("+")

            string2 += str(float(abs(reaction.metabolites[i])))
            string2 += " "
            string2 += str(i.id)
            string2 += " + "

        else:
            list_met.append(i.name)
            string += str(float(abs(reaction.metabolites[i])))
            string += " "
            string += str(i.name)
            string += " "

            string2 += str(float(abs(reaction.metabolites[i])))
            string2 += " "
            string2 += str(i.id)
            string2 += " "

    list_met.append(fleche)
    string += fleche
    string += " "

    string2 += fleche
    string2 += " "


    for i in liste_produit:
        if i != liste_produit[-1]:
            list_met.append(i.name)
            string += str(float(abs(reaction.metabolites[i])))
            string += " "
            string += str(i.name)
            string += " + "
            list_met.append("+")

            string2 += str(float(abs(reaction.metabolites[i])))
            string2 += " "
            string2 += str(i.id)
            string2 += " + "

        else:
            list_met.append(i.name)
            string += str(float(abs(reaction.metabolites[i])))
            string += " "
            string += str(i.name)

            string2 += str(float(abs(reaction.metabolites[i])))
            string2 += " "
            string2 += str(i.id)
    if v :
        print(string)
        print("")
        print(string2)
    else :
        return(string, string2, list_met)
