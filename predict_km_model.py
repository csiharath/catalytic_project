import pickle
import torch
import pandas as pd
from KM_prediction import KM_predicton

model, alphabet = torch.hub.load("facebookresearch/esm:v0.4.0", "esm1b_t33_650M_UR50S")

with open("data/enzyme.p", "rb") as e:
    enzymes = pickle.load(e)

with open("data/substrat.p", "rb") as s:
    substrats = pickle.load(s)

with open("data/compound.p", "rb") as c:
    dict_compounds = pickle.load(c)


df_km = KM_predicton(substrate_list = substrats, enzyme_list= enzymes)
df_km.to_csv("predicted_km.csv")