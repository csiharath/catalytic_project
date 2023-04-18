import pickle
from kcat_prediction import kcat_predicton

with open("data/kcat_enzyme.p", "rb") as e:
    enzymes = pickle.load(e)

with open("data/kcat_substrat.p", "rb") as s:
    substrates = pickle.load(s)

with open("data/kcat_product.p", "rb") as p:
    products = pickle.load(p)

df_kcat = kcat_predicton(substrates = substrates,
               products = products,
               enzymes = enzymes)

df_kcat.to_csv("predicted_km.csv")