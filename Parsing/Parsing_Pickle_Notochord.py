import pickle as pkl
import urllib.request as url
import json

notochord = set()

cell = url.urlopen('https://www.aniseed.cnrs.fr/api/all_genes_by_territory:?cell=B8.6_cell_pair&organism_id=112')
inp = json.load(cell)
for dictionary in inp:
    for what in dictionary['genes']:
        notochord.add(what['gene_model'].split(':')[1])

cell = url.urlopen('https://www.aniseed.cnrs.fr/api/all_genes_by_territory:?cell=A7.3_cell_pair&organism_id=112')
inp = json.load(cell)
for dictionary in inp:
    for what in dictionary['genes']:
        notochord.add(what['gene_model'].split(':')[1])

cell = url.urlopen('https://www.aniseed.cnrs.fr/api/all_genes_by_territory:?cell=A7.7_cell_pair&organism_id=112')
inp = json.load(cell)
for dictionary in inp:
    for what in dictionary['genes']:
        notochord.add(what['gene_model'].split(':')[1])

with open('/home/jeff/Desktop/Thesis/Output/obj/notochord.pkl', 'wb') as fo:
    pkl.dump(notochord, fo, pkl.HIGHEST_PROTOCOL)
