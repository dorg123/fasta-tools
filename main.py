import fasta
from pprint import pprint as pp

david =\
    fasta.DavidReader('DAVID_data_slim.david')
reactome = fasta.SimplifiedTabReader('ReactomePathways.txt')
caption = ''
lst = []
tissues_list = []
reactome_list = []
kegg_list = []
gad_list = []
symbol_list = []
for acc in lst:
    tissues_list += david[acc].get('UP_TISSUE', [])
    reactome_list += david[acc].get('REACTOME_PATHWAY', [])
    kegg_list += david[acc].get('KEGG_PATHWAY', [])
    gad_list += david[acc].get('GAD_DISEASE_CLASS', [])
    symbol_list += david[acc].get('OFFICIAL_GENE_SYMBOL', [])
print(caption)
print('Tissues: ')
pp(sorted(fasta.distribute(tissues_list), key=lambda x: x[1]))
print('Reactome Pathways: ')
pp(sorted(((reactome[p.split(':')[0]], _) for p, _ in fasta.distribute(reactome_list)), key=lambda x: x[1]))
print('Kegg Pathways: ')
pp(sorted(fasta.distribute(kegg_list), key=lambda x: x[1]))
print('GAD Disease Class: ')
pp(sorted(fasta.distribute(gad_list), key=lambda x: x[1]))
print('Gene Symbols: ')
pp(sorted(set(symbol_list)))
print('done')
