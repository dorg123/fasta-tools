import fasta
from pprint import pprint as pp

david =\
    fasta.DavidReader('DAVID_data_slim.david')
reactome = fasta.SimplifiedTabReader('ReactomePathways.txt')
caption = 'hsv2-miR-H20 MIMAT0014702 Herpes Simplex miR-H20'
lst = ['NM_006468', 'NM_001303456', 'NR_024529', 'NR_024528', 'NR_000029', 'NM_001198621', 'NM_020754', 'NR_024615',
       'NR_024241', 'NR_024252', 'NM_014243', 'NM_198281', 'NR_130740', 'NR_130741', 'NR_130742', 'NM_001242412',
       'NM_020731', 'NM_001135575', 'NM_031468', 'NM_001037165', 'NM_001134335', 'NM_014413', 'NR_038120',
       'NM_001017440', 'NM_001083537', 'NR_003494', 'NR_003572', 'NM_001205266', 'NM_004942', 'NM_001286657',
       'NM_001286661', 'NM_152417', 'NM_030780', 'NR_102337', 'NR_102338', 'NM_001135776', 'NM_014007', 'NM_001127898',
       'NM_001127899', 'NM_000084', 'NM_001282163', 'NM_001185075', 'NM_001185076', 'NM_001185081', 'NM_001185082',
       'NM_002024', 'NM_173348', 'NM_002458', 'NR_024249', 'NM_001099653', 'NM_018172', 'NM_152563', 'NR_024254',
       'NM_139071', 'NM_003076', 'NM_001252036', 'NM_001252037', 'NM_002868', 'NR_104647', 'NM_001322238',
       'NM_001322227', 'NM_001322228', 'NM_001322229', 'NM_001322230', 'NM_001322231', 'NM_001322232', 'NM_001322233',
       'NM_001322234', 'NM_001322236', 'NM_004755', 'NM_001322237', 'NM_001322235', 'NM_001286414', 'NM_014444',
       'NM_001289029', 'NM_201400', 'NM_201598', 'NM_024086', 'NR_027620', 'NM_007238', 'NM_183397', 'NM_001205266',
       'NR_135765']
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
pp(sorted(symbol_list))
print('done')
