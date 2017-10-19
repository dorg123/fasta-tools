from fasta import *
from pprint import pprint as pp

# code goes here
reader = FastaReader('fluSeqs1.txt')
pp(blast(reader['ABN59412|H1N1|1935/0/0|USA'], reader['AAO46268|H2N2|1957/0/0|Chile']))
print('done')
