from read import *
from blast import *
from pprint import pprint as pp

# code goes here
# reader = PdbReader('1ECD.pdb')
# pp(reader.data)
r = FastaReader('eztaxon_qiime_full.fasta', '130648')
seq = r['130648']

print('done')
