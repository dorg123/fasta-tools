from read import *
from time import time
from blast import *
from pprint import pprint as pp

# code goes here
# reader = PdbReader('1ECD.pdb')
# pp(reader.data)
start = time()
r = FastaReader('eztaxon_qiime_full.fasta')
end = time()
print(1000 * (end - start))
seq = r['130648']
start = time()
subseqs = list(subsequences(seq, 18)) + list(subsequences(seq, 19)) + list(subsequences(seq, 20))\
          + list(subsequences(seq, 21)) + list(subsequences(seq, 22))
end = time()
print(1000 * (end - start))
print(len(subseqs))
d = dict()
for subseq in subseqs:
    start = time()
    d[subseq] = sum(1 if subseq in seq else 0 for seq in r.data.values())
    end = time()
    print(1000 * (end - start))
pp(sorted(d.values()))
print('done')
