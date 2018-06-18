codon_to_aa_dict = \
    {'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C', 'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',
     'UUA': 'L', 'UCA': 'S', 'UAA': '*', 'UGA': '*', 'UUG': 'L', 'UCG': 'S', 'UAG': '*', 'UGG': 'W',
     'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R', 'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
     'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', 'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
     'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S', 'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
     'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', 'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
     'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G', 'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
     'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', 'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}
aa_to_codons_dict = \
    {'*': ['UAA', 'UGA', 'UAG'], 'A': ['GCU', 'GCC', 'GCA', 'GCG'], 'C': ['UGU', 'UGC'],
     'D': ['GAU', 'GAC'], 'E': ['GAA', 'GAG'], 'F': ['UUU', 'UUC'], 'G': ['GGU', 'GGC', 'GGA', 'GGG'],
     'H': ['CAU', 'CAC'], 'I': ['AUU', 'AUC', 'AUA'], 'K': ['AAA', 'AAG'],
     'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'], 'M': ['AUG'], 'N': ['AAU', 'AAC'],
     'P': ['CCU', 'CCC', 'CCA', 'CCG'], 'Q': ['CAA', 'CAG'], 'Y': ['UAU', 'UAC'],
     'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
     'T': ['ACU', 'ACC', 'ACA', 'ACG'], 'V': ['GUU', 'GUC', 'GUA', 'GUG'], 'W': ['UGG']}


def compare(*seqs):
    differ_sign = '-'

    def differ(*objects):
        if all(objects[0] == obj for obj in objects):
            return objects[0]
        else:
            return differ_sign

    new_seq = ''.join(map(differ, *seqs)) + differ_sign * (max(map(len, seqs)) - min(map(len, seqs)))
    num = new_seq.count(differ_sign)
    return new_seq, num, 100 - num / max(map(len, seqs)) * 100


def average(*seqs):
    max_len = max(map(len, seqs))
    new_seqs = list(map(lambda seq: list(seq + '-' * (max_len - len(seq))), seqs))
    histogram_list = list(max(histogram(pos), key=lambda x: x[1])[0] for pos in zip(*new_seqs))
    while '-' in histogram_list:
        histogram_list.remove('-')
    return ''.join(histogram_list).rstrip('-')


def mutate(sequence):
    import random

    def mutation(f):
        def wrap(seq):
            seq = list(seq)
            f(seq)
            return ''.join(seq)

        return wrap

    @mutation
    def remove(seq):
        del seq[random.randrange(len(seq))]

    @mutation
    def insert(seq):
        seq.insert(random.randrange(len(seq)), random.choice(nuc))

    @mutation
    def replace(seq):
        y = random.randrange(len(seq))
        seq[y] = random.choice(nuc)

    nuc = sorted(list(set(sequence)))
    z = random.randrange(3)
    return {0: remove, 1: insert, 2: replace}[z](sequence)


def translate(sequence, cut=True):
    start = sequence.find('AUG')
    aas = ''.join(codon_to_aa_dict.get(sequence[i:i + 3], '') for i in range(start, len(sequence), 3))
    return aas[:aas.find('*')] if cut else aas


def complementary(dna_seq):
    switch_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(switch_dict[b] for b in dna_seq)


def transcript(dna_seq):
    return dna_seq.replace('T', 'U')


def reverse_transcript(rna_seq):
    return rna_seq.replace('U', 'T')


def reverse(seq):
    return seq[::-1]


def best_match(s, sub, start=-1, end=-1):
    main = s[start if start != -1 else 0:end if end != -1 else len(s)]
    if len(main) < len(sub):
        return 0, '', 0, 0
    elif len(main) == len(sub):
        a, b, c = compare(main, sub)
        return 0, a, b, c
    else:
        comps = []
        for i in range(len(main) - len(sub)):
            a, b, c = compare(main[i:i + len(sub)], sub)
            comps.append((i, a, b, c))
        return max(comps, key=lambda x: x[3])


def match(s, sub, similarity_range, start=-1, end=-1):
    main = s[start if start != -1 else 0:end if end != -1 else len(s)]
    if len(main) < len(sub):
        return 0, []
    elif len(main) == len(sub):
        a, b, c = compare(main, sub)
        return 1, [(0, a, b, c)] if similarity_range[0] <= c <= similarity_range[1] else 0, []
    else:
        comps = []
        for i in range(len(main) - len(sub)):
            a, b, c = compare(main[i:i + len(sub)], sub)
            comps.append((i, a, b, c))
        comps = list(filter(lambda x: similarity_range[0] <= x[3] <= similarity_range[1], comps))
        return len(comps), comps


def histogram(lst):
    return list((key, lst.count(key)) for key in sorted(set(lst)))


def subsequences(sequence, length):
    if len(sequence) <= length:
        return [sequence] if len(sequence == length) else []
    else:
        for i in range(len(sequence) - length):
            yield sequence[i:i + length]
