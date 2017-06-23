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


class FileReader:
    def __init__(self, filename):
        self._filename = filename
        self._read()

    @property
    def filename(self):
        return self._filename

    @property
    def data(self):
        return self._data


class FastaReader(FileReader):
    def _read(self):
        with open(self._filename, 'r') as f:
            lines = list(line.rstrip('\n') for line in f.readlines())
        info = None
        self._data = dict()
        for line in lines:
            if line.startswith('>'):
                info = line.lstrip('>')
                self._data[info] = ''
            else:
                self._data[info] += line.upper()

    @staticmethod
    def read_head(head, raw=False):
        parts = head.split('|')
        if raw:
            return dict(zip(range(len(parts)), parts))

        def find_all(s, sub):
            index = [s.find(sub)]
            while index[-1] != -1:
                index.append(s.find(sub, index[-1] + 1))
            index.pop()
            return index

        last = parts[-1]
        d = dict((parts[i], parts[i + 1]) for i in range(0, len(parts) - 1, 2))
        if '=' in last:
            separators = [0] + list(last.rfind(' ', 0, i) for i in find_all(last, '=')) + [len(last)]
            parts = list(last[separators[i - 1]:separators[i]].strip() for i in range(1, len(separators)))
            d.update(dict((part[:part.find('=')], part.split('=')[1]) if '=' in part else (0, part) for part in parts))
        return d

    def __str__(self):
        return 'fasta.FastaReader: {}, {} entries' \
            .format(self.filename, len(self.data))

    def __repr__(self):
        return 'fasta.FastaReader({})<length: {}>' \
            .format(self.filename, len(self.data))


class TabReader(FileReader):
    @property
    def head(self):
        return self._head

    def _read(self):
        with open(self._filename, 'r') as f:
            lines = f.readlines()
        it = iter(lines)
        self._head = next(it).rstrip('\n').split('\t')
        self._data = list(dict(zip(self._head, line.rstrip('\n').split('\t'))) for line in it)

    def __str__(self):
        return 'fasta.TabReader: {} (head: {}), {} entries' \
            .format(self.filename, ', '.join(self.head), len(self.data))

    def __repr__(self):
        return 'fasta.TabReader({})<length: {}; head: {}>' \
            .format(self.filename, len(self.data), self.head)


class FastqReader(FileReader):
    def _read(self):
        with open(self.filename, 'r') as f:
            lines = list(line.rstrip('\n') for line in f.readlines())
        self._data = dict()
        head = ''
        reading = False
        for line in lines:
            if line.startswith('@'):
                head = line.lstrip('@')
                reading = True
                self._data[head] = '', ''
            elif line.startswith('+'):
                reading = False
            else:
                s, a = self._data[head]
                if reading:
                    s += line
                else:
                    a += line
                self._data[head] = s, a

    def __str__(self):
        return 'fasta.FastqReader: {}, {} entries' \
            .format(self.filename, len(self.data))

    def __repr__(self):
        return 'fasta.FastqReader({})<length: {}>' \
            .format(self.filename, len(self.data))


class GenBankReader(FileReader):
    def _read(self):
        with open(self.filename, 'r') as f:
            lines = list(line.strip('\n').strip('\n').strip(' ') for line in f.readlines())
        self._data = []
        obj = GenBankRecord()

        def is_cap(cap_in):
            ca = cap_in
            while '  ' in ca:
                ca = ca.replace('  ', ' ')
            if ca.count(' ') != 1:
                return False
            ps = ca.split(' ')
            s = 0
            g1 = list(str(i) for i in range(10))
            for cc in ps[1]:
                if s == 0:
                    if cc in ['.', '>', '<']:
                        s = 1
                    elif cc not in g1:
                        return None
                if s == 1:
                    if cc in g1:
                        s = 2
                    elif cc not in ['.', '>', '<']:
                        return None
                if s == 2:
                    if cc not in g1:
                        return None
            return ca
        state = ''
        for line in lines:
            if line.startswith('LOCUS'):
                obj.locus = line[5:].strip()
            elif line.startswith('DEFINITION') or state == 'DEFINITION':
                if state == '':
                    definition = line[10:].strip()
                    state = 'DEFINITION'
                    continue
                if line.startswith('ACCESSION'):
                    obj.definition = definition
                    state = ''
                else:
                    definition += ' ' + line
            if line.startswith('ACCESSION'):
                obj.accession = line[9:].strip()
            elif line.startswith('VERSION'):
                obj.version = line[7:].strip()
            elif line.startswith('SOURCE'):
                obj.source = line[6:].strip()
            elif line.startswith('ORGANISM') or state == 'ORGANISM':
                if state == '':
                    organism = line[8:].strip()
                    taxonomy = ''
                    state = 'ORGANISM'
                    continue
                if line.endswith('.'):
                    taxonomy += ' ' + line.strip('.') + '; ' + organism
                    obj.taxonomy = taxonomy.strip().split('; ')
                    state = ''
                else:
                    taxonomy += ' ' + line
            elif line.startswith('PUBMED'):
                obj.pubmed.append(line[6:].strip())
            elif line.startswith('COMMENT') or state == 'COMMENT':
                if state == '':
                    com = line[7:].strip()
                    state = 'COMMENT'
                    continue
                if len(line) > 0 and all(cc.isupper() for cc in line[:line.find(' ')]):
                    obj.comment = com
                    state = ''
                else:
                    com += '\n' + line
            if line.startswith('FEATURES') or state == 'FEATURES':
                if state == '':
                    features = dict()
                    state = 'FEATURES'
                    continue
                if line.startswith('ORIGIN'):
                    obj.features = features
                    state = ''
                else:
                    c = is_cap(line)
                    if c:
                        cap = c
                        features[cap] = []
                    elif line.startswith('/'):
                        eq_sign = line.find('=')
                        features[cap].append((line[:eq_sign].strip('/'), line[eq_sign + 1:].strip('"')))
                    else:
                        if len(features[cap]) == 0:
                            cap = line
                            features[cap] = []
                        else:
                            a, b = features[cap][-1]
                            features[cap][-1] = a, b + line.strip('"')
            if line.startswith('ORIGIN') or state == 'ORIGIN':
                if state == '':
                    origin = ''
                    state = 'ORIGIN'
                    continue
                if line == '//':
                    obj.sequence = origin
                    state = ''
                else:
                    origin += line[line.find(' '):].replace(' ', '').upper()
            if line.startswith('//'):
                self._data.append(obj)
                obj = GenBankRecord()

    def __str__(self):
        return 'fasta.GenBankReader: {}, {} entries' \
            .format(self.filename, len(self.data))

    def __repr__(self):
        return 'fasta.GenBankReader({})<length: {}>' \
            .format(self.filename, len(self.data))


class DavidReader(FileReader):
    def _read(self):
        with open(self._filename, 'r') as f:
            lines = list(l.strip('\n') for l in f.readlines())
        flag = True
        self._data = list()
        self._accessions = dict()
        index = -1
        for line in lines:
            if flag:
                index += 1
                caption, _, data = line.partition('\t')
                self._data.append({caption: data})
                self._accessions.update((acc, index) for acc in caption.split(', '))
                flag = False
            elif line == '':
                flag = True
            else:
                title, _, data = line.partition('\t')
                self._data[index][title] = data

    def get_entry(self, accession):
        try:
            return self._data[self._accessions[accession]]
        except KeyError:
            return dict()

    def __str__(self):
        return 'fasta.DavidReader: {}, {} entries' \
            .format(self.filename, len(self.data))

    def __repr__(self):
        return 'fasta.DavidReader({})<length: {}>' \
            .format(self.filename, len(self.data))


class GenBankRecord:
    def __init__(self):
        self._locus = ''
        self._definition = ''
        self._accession = ''
        self._version = ''
        self._source = ''
        self._taxonomy = []
        self._pubmed = []
        self._comment = ''
        self._features = dict()
        self._sequence = ''

    @property
    def locus(self):
        return self._locus

    @locus.setter
    def locus(self, locus):
        self._locus = locus

    @property
    def definition(self):
        return self._definition

    @definition.setter
    def definition(self, definition):
        self._definition = definition

    @property
    def accession(self):
        return self._accession

    @accession.setter
    def accession(self, accession):
        self._accession = accession

    @property
    def version(self):
        return self._version

    @version.setter
    def version(self, version):
        self._version = version

    @property
    def source(self):
        return self._source

    @source.setter
    def source(self, source):
        self._source = source

    @property
    def taxonomy(self):
        return self._taxonomy

    @taxonomy.setter
    def taxonomy(self, taxonomy):
        self._taxonomy = taxonomy

    @property
    def pubmed(self):
        return self._pubmed

    @pubmed.setter
    def pubmed(self, pubmed):
        self._pubmed = pubmed

    @property
    def comment(self):
        return self._comment

    @comment.setter
    def comment(self, comment):
        self._comment = comment

    @property
    def features(self):
        return self._features

    @features.setter
    def features(self, features):
        self._features = features

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, sequence):
        self._sequence = sequence

    def __str__(self):
        return \
            'GenBank Record {}\n{}\n\tLocus:\t\t{}\n\tDefinition:\t{}\n\tSource:\t\t{}\n\t\t\t\t{}'\
            '\n\tPubMed:\t\t{}\n\tComment:\t{}\n\tFeatures:\t{}\n\tSequence:\t{}'\
            .format(self.version, (15 + len(self.version)) * '-', self.locus, self.definition, self.source,
                    '\n\t\t\t\t'.join('; '.join(self.taxonomy[i:i + 5]) for i in range(0, len(self.taxonomy), 5)),
                    ', '.join(self.pubmed), self.comment.replace('\n', '\n\t\t\t\t'),
                    '\n\t\t\t\t'.join(head + '\n\t\t\t\t\t' + '\n\t\t\t\t\t'
                                .join(p1 + '="' + p2 + '"' for p1, p2 in body) for head, body in self.features.items()),
                    '\n\t\t\t\t'.join(self.sequence[i:i + 90] for i in range(0, len(self.sequence), 90)))

    def __repr__(self):
        return 'fasta.GenBankRecord()<accession = {}; source = {}; sequence = {}>'\
            .format(self.version, self.source, len(self.sequence))


class FileWriter:
    def __init__(self, filename, open_file=False):
        self._open = False
        self._filename = filename
        self._file = None
        if open_file:
            self.open()

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def open(self):
        self._open = True
        self._file = open(self._filename, 'w')

    def close(self):
        try:
            self._file.close()
        except AttributeError:
            raise ValueError('File is not open.')
        finally:
            self._open = False

    @property
    def filename(self):
        return self._filename

    @property
    def is_open(self):
        return self._open


class FastaWriter(FileWriter):
    def mass_write(self, entries, lim=60):
        try:
            for head, body in entries:
                self.write(head, body, lim)
        except ValueError:
            raise ValueError("The given input is either not an iterable or its elements aren't all 2-part tuples.")

    def write(self, head, body, lim=60):
        self._file.write('>{}\n'.format(head))
        l = len(body)
        for i in range(0, l, lim):
            if l - i > lim:
                self._file.write('{}\n'.format(body[i:i + lim]))
            else:
                self._file.write('{}\n'.format(body[i:]))

    def __str__(self):
        return 'fasta.FastaWriter: {}, currently {}' \
            .format(self.filename, 'open' if self.is_open else 'close')

    def __repr__(self):
        return 'fasta.FastaWriter({})<{}>' \
            .format(self.filename, 'open' if self.is_open else 'close')


class TabWriter(FileWriter):
    def __init__(self, filename, head, open_file=False):
        self._head = head
        super().__init__(filename, open_file)

    @property
    def head(self):
        return self._head

    def open(self):
        super().open()
        self._file.write('\t'.join(self._head) + '\n')

    def write(self, line):
        self._file.write('\t'.join(line[col] for col in self._head) + '\n')

    def mass_write(self, *lines):
        for line in lines:
            self.write(line)

    def __str__(self):
        return 'fasta.TabWriter: {} ({}), currently {}' \
            .format(self.filename, ', '.join(self.head), 'open' if self.is_open else 'close')

    def __repr__(self):
        return 'fasta.TabWriter({}, {})<{}>' \
            .format(self.filename, self.head, 'open' if self.is_open else 'close')


class FastqWriter(FileWriter):
    def mass_write(self, entries, lim=60):
        try:
            for head, seq, annotation in entries:
                self.write(head, seq, annotation, lim)
        except ValueError:
            raise ValueError("The given input is either not an iterable or its elements aren't all 3-part tuples.")

    def write(self, head, seq, annotation, lim=60):
        self._file.write('@{}\n'.format(head))
        l = len(seq)
        for i in range(0, l, lim):
            if l - i > lim:
                self._file.write('{}\n'.format(seq[i:i + lim]))
            else:
                self._file.write('{}\n'.format(seq[i:]))
        self._file.write('+\n{}\n'.format(annotation))

    def __str__(self):
        return 'fasta.FastqWriter: {}, currently {}' \
            .format(self.filename, 'open' if self.is_open else 'close')

    def __repr__(self):
        return 'fasta.FastqWriter({})<{}>' \
            .format(self.filename, 'open' if self.is_open else 'close')


def compare_sequences(*seqs):
    differ_sign = '-'

    def differ(*objects):
        if all(objects[0] == obj for obj in objects):
            return objects[0]
        else:
            return differ_sign

    new_seq = ''.join(map(differ, *seqs)) + differ_sign * (max(map(len, seqs)) - min(map(len, seqs)))
    num = new_seq.count(differ_sign)
    return new_seq, num, 100 - num / max(map(len, seqs)) * 100


def get_segments(sequence):
    dash = True
    last = 0
    for i in range(len(sequence)):
        if (sequence[i] == '-') == dash:
            yield sequence[last:i]
            last = i
            dash = not dash


def average_sequence(*seqs):
    max_len = max(map(len, seqs))
    new_seqs = list(map(lambda seq: list(seq + '-' * (max_len - len(seq))), seqs))
    distribution_list = list(max(list((aa, pos.count(aa)) for aa in set(pos)),
                                 key=lambda x: x[1])[0] for pos in zip(*new_seqs))
    return ''.join(distribution_list).rstrip('-')


def create_mass_sequences(seqs, func):
    return list(map(func, seqs))


def mutate_sequence(sequence):
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


def translate_sequence(sequence, cut=True):
    start = sequence.find('AUG')
    aas = ''.join(codon_to_aa_dict.get(sequence[i:i + 3], '') for i in range(start, len(sequence), 3))
    return aas[:aas.find('*')] if cut else aas


def complementary_dna(dna_seq):
    switch_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(switch_dict[b] for b in dna_seq)


def transcript(dna_seq):
    return dna_seq.replace('T', 'U')


def reverse_transcript(rna_seq):
    return rna_seq.replace('U', 'T')


def reverse_sequence(seq):
    return seq[::-1]


def find_best_match(s, sub, start=-1, end=-1):
    main = s[start if start != -1 else 0:end if end != -1 else len(s)]
    if len(main) < len(sub):
        return 0, '', 0, 0
    elif len(main) == len(sub):
        a, b, c = compare_sequences(main, sub)
        return 0, a, b, c
    else:
        comps = []
        for i in range(len(main) - len(sub)):
            a, b, c = compare_sequences(main[i:i + len(sub)], sub)
            comps.append((i, a, b, c))
        return max(comps, key=lambda x: x[3])


def count_matches(s, sub, similarity_range, start=-1, end=-1):
    main = s[start if start != -1 else 0:end if end != -1 else len(s)]
    if len(main) < len(sub):
        return 0, []
    elif len(main) == len(sub):
        a, b, c = compare_sequences(main, sub)
        return 1, [(0, a, b, c)] if similarity_range[0] <= c <= similarity_range[1] else 0, []
    else:
        comps = []
        for i in range(len(main) - len(sub)):
            a, b, c = compare_sequences(main[i:i + len(sub)], sub)
            comps.append((i, a, b, c))
        comps = list(filter(lambda x: similarity_range[0] <= x[3] <= similarity_range[1], comps))
        return len(comps), comps


def distribute(lst):
    return list((key, lst.count(key)) for key in sorted(set(lst)))
