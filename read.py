class FileReader:
    def __init__(self, filename, accession=None):
        self._filename = filename
        self._read(accession)

    @property
    def filename(self):
        return self._filename

    @property
    def data(self):
        return self._data

    def __getitem__(self, item):
        try:
            return self._data[item]
        except KeyError:
            return None
        except IndexError:
            return None

    def __str__(self):
        return 'fasta.FileReader: {}, {} entries' \
            .format(self._filename, len(self._data))

    def __repr__(self):
        return 'fasta.FastaReader({})<length: {}>' \
            .format(self._filename, len(self._data))


class FastaReader(FileReader):
    def _read(self, accession=None, read_head=lambda x: x):
        with open(self._filename, 'r') as f:
            info = None
            self._data = dict()
            flag = False
            for line in f:
                if line.startswith('>'):
                    if flag:
                        break
                    info = read_head(line.rstrip('\n').lstrip('>'))
                    if info == accession:
                        flag = True
                else:
                    if accession is None or flag:
                        self._data[info] = self._data.get(info, '') + line.rstrip('\n')

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


class MapReader(FileReader):
    def _read(self):
        with open(self._filename, 'r') as f:
            lines = f.readlines()
        self._data = dict((a, b.replace('\t', ' - ')) for a, _, b in (l.rstrip('\n').partition('\t') for l in lines))

    def __str__(self):
        return 'fasta.SimplifiedTabReader: {}, {} entries' \
            .format(self.filename, len(self.data))

    def __repr__(self):
        return 'fasta.SimplifiedTabReader({})<length: {}>' \
            .format(self.filename, len(self.data))


class FastqReader(FileReader):
    def _read(self):
        with open(self._filename, 'r') as f:
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
                pieces = data.rstrip(',').split('$')[1:]
                pieces = list(p[:p.rfind(',')] for p in pieces[:-1]) + [pieces[-1]]
                self._data[index][title] = pieces

    def get_entry(self, accession):
        try:
            return self._data[self._accessions[accession]]
        except KeyError:
            return dict()

    def __getitem__(self, item):
        return self.get_entry(item)

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


class PdbReader(FileReader):
    def _read(self):
        with open(self._filename, 'r') as f:
            lines = list(line.strip().split() for line in f.readlines())
        d = dict(map(lambda x: (x, []), set(map(lambda x: x[0], lines))))
        for line in lines:
            d[line[0]].append(line[1:])
        self._data = d
