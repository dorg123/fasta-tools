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
        if lim == -1:
            self._file.write('{}\n'.format(body))
            return
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
