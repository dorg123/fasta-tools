from sequence_funcs import *


def align(seq1, seq2):
    # Set up the matrix
    matrix = [[(None, None, None) for i in range(len(seq2) + 1)] for j in range(len(seq1) + 1)]
    # fill initial values
    matrix[0][0] = 0, None, None
    for i in range(1, len(seq1) + 1):
        matrix[i][0] = -2 * i, i - 1, 0
    for i in range(1, len(seq2) + 1):
        matrix[0][i] = -2 * i, 0, i - 1
    # fill the matrix
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            # calculate options
            up = matrix[i][j - 1][0] - 2, i, j - 1
            left = matrix[i - 1][j][0] - 2, i - 1, j
            across = matrix[i - 1][j - 1][0] + (1 if seq1[i - 1] == seq2[j - 1] and
                                                     seq1[i - 1] != '-' and
                                                     seq2[j - 1] != '-' else -1), i - 1, j - 1
            # and choose the best one
            matrix[i][j] = max(up, left, across, key=lambda _: _[0])
    # Set up run
    x = len(seq1)
    y = len(seq2)
    nseq1 = ''
    nseq2 = ''
    # Run through the matrix to find the sequences
    while x != 0 and y != 0:
        _, nx, ny = matrix[x][y]
        if x == nx:
            nseq1 = '-' + nseq1
            nseq2 = seq2[y - 1] + nseq2
        elif y == ny:
            nseq1 = seq1[x - 1] + nseq1
            nseq2 = '-' + nseq2
        else:
            nseq1 = seq1[x - 1] + nseq1
            nseq2 = seq2[y - 1] + nseq2
        x = nx
        y = ny
    # Return the results
    return nseq1, nseq2, matrix[len(seq1)][len(seq2)][0]


def msa(*seqs):
    pass


def align_nd(*seqs):
    # Set up variables
    n = len(seqs)
    ls = list(map(len, seqs))  # lengths

    def create_matrix(d=0):  # This recursive function creates the n-dimensional matrix.
        if d == n:
            mat = [None] * (n + 1)
        else:
            mat = create_matrix(d + 1)  # RECURSIVE CALL
        return [mat for _ in range(ls[d - 1])]

    # Create the n-dimensional matrix
    matrix = create_matrix()

    def matrix_get(mat, *j):  # Access the n-dimensional matrix (recursive)
        if len(j) == 1:
            return mat[j[0]]
        else:
            return matrix_get(mat[j[0]], *j[1:])  # RECURSIVE CALL

    def matrix_set(mat, new_val, *j):  # Set value in the n-dimensional matrix (recursive)
        if len(j) == 1:
            mat[j[0]] = new_val
        else:
            matrix_set(mat[j[0]], new_val, *j[1:])  # RECURSIVE CALL

    def increment(c):  # This function increments the counter by 1.
        for i in range(n - 1, -1, -1):
            c[i] += 1
            if c[i] == ls[i]:
                c[i] = 0
            else:
                break
        return c

    print(matrix)
    # Set up the matrix with initial values
    matrix_set(matrix, [0] + [None] * n, [0] * n)
    for i in range(n):
        index = [0] * n
        prev_index = [0] * n
        for j in range(1, ls[i]):
            index[i] = j
            val = [-2 * (n - 1) * i] + prev_index
            matrix_set(matrix, val, index)
            prev_index = index

    def score(idx1, idx2):
        pairs = list(zip(range(len(idx1)), idx1, idx2))
        sc = 0
        letters = []
        for i, a, b in pairs:
            if a == b:
                sc += -2
            else:
                letters.append(seqs[i][idx1[i]])
        letters = histogram(letters)



