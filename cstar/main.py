MatchScore = 1
MisMatchScore = -1
GapScore = -4

from itertools import combinations
from multiprocessing import Pool
class Cell:
    def __init__(self, score):
        self.score = score
        self.tracebacks_array = []
def initialize_cells(row, col):
    return [[Cell(0) for i in range(row)]
            for j in range(col)]
class Algorithm:
    def __init__(self, s1, s2, scores):
        self.sequence1 = s1
        self.sequence2 = s2
        self.match, self.mismatch, self.gap = list(map(int, scores))
        self.cells = initialize_cells(len(s2) + 1, len(s1) + 1)
    def nw(self, onlyOne=False):
        self.prepare()

        for i in range(1, len(self.sequence1) + 1):
            for j in range(1, len(self.sequence2) + 1):
                nP = self.cells[i - 1][j].score + self.gap
                wP = self.cells[i][j - 1].score + self.gap
                nwP = self.cells[i - 1][j - 1].score
                if self.sequence1[i - 1] == self.sequence2[j - 1]:
                    nwP += self.match
                else:
                    nwP += self.mismatch
                max_ = max(nwP, wP, nP)
                self.cells[i][j].score = max_
                if max_ == nP:
                    self.cells[i][j].tracebacks_array.append('n')
                if max_ == wP:
                    self.cells[i][j].tracebacks_array.append('w')
                if max_ == nwP:
                    self.cells[i][j].tracebacks_array.append('nw')
        return {'nw': self.tracebacks(onlyOne), 'score': self.cells[-1][-1].score}
    def tracebacks(self, onlyOne):
        stack = [(len(self.sequence1), len(self.sequence2), '', ''), ]
        results = []
        while True:
            i, j, s1, s2 = stack.pop()
            while i > 0 or j > 0:
                if len(self.cells[i][j].tracebacks_array) > 1:
                    stack.append((i, j, s1, s2))
                    pole = self.cells[i][j].tracebacks_array.pop()
                else:
                    pole = self.cells[i][j].tracebacks_array[0]
                if pole == 'n':
                    s2 += '-'
                    i -= 1
                    s1 += self.sequence1[i]
                elif pole == 'nw':
                    j -= 1
                    s2 += self.sequence2[j]
                    i -= 1
                    s1 += self.sequence1[i]
                else:
                    j -= 1
                    s2 += self.sequence2[j]
                    s1 += '-'
            results.append((s1[::-1], s2[::-1]))
            if onlyOne or not stack:
                return results

    def prepare(self):
        for i in range(len(self.cells)):
            self.cells[i][0].score = self.gap * i
            self.cells[i][0].tracebacks_array.append('n')

        for i in range(len(self.cells[0])):
            self.cells[0][i].score = self.gap * i
            self.cells[0][i].tracebacks_array.append('w')


def align_similar(s1, s2):
    change1, change2 = list(), list()
    i = 0
    while s1 != s2:
        if i > len(s1) - 1:
            s1 += s2[i:]
            change1.extend(range(i, i + len(s2[i:])))
            continue
        if i > len(s2) - 1:
            s2 += s1[i:]
            change2.extend(range(i, i + len(s1[i:])))
            continue
        if s1[i] != s2[i]:
            if s1[i] == '-':
                s2 = s2[0:i] + '-' + s2[i:]
                change2.append(i)
            else:
                s1 = s1[0:i] + '-' + s1[i:]
                change1.append(i)
        i += 1
    return sorted(change1), sorted(change2)
def adjust(string_list, indices):
    for i, string in enumerate(string_list):
        for index in indices:
            string = string[:index] + '-' + string[index:]
        string_list[i] = string
def worker(it):
    ((i, string_i), (j, string_j)), scores = it
    model = Algorithm(string_i, string_j, scores).nw(True)
    (string_ai, string_aj), score = model['nw'][0], model['score']
    return (i, string_ai), (j, string_aj), score
class CenterStar:

    def __init__(self, scores, strings):
        self.score_rate = scores
        self.sequences = strings
        self.matrix = [[0] * (len(strings) + 1) for _ in range(len(strings))]

    def multiple_sequence_alignment(self):
        msa_result = []
        max_row, max_value = 0, 0
        len_strings = len(self.sequences)
        dos = tuple(combinations(zip(range(len_strings), self.sequences), 2))
        dos = zip(dos, (self.score_rate for _ in range(len(dos))))
        pool = Pool()
        result = pool.map(worker, dos)
        for elem in result:
            (i, string_i), (j, string_j), score = elem
            self.matrix[i][j] = (string_i, string_j, score)
            self.matrix[j][i] = (string_j, string_i, score)
            self.matrix[i][-1] += score
            self.matrix[j][-1] += score

            if self.matrix[j][-1] > max_value:
                max_row = j
                max_value = self.matrix[j][-1]
            if self.matrix[i][-1] > max_value:
                max_row = i
                max_value = self.matrix[i][-1]

        for i in range(len_strings):
            if i == max_row:
                continue
            if not msa_result:
                msa_result.extend(self.matrix[max_row][i][0: 2])
                continue

            new = list(self.matrix[max_row][i][0: 2])
            ch_index1, ch_index2 = align_similar(msa_result[0], new[0])

            adjust(msa_result, ch_index1)
            adjust(new, ch_index2)
            msa_result.extend(new[1:])

        return msa_result
lines = []
for i in range(4):
    lines.append(input())
scores = [MatchScore, MisMatchScore, GapScore]
msa = CenterStar(scores, lines).multiple_sequence_alignment()
for item in msa:
    print(item)
