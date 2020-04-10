from itertools import combinations
from multiprocessing import Pool
class Cell:
    def __init__(self, score):
        self.score = score
        self.tracebacks_array = []


def initial_cells(n, m):
    return [[Cell(0) for _ in range(n)] for _ in range(m)]


class NeedlemanWunsch:

    def __init__(self, string1, string2, scores):
        self.string1 = string1
        self.string2 = string2
        self.match, self.mismatch, self.gap = list(map(int, scores))
        self.cells = initial_cells(len(string2) + 1, len(string1) + 1)

    def nw(self, onlyOne=False):
        self.__prepare_matrix()

        for i in range(1, len(self.string1) + 1):
            for j in range(1, len(self.string2) + 1):
                nP = self.cells[i - 1][j].score + self.gap
                wP = self.cells[i][j - 1].score + self.gap
                nwP = self.cells[i - 1][j - 1].score

                if self.string1[i - 1] == self.string2[j - 1]:
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

        return {'nw': self.tracebacks(onlyOne),
                'score': self.cells[-1][-1].score}

    def tracebacks(self, onlyOne):
        
        stack = [(len(self.string1), len(self.string2), '', ''), ]
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
                    s1 += self.string1[i]
                elif pole == 'nw':
                    j -= 1
                    s2 += self.string2[j]
                    i -= 1
                    s1 += self.string1[i]
                else:
                    j -= 1
                    s2 += self.string2[j]
                    s1 += '-'
            results.append((s1[::-1], s2[::-1]))
            if onlyOne or not stack:
                return results

    def __prepare_matrix(self):
        for i, row in enumerate(self.cells):
            row[0].score = self.gap * i
            row[0].tracebacks_array.append('n')

        for i, cell in enumerate(self.cells[0]):
            cell.score = self.gap * i
            cell.tracebacks_array.append('w')


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
    model = NeedlemanWunsch(string_i, string_j, scores).nw(True)
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

        tasks = tuple(combinations(zip(range(len_strings), self.sequences), 2))
        tasks = zip(tasks, (self.score_rate for _ in range(len(tasks))))

        with Pool() as pool:
            result = pool.map(worker, tasks)
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


input_file = "input.txt"
output_file = "out.txt"
with open(input_file, 'r') as f:
    lines = [line.strip() for line in f.readlines() if line.strip()]
    scores = lines.pop(0).split('|')
    msa = CenterStar(scores, lines).multiple_sequence_alignment()

    with open(output_file, 'w') as out:
        out.writelines(map(lambda x: x + '\n', msa))
