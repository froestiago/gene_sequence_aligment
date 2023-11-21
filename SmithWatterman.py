import numpy as np

class SmithWatterman():
    def __init__(self,
                 seq_1,
                 seq_2,
                 match_score,
                 mismatch_score,
                 gap_score):
        self.seq_1 = seq_1
        self.seq_2 = seq_2
        self.matrix = self.__fill_matrix()
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_score = gap_score

    def __init_matrix(self):
        return np.zeros((len(self.seq_1)+1, len(self.seq_2)+1))

    def get_score(self, matrix, i, j):
        match_score = matrix[i-1][j-1] + self.match_score
        mismatch_score = matrix[i-1][j-1] - self.mismatch_score
        gap_horizontal_score = matrix[i][j-1] - self.gap_score
        gap_vertical_score = matrix[i-1][j] - self.gap_score
        return np.max([match_score,
                       mismatch_score,
                       gap_horizontal_score,
                       gap_vertical_score])

    def __fill_matrix(self):
        matrix = self.__init_matrix()
        for i in range(len(matrix[0][:]) - 1):
            for j in range (len(matrix[:][0]) - 1):
                x = self.get_score(matrix, i, j)
                print(f"WTF bernardo?")
                # print(i, j)
        return matrix
    



lero = SmithWatterman(seq_1= 'UAU',
                      seq_2= 'GAG',
                      match_score= 1,
                      mismatch_score= -1,
                      gap_score= -2)

print(lero.matrix)
# lero.fill_matrix()
# print(lero.matrix)