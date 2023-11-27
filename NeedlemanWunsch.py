import numpy as np

class NeedlemanWunsch():
    def __init__(self,
                 seq_1:str,
                 seq_2:str,
                 match_score:int = 1,
                 mismatch_score:int = -1,
                 gap_score:int = -2):
        self.seq_1 = " " + seq_1
        self.seq_2 = " " + seq_2
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_score = gap_score
        self.matrix = self.__fill_matrix()

    def __get_score(self, matrix, i, j):
        """
        Check if bases are the same or not and returns the score

        Args:
            matrix (ndarray):
            i (int):
            j (int):

        Returns:
            int: Max score
        """
        if(self.seq_1[i] == self.seq_2[j]):
            score_list = [matrix[i-1][j-1] + self.match_score,
                          matrix[i][j-1] + self.gap_score,
                          matrix[i-1][j] + self.gap_score]
        else:
            score_list = [matrix[i-1][j-1] + self.mismatch_score,
                          matrix[i][j-1] + self.gap_score,
                          matrix[i-1][j] + self.gap_score]

        return np.max(score_list)

    def __fill_matrix(self):
        """
        Initialize matrix with size of seq_1 by seq_2 and assigns scores

        Returns:
            ndarray: Matrix with alignment scores
        """
        matrix = np.zeros((len(self.seq_1), len(self.seq_2)))
        for _ in range(matrix.shape[0]): matrix[_][0] = -2*_
        for _ in range(matrix.shape[1]): matrix[0][_] = -2*_

        for i in range(matrix.shape[0] - 1):
            i += 1
            for j in range (matrix.shape[1] - 1):
                j += 1
                x = self.__get_score(matrix, i, j)
                matrix[i][j] = x
        return matrix
    
    def __find_max_positions(self):
        """
        Find the indices of all occurrences of the maximum value

        Returns:
            list: All matrix indices with the maximum value
        """
        
        flat_indices = np.flatnonzero(self.matrix == np.max(self.matrix))
    
        shape = self.matrix.shape
        positions = [(index // shape[1], index % shape[1]) for index in flat_indices]
        return positions
    
    # # # # # AQUI # # # # # 
    # assign char (*, | or _) based on moviment
    def __get_global_align(self, coord: tuple[int, int]):
        sequence = []
        x, y = coord[0], coord[1]
        def __get_char(max_value):
            print(self.seq_1[x], ' - ', self.seq_2[y])
            if (max_value == 0 and self.seq_1[x] == self.seq_2[y]): # match
                return "match"
            elif (max_value == 0 and self.seq_1[x] != self.seq_2[y]): # mismatch
                return "mismatch"
            elif (max_value == 1):
                return "gap_h"
            elif (max_value == 2):
                return "gap_v"

        while(self.matrix[x][y]):
            indices_to_check = [index for index in [(x - 1, y - 1), (x, y - 1), (x - 1, y)] if index is not None]
            values_to_check = [self.matrix[index] for index in indices_to_check]
            max_value = np.argmax(values_to_check)
            max_position = indices_to_check[max_value]
            sequence.append(__get_char(max_value))
            x, y = max_position
        return sequence[::-1]

    def get_align(self):
        alignment = ''
        max_position = self.__find_max_positions()
        sequence = self.__get_global_align(max_position[-1])
        for _ in range(len(sequence)):
            if sequence[_] == 'match':
                alignment += '*'
            elif sequence[_] == 'mismatch':
                alignment += '|'
            elif sequence[_] == 'gap_v':
                alignment += ' '
                self.seq_2 = self.seq_2[:_] + '_' + self.seq_2[_:]
            elif sequence[_] == 'gap_h':
                alignment += ' '
                self.seq_1 = self.seq_1[:_] + '_' + self.seq_1[_:]

        print(f"{self.seq_1[1:]}\n{alignment}\n{self.seq_2[1:]}")
        return alignment



instance = NeedlemanWunsch(seq_1= 'AATCG',   # 'AATCG' - 'GGTTGACTA'
                          seq_2= 'AACG')    # 'AACG' - 'TGTTACGG'

print(instance.matrix)
sequences = instance.get_align()
# print(sequences)