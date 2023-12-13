import numpy as np

class SmithWatterman():
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
                          matrix[i-1][j] + self.gap_score,
                          0]
        else:
            score_list = [matrix[i-1][j-1] + self.mismatch_score,
                          matrix[i][j-1] + self.gap_score,
                          matrix[i-1][j] + self.gap_score,
                          0]

        return np.max(score_list)

    def __fill_matrix(self):
        """
        Initialize matrix with size of seq_1 by seq_2 and assigns scores

        Returns:
            ndarray: Matrix with alignment scores
        """
        matrix = np.zeros((len(self.seq_1), len(self.seq_2)))
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

    def __get_local_align(self, coord: tuple[int, int]):
        sequence = []
        x, y = coord[0], coord[1]
        def __get_char(max_value):
            if (max_value == 0 and self.seq_1[x] == self.seq_2[y]): # match
                return "match"
            elif (max_value == 0 and self.seq_1[x] != self.seq_2[y]): # mismatch
                return "mismatch"
            elif (max_value == 1):
                return "gap_h"
            elif (max_value == 2):
                return "gap_v"

        while(x > 0 and y > 0):
            indices_to_check = [index for index in [(x - 1, y - 1), (x, y - 1), (x - 1, y)] if index is not None]
            values_to_check = [self.matrix[index] for index in indices_to_check]
            max_value = np.argmax(values_to_check)
            max_position = indices_to_check[max_value]
            # print(f"max_posistion: {max_position} - {self.matrix[x - 1][y - 1]}")
            if(self.matrix[x - 1][y - 1] == 0):
                sequence.append(__get_char(0))
                break
            else:
                sequence.append(__get_char(max_value))
                x, y = max_position
        return sequence[::-1]

    def get_align(self):
        aligns = []
        max_position = self.__find_max_positions()
        for count, coord in enumerate(max_position):
            seq_1 = self.seq_1
            seq_2 = self.seq_2
            alignment = ""
            sequence = self.__get_local_align(coord)
            # print(sequence)
            deslocamento = min(coord)
            # print(coord)
                    # add gaps
            for i, action in enumerate(sequence):
                if sequence[i] == 'gap_v':
                    self.seq_2 = self.seq_2[:i + deslocamento] + '_' + self.seq_2[i + deslocamento:]
                if sequence[i] == 'gap_h':
                    self.seq_1 = self.seq_1[:i + deslocamento] + '_' + self.seq_1[i + deslocamento:]

            seq_1 = seq_1.replace(" ", "")
            seq_2 = seq_2.replace(" ", "")
            # print(sequence)
            for _ in range(len(sequence)):
                if seq_1[_] == seq_2[_]:
                    alignment += '*'
                elif (seq_1[_] != seq_2[_]) and ((seq_1[_] != "_") and (seq_2[_] != "_")):
                    alignment += '|'
                else:
                    alignment += ' '

            print(f"{seq_1}\n{alignment}\n{seq_2}")
            file_name = 'lero' + str(count) + '.txt'
            # with open(file_name, 'w') as f: f.write(f"{seq_1}\n{alignment}\n{seq_2}")
            x, y = self.matrix.shape
            with open(file_name, 'w') as f: f.write(f"{self.seq_1}\n{alignment}\n{self.seq_2}\nscore: {self.matrix[x-1][y-1]}")

        return aligns



instance = SmithWatterman(seq_1= 'GGTTGACTA',   # 'AATCG' - 'GGTTGACTA'
                          seq_2= 'TGTTACGG')    # 'AACG' - 'TGTTACGG'

print(instance.matrix)
sequences = instance.get_align()