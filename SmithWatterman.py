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
    
    # # # # # AQUI # # # # # 
    # assign char (*, | or _) based on moviment
    def __get_local_align(self, coord: tuple[int, int]):
        sequence = ""
        current_value, next_value, max_value = None, None, 0
        x, y = coord[0], coord[1]
        def __get_char(max_value):
            print(self.seq_1[x], ' - ', self.seq_2[y])
            if (max_value == 0 and self.seq_1[x] == self.seq_2[y]): # match
                return "*"
            elif (max_value == 0 and self.seq_1[x] != self.seq_2[y]): # mismatch
                return "|"
            elif (max_value == 1 or max_value == 2): # gap
                return "_"

        while(self.matrix[x][y]):
            if(self.matrix[x - 1][y - 1] == 0):
                sequence += __get_char(0)
                return sequence[::-1]
            else:
                indices_to_check = [index for index in [(x - 1, y - 1), (x, y - 1), (x - 1, y)] if index is not None]
                values_to_check = [self.matrix[index] for index in indices_to_check]
                max_value = np.argmax(values_to_check)
                max_position = indices_to_check[max_value]
                sequence += __get_char(max_value)
                x, y = max_position
        return sequence[::-1]

    def get_align(self):
        aligns = []
        max_position = self.__find_max_positions()
        for coord in max_position:
            aligns.append(self.__get_local_align(coord))

        return aligns



instance = SmithWatterman(seq_1= 'AATCG',   # 'AATCG' - 'GGTTGACTA'
                          seq_2= 'AACG')    # 'AACG' - 'TGTTACGG'

print(instance.matrix)
sequences = instance.get_align()
for alignment in sequences:
    print(f"{instance.seq_1[:len(alignment)+1]}\n {alignment}\n{instance.seq_2[:len(alignment)+1]}")