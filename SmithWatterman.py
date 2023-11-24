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

    def __init_matrix(self):
        """
        Initilize matrix with the size seq_1 by seq_2 with zeros

        Returns:
            ndarray: Matrix with size seq_1 by seq_2 with zeros
        """
        return np.zeros((len(self.seq_1), len(self.seq_2)))

    def get_score(self, matrix, i, j):
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
        matrix = self.__init_matrix()
        for i in range(matrix.shape[0] - 1):
            i += 1
            for j in range (matrix.shape[1] - 1):
                j += 1
                x = self.get_score(matrix, i, j)
                matrix[i][j] = x
        return matrix
    
    def __find_max_positions(self):
        """
        Find the indices of all occurrences of the maximum value

        Returns:
            list: All matrix indices with the maximum value
        """
        flat_indices = np.flatnonzero(self.matrix == np.max(self.matrix))
    
        # Convert the flattened indices to 2D indices (row, column)
        shape = self.matrix.shape
        positions = [(index // shape[1], index % shape[1]) for index in flat_indices]
    
        return positions
    
    # # # # # AQUI # # # # # 
    def __get_local_align(self, coord):
        sequence = ""
        current_value, next_value = None, None
        x, y = coord[0], coord[1]
        while(self.matrix[x][y] != 0):
            sequence = sequence + self.seq_1[x]
            x -= 1
            y -= 1
        return sequence[::-1]

    def get_align(self):
        aligns = []
        max_position = self.__find_max_positions()
        # print(f"max_position: {max_position}")
        for coord in max_position:
            # print(f"coord - {coord}")
            aligns.append(self.__get_local_align(coord))

        return aligns



instance = SmithWatterman(seq_1= 'AUGUT', 
                          seq_2= 'AUCUT')

print(instance.matrix)
sequences = instance.get_align()
print(f"sequences: {sequences}")