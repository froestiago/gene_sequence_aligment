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
        return np.zeros((len(self.seq_1), len(self.seq_2)))
    
    def __match_score(self, matrix, i, j):
        return matrix[i-1][j-1] + self.match_score
    
    def __mismatch_score(self, matrix, i, j):
        return matrix[i-1][j-1] + self.mismatch_score

    def __gap_horizontal_score(self, matrix, i, j):
        return matrix[i][j-1] + self.gap_score
    
    def __gap_vertical_score(self, matrix, i, j):
        return matrix[i-1][j] + self.gap_score

    def get_score(self, matrix, i, j):
        print(f"comparing {self.seq_1[i]} - {self.seq_2[j]}")
        if(self.seq_1[i] == self.seq_2[j]):
            print(f"\tigual")
            score_list = [self.__match_score(matrix, i, j),
                          self.__gap_horizontal_score(matrix, i, j),
                          self.__gap_vertical_score(matrix, i, j),
                          0]
        else:
            print(f"\tdiferente")
            score_list = [self.__mismatch_score(matrix, i, j),
                          self.__gap_horizontal_score(matrix, i, j),
                          self.__gap_vertical_score(matrix, i, j),
                          0]

            
        print(f"\tscore_list: {score_list}")
        print(f"\tscore: {np.max(score_list)}\n")
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
        # Find the indices of all occurrences of the maximum value in the flattened matrix
        flat_indices = np.flatnonzero(self.matrix == np.max(self.matrix))
    
        # Convert the flattened indices to 2D indices (row, column)
        shape = self.matrix.shape
        positions = [(index // shape[1], index % shape[1]) for index in flat_indices]
    
        return positions

    # def __find_next_element(self, coord):
    #     line = coord[0]
    #     column = coord[1]
    #     # sequence = self.matrix[line][column]
    #     up_left = [(line - 1, column - 1), self.matrix[(line - 1,column - 1)]]
    #     left = [(line, column -1), self.matrix[(line, column -1)]]
    #     up = [(line -1, column), self.matrix[(line - 1, column)]]
    #     positions_values = [up_left, left, up]

    #     max_value = np.max(positions_values)
    #     # max_value = np.max([up_left[1], left[1], up[1]])
    #     print(f"max_value: {max_value}")


    #     print(f"up_left -> coord: {up_left[0]} | value: {up_left[1]}")
    #     print(f"left -> coord: {left[0]} | value: {left[1]}")
    #     print(f"up -> coord: {up[0]} | value: {up[1]}")
    #     # print(f"{up_left} - {left} - {up}")

    def find_max_neighbor(self, coord):
        
        next_value = None
        sequence = ""
        while(next_value != 0):
            # Get the indices for left, up, and diagonal-up-left
            x, y = coord[0], coord[1]
            left_index = (x, y - 1) if y > 0 else None
            up_index = (x - 1, y) if x > 0 else None
            diagonal_up_left_index = (x - 1, y - 1) if x > 0 and y > 0 else None

            # Collect the indices in a list
            indices_to_check = [index for index in [left_index, up_index, diagonal_up_left_index] if index is not None]

            # Get values at the specified indices and calculate the maximum
            values_to_check = [self.matrix[index] for index in indices_to_check]
            max_value = np.max(values_to_check)

            # Find the first occurrence of the maximum value in the list
            max_position = indices_to_check[values_to_check.index(max_value)]
            next_value = 0

        return sequence
        

    def get_align(self):
        alings = []
        max_position = self.__find_max_positions()
        print(max_position)
        for coord in max_position:
            print(self.find_max_neighbor(coord))



lero = SmithWatterman(seq_1= 'WHAT',
                      seq_2= 'WHY')

print(lero.matrix)
max = lero.get_align()
# print(max)