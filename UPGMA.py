import numpy as np

matrix = np.asanyarray(
         [[0, 9, 9, 9, 9],
          [2, 0, 9, 9, 9],
          [4, 4, 0, 9, 9],
          [6, 6, 6, 0, 9],
          [6, 6, 6, 4, 0],
          [8, 8, 8, 8, 8]])

def find_min_distance_positions(matrix):
    """
    Find the indices of all occurrences of the maximum value
    Returns:
        list: All matrix indices with the maximum value
    """
    
    flat_indices = np.flatnonzero(matrix == np.min(matrix[matrix>0]))

    shape = matrix.shape
    positions = [(index // shape[1], index % shape[1]) for index in flat_indices]
    return positions[0], positions[1]

# def clean_matrix(matrix):
#     i, j = matrix.shape
#     for i in range(matrix.shape[1]):
#         matrix[]

result = find_min_distance_positions(matrix)
print(result)
