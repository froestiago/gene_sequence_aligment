import numpy as np
from icecream import ic

def get_key(my_dict, val):
   
    for key, value in my_dict.items():
        if val == value:
            return key
    return "key doesn't exist"

def neighbor_joining(D, taxa):
    taxa_dict = dict()
    n = len(taxa)
    for _ in range(n): taxa_dict[taxa[_]] = _

    if n == 2:
        tree = f'({taxa[0]}, {taxa[1]})'
        return (tree)


    # get min value ignoring values less than 0
    nonzero_indices = np.argwhere(D > 0)
    min_nonzero_index = np.argmin(D[nonzero_indices[:, 0], nonzero_indices[:, 1]])
    i, j = nonzero_indices[min_nonzero_index]
    print(f'pair: ({taxa[i]}, {taxa[j]})')
    branch_length = 0.5 * (D[i, j])

    # Update distance D
    new_D = np.zeros((n-1, n-1))
    new_index = 0
    for key, item in taxa_dict.items():
        if(item != i and item != j):
            new_D[new_index,0:n-2] = np.delete(D[item], [i, j])
            # use len taxa without ( , ) as multiplier and divider
            a = len(''.join(char for char in taxa[i] if char.isalpha()))
            b = len(''.join(char for char in taxa[j] if char.isalpha()))
            nova_dist = ((D[item, i])*a + (D[item, j])*b)/(a + b)
            new_D[new_index, -1] = nova_dist
            new_D[-1, new_index] = new_D[new_index, -1]
            new_index += 1
    ic(new_D)

    # concat nodes
    i, j = get_key(taxa_dict, i), get_key(taxa_dict, j)
    taxa_dict.pop(i)
    taxa_dict.pop(j)
    new_taxa = list()
    for key, value in taxa_dict.items(): new_taxa.append(key)
    
    # newick tree
    new_taxa.append(f'({i}, {j})')

    return neighbor_joining(new_D, new_taxa)


# Exemplo prof
taxa = ['A', 'B', 'C', 'D', 'E', 'F']
D = np.array([[0,   5,   4,   7,   6,   8],
              [5,   0,   7,  10,   9,  11],
              [4,   7,   0,   7,   6,   8],
              [7,  10,   7,   0,   5,   9],
              [6,   9,   6,   5,   0,   8],
              [8,  11,   8,   9,   8,   0]])

# Exemplo prof
taxa =       ['A', 'B', 'C', 'D', 'E', 'F']
D = np.array([[0,   2,   4,   6,   6,   8],
              [2,   0,   4,   6,   6,   8],
              [4,   4,   0,   6,   6,   8],
              [6,   6,   6,   0,   4,   8],
              [6,   6,   6,   4,   0,   8],
              [8,   8,   8,   8,   8,   0]])


# taxa = ['A', 'B', 'C', 'D']
# D = np.array([[0, 4, 5, 10],
#               [4, 0, 7, 12],
#               [5, 7, 0, 9],
#               [10, 12, 9, 0]])

# taxa = ['A', 'B', 'C', 'D']
# D = np.array([[0, 2, 2, 2],
#               [2, 0, 3, 2],
#               [2, 3, 0, 2],
#               [2, 2, 2, 0]])

# exemplo wikipedia
# taxa = ['A', 'B', 'C', 'D', 'E']
# D = np.array([[0, 5, 9, 9, 8],
#               [5, 0, 10, 10, 9],
#               [9, 10, 0, 8, 3],
#               [9, 10, 8, 0, 3],
#               [8, 9, 7, 3, 0]])

# exemplo https://www.youtube.com/watch?v=7tn90VWGmV4
# taxa = ['A', 'B', 'C', 'D', 'E']
# D = np.array([[0, 17, 21, 31, 23],
#               [17, 0, 30, 34, 21],
#               [21, 30, 0, 28, 39],
#               [31, 34, 28, 0, 43],
#               [23, 21, 39, 43, 0]])

ic(D)
phylogenetic_tree = neighbor_joining(D, taxa)
print(f'\n\nfinal tree: {phylogenetic_tree}')