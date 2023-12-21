import numpy as np
from icecream import ic

def get_key(my_dict, val):
   
    for key, value in my_dict.items():
        if val == value:
            return key
        
    return "key doesn't exist"

def neighbor_joining(D, taxa):
    
    # setting dict for taxa
    taxa_dict = dict()
    for _ in range(len(taxa)): taxa_dict[taxa[_]] = _

    # recursion finish condition
    if len(taxa) == 2:
        tree = f'({taxa[0]}, {taxa[1]})'
        return (tree)

    # create Q value matrix
    n = len(taxa)
    Q = np.zeros((n, n))

    # calculate Q value for each pair
    for i in range(n):
        for j in range(i + 1, n):
            
            # formula
            Q[i, j] = D[i, j] * (n - 2) - np.sum(D[i, :]) - np.sum(D[:, j])
            
            # mirror value on matrix
            Q[j, i] = Q[i, j]
    ic(Q)
    
    # get smallest value of the matrix
    min_Q = np.argmin(Q)
    min_Q_ij = divmod(min_Q, n)

    # closest pair
    i, j = min_Q_ij
    print(f'({taxa[i]}, {taxa[j]})')

    # distance to new node (if needed)
    # delta = (np.sum(D[i, :]) - np.sum(D[:, j])) / (n - 2)
    # node_to_i = 0.5 * (D[i, j] + delta)
    # node_to_j = 0.5 * (D[i, j] - delta)

    # calculate new distances & update matrix
    new_D = np.zeros((n-1, n-1))
    new_index = 0
    for key, item in taxa_dict.items():
        if(item != i and item != j):
            new_D[new_index,0:n-2] = np.delete(D[item], [i, j])
            
            # formula
            new_distance = 0.5*(D[item, i] + D[item, j] - D[i, j])
            
            new_D[new_index, -1] = new_distance

            # mirror value on matrix
            new_D[-1, new_index] = new_D[new_index, -1]
            new_index += 1
    
    ic(new_D)

    # concat nodes
    i, j = get_key(taxa_dict, i), get_key(taxa_dict, j)
    taxa_dict.pop(i)
    taxa_dict.pop(j)
    new_taxa = list()
    for key, value in taxa_dict.items(): new_taxa.append(key)
    
    # update newick tree
    new_taxa.append(f'({i}, {j})')

    return neighbor_joining(new_D, new_taxa)

# assigment distance matrix
taxa = ['A', 'B', 'C', 'D']
original_names = {'A': 'L.braziliensis', 'B': 'T. rangeli', 'C': 'T. cruzi', 'D': 'T. gambiae'}
D = np.array([[0.000,	0.010,	0.300,	0.280],
              [0.010,	0.000,	0.280,	0.270],
              [0.300,	0.280,	0.000,	0.015],
              [0.280,	0.270,	0.015,	0.000]])
ic(D)
phylogenetic_tree = neighbor_joining(D, taxa)

ic(phylogenetic_tree)

for key, value in original_names.items():
    phylogenetic_tree = phylogenetic_tree.replace(key, value)

ic(phylogenetic_tree)
