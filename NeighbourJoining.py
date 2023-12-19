import numpy as np

def get_key(my_dict, val):
   
    for key, value in my_dict.items():
        if val == value:
            return key
 
    return "key doesn't exist"

def neighbor_joining(D, taxa):
    
    taxa_dict = dict()
    for _ in range(len(taxa)):
        taxa_dict[taxa[_]] = _

    if len(taxa) == 2:
        tree = {'name': taxa[0], 'children': [{'name': taxa[1], 'distance': D[0, 1]}]}
        return tree

    n = len(taxa)
    Q = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            Q[i, j] = D[i, j] * (n - 2) - np.sum(D[i, :]) - np.sum(D[:, j])
            Q[j, i] = Q[i, j]

    print(f'Q:\n{Q},\n')
    min_Q = np.argmin(Q)
    min_Q_ij = divmod(min_Q, n)

    i, j = min_Q_ij
    print(f'pair: {taxa[i]}, {taxa[j]}')
    print(i, j)
    delta = (np.sum(D[i, :]) - np.sum(D[:, j])) / (n - 2)
    node_to_i = 0.5 * (D[i, j] + delta)
    print(f'node_to_i: {node_to_i}')
    node_to_j = 0.5 * (D[i, j] - delta)
    print(f'node_to_j: {node_to_j}')

    # Update distance matrix
    new_D = np.zeros((n-1, n-1))
    new_index = 0
    for key, item in taxa_dict.items():
        if(item != i and item != j):
            print(item)
            new_D[new_index,0:n-2] = np.delete(D[item], [i, j])
            # print(new_D,'\n')
            nova_dist = 0.5*(D[item, i] + D[item, j] - D[i, j])
            new_D[new_index, -1] = nova_dist
            new_D[-1, new_index] = new_D[new_index, -1]
            new_index += 1
    
    print(new_D,'\n')

    # concat nodes
    i, j = get_key(taxa_dict, i), get_key(taxa_dict, j)
    taxa_dict.pop(i)
    taxa_dict.pop(j)
    new_taxa = list()
    for key, value in taxa_dict.items(): new_taxa.append(key)
    new_taxa.append(str(i+j))


    return neighbor_joining(new_D, new_taxa)


# Exemplo prof
# taxa = ['A', 'B', 'C', 'D', 'E', 'F']
# D = np.array([[0,   5,   4,   7,   6,   8],
#               [5,   0,   7,  10,   9,  11],
#               [4,   7,   0,   7,   6,   8],
#               [7,  10,   7,   0,   5,   9],
#               [6,   9,   6,   5,   0,   8],
#               [8,  11,   8,   9,   8,   0]])

taxa = ['A', 'B', 'C', 'D']
D = np.array([[0, 4, 5, 10],
              [4, 0, 7, 12],
              [5, 7, 0, 9],
              [10, 12, 9, 0]])

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

print(f'D:\n{D}\n\n')
phylogenetic_tree = neighbor_joining(D, taxa)
print(phylogenetic_tree)


# constructor = DistanceTreeConstructor()
# NJTree = constructor.nj(D)
# print(NJTree)