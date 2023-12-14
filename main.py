from needleman import NW

def read_fasta(file_path):
    with open(file_path, 'r') as file:
        first_line = file.readline().strip()
        remaining_text = file.read().replace('\n', '')
    return first_line, remaining_text

if __name__ == "__main__":
    seq1 = read_fasta('../dataset_1/seq_1.fasta')[1]
    seq2 = read_fasta('../dataset_1/seq_2.fasta')[1]
    match = 1
    mismatch = -1
    gap = -2
    costo_extension = -2


    NW(seq1, seq2, match, mismatch, gap,costo_extension)