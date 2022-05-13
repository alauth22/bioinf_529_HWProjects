def blosum_score(amino_acid_A, amino_acid_B):
    ''' Function to return BLOSUM62 scores for amino acid substitutsions
    
    Args:
        amino_acid_A (str): amino acid in sequence 1
        amino_acid_B (str): amino acid in sequence 2
    
    '''
    with open('data/BLOSUM62.txt') as f:
        fileContent = f.readlines()
        
    alignment_score = []
    for i in range (1, len(fileContent)):
        alignment_score = alignment_score + [ fileContent[i].split() ]
        
    alignment_score = [line[1:] for line in alignment_score]
    
    alignmentIndices = fileContent[0].split()
    
    indice_a = -1
    indice_b = -1
    for i in range(0, len(alignmentIndices)):
        if (amino_acid_A == alignmentIndices[i]):
            indice_a = i
        if (amino_acid_B == alignmentIndices[i]):
            indice_b = i
            
    return int(alignment_score[indice_a][indice_b])



