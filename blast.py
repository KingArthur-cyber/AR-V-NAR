import time
st = time.time()
import numpy as np
def reverse_complement(seq):
    rev_comp_seq = ''
    seq = seq[::-1]
    for i in range(0,len(seq)):
        if seq[i] == 'A':
            rev_comp_seq += 'T'
        elif seq[i] == 'T':
            rev_comp_seq += 'A'
        elif seq[i] == 'G':
            rev_comp_seq += 'C'
        else:
            rev_comp_seq += 'G'
    return rev_comp_seq
    
        
def fasta_reader(file_name):
    with open(file_name) as f:
        new_dict = {}
        
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                
                key = line.lstrip('>')
                new_dict[key] = []
            else:
                
                new_dict[key].append(line)
    for key,value in new_dict.items():
        value = ''.join(value)
        new_dict[key] = value
    value = []
    for values in new_dict.values():
        value.append(values)
    return value
    
        
def edit_dist(query,seq):
    matrix = np.zeros((len(query)+1,len(seq)+1))
    for i in range(0,len(seq)+1):
        matrix[i][0] = i 
        matrix[0][i] = i 
    for i in range(1,len(seq)+1):
        for j in range(1,len(query)+1):
            left = matrix[i][j-1]
            right =  matrix[i-1][j]
            diag = matrix[i-1][j-1]
            minim  =  min(left,right,diag)
            if seq[i-1] == query[j-1]:
                matrix[i][j] = minim
            elif seq[i-1] != query[j-1]:
                matrix[i][j] = minim + 1
    return matrix[len(seq)][len(query)]

def kmer(query,w=11):
    hash_map = {}
    for i in range(0,len(query)-w+1):
        sub_seq = query[i:i+w]
        if sub_seq in hash_map:
            hash_map[sub_seq].append(i)
        else:
            hash_map[sub_seq] = [i]
    return hash_map
def hits(seq,query_map,w=11):
    k_index = []
    for i in range(0,len(seq)- w + 1):
        seq_mer = seq[i:i+w]
        if seq_mer in query_map:
            l = query_map[seq_mer]
            #print(l)
            for ind in l:
                k_index.append((ind,i))
    return k_index
    
def hamming_dist(seq,query):
    count = 0
    length = len(query)
    for i in range(0,len(seq)):
        if seq[i] == query[i]:
            pass
        elif seq[i] != query[i]:
            count += 1
    percent = (count/length) * 100
    return percent
    
def extend_hits(query,seq,indexes,w=11):
    index_q = []
    index_s = []
    answers = []
    for i in range(0,len(indexes)):
        index_q.append(indexes[i][0])
        index_s.append(indexes[i][1])
    
    for i in range(0,len(index_q)):
        index_f_q = index_q[i]
        index_b_q = index_q[i]
        
        index_f_s = index_s[i]
        index_b_s = index_s[i]
        score = 0
        count = 0
        threshold_score = 95
        #print(i)
        if index_q[i] >-1:
            while (index_f_q < (len(query)-1) and index_f_s < (len(seq)-1)) and (score < threshold_score):
                score = hamming_dist(seq[index_b_s:index_f_s+1],query[index_b_q:index_f_q+1])
                
                index_f_q += 1
                index_f_s += 1
                if index_b_q > 0 and index_b_s > 0:
                    index_b_q -= 1
                    index_b_s -= 1
                elif index_b_q <= 0:
                    index_b_q = 0
                if index_b_s <= 0:
                    index_b_s = 0
                
                count += 1
                #print(count)
        if (score < 95) and (abs(index_q[i] - index_f_q) > 600):
            answers.append((index_q[i],index_s[i],index_f_q,index_f_s))
    return answers
            
       
            
            
        
        
            
            
db_seq = fasta_reader('sequence.fasta')[0]
rev_seq = reverse_complement(db_seq)
query = fasta_reader('fusb.fasta')[0]


hash_map = kmer(query)
#print(db_seq)

#print(len(hits(db_seq,hash_map)))
hitu = hits(rev_seq,hash_map)
print(len(rev_seq))
print(hitu)
print(extend_hits(query,rev_seq,hitu))
et = time.time()
elapsed_time = et-st
print("Execution time: ",elapsed_time,"seconds")
