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

def kmer(query,w=11):
    mapu = {}
    for i in range(0,len(query)-w+1):
        sub_seq = query[i:i+w]
        if sub_seq in mapu:
            mapu[sub_seq].append(i)
        else:
            mapu[sub_seq] = [i]
    return mapu

def get_hits(seq,m,w=11):
    res = []
    for i in range(0,len(seq)-w+1):
        sub_seq = seq[i:i+w]
        if sub_seq in m:
            l = m[sub_seq]
            for ind in l:
                res.append((ind,i))
    return res
    
def extend_hits(seq,hit,query,w=11)
    stq,sts = hit[0],hit[1]
    matfw = 0
    k = 0
    bestk = 0
    while 2*matfw >= k and stq + w + k < len(query) and sts + w + k < len(seq):
        if query[sts+w+k] == seq[sts+w+k]:
            matfw += 1
            bestk = k + 1
        k += 1
    size = w + bestk
    k = 0
    matbw = 0
    bestk = 0
    while 2*matbw >= k and stq > k and sts > k:
        if query[stq-k-1] == seq[sts-k-1]:
            matbw +=1 
            bestk = k+1
        k+=1
    size += bestk
    return (stq-bestk,sts-bestk,size,w+matfw+matbw)
    
def hit_best_score(seq,query,m,w):
    hits = get_hits(seq,m,w)
    best_score = -1.0
    for h in hits:
        ext = extends_hit(seq,h,query,w)
        score = ext[3]
        if score > best_score or
seq = fasta_reader('sequence.fasta')[0]
query = fasta_reader('sequence.fasta')[1]
m = kmer(query,11)
get_hit = get_hits(seq,m,11)
print(get_hit)
