import re 
k=input('Enter k numbers: ')
x=input('Enter at least number: ')
file=input('Enter file name : ')
mylist=[]
rev={'A':'T','T':'A','G':'C','C':'G'}
mydict={}
possibleKmer=[]
with open(file,'r') as f:
        for line in f:
            for i in range(len(line)-int(k)):
                pattern=line[i:i+int(k)]
                p = re.compile(pattern)
                kmer=p.findall(line)
                mylist.append(kmer)

for g in range(len(mylist)):
    mydict[mylist[g][0]]=len(mylist[g])

for keys in mydict: #keys -> string
    z='' # revers of each k mer
    if mydict.get(keys) >= int(x):
        print(keys +' : ' +str(mydict.get(keys)))
        #get reverse comlement of possible k-mer
        for m in range(len(keys)-1,-1,-1):
            z +=rev.get(keys[m]) #z is reverse complement of each k-mer
        
        possibleKmer.append(z) # possible kmer that satisfied the conditions
        
print('Reverse complements: ')

for i in range(len(possibleKmer)):
    if possibleKmer[i] in mydict.keys():
        print(possibleKmer[i] + ' : ' + str(mydict.get(possibleKmer[i]))) 

