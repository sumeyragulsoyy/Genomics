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
                pattern=line[i:i+int(k)] #get pattern
                p = re.compile(pattern) 
                kmer=p.findall(line) #find all existence in line ,return list
                mylist.append(kmer) #add to list each kmer 

for g in range(len(mylist)):
    mydict[mylist[g][0]]=len(mylist[g]) #'ACA' -> 3 ,kmer and count in line as a dictionary
#Print all possible kmer that satisfied the given condition ,at least x times
flag1=0 #to print kmer that satisfied the given condition
flag2=0 #to print reverse complement of possible k mers
for keys in mydict: #keys -> string
    z='' # revers of each k mer
    if mydict.get(keys) >= int(x):
        flag1=1
        print(keys +' : ' +str(mydict.get(keys)))
        #get reverse comlement of possible k-mer
        for m in range(len(keys)-1,-1,-1):
            z +=rev.get(keys[m]) #z is reverse complement of each k-mer
        
        possibleKmer.append(z) # possible kmer that satisfied the conditions

if flag1==0:
    print("There is no %s mer, at least %s times " % (k, x))        
print('Reverse complements: ')

for i in range(len(possibleKmer)):
    if possibleKmer[i] in mydict.keys():
        flag2=1
        print(possibleKmer[i] + ' appearing ' + str(mydict.get(possibleKmer[i])) +' times') 
 
if flag2==0:
    print("There is no reverse of possible %s mer"%(k))