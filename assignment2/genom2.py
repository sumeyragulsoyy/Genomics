import random
import numpy as np
#Generate 10 DNA strings that has 
dictt={0:'A',1:'T',2:'C',3:'G'} #to generate random bases
def randomBases():
    listt=[] 
    for i in range(10): #10 string
        sublist=""
        for j in range(500):#500 bases each string
            pos=random.randint(0,3)
            sublist +=dictt[pos]
        listt.append(sublist)
    return listt

#10-mer that will be implanted to each 10 DNA strings
def Implant():
    Implant=""
    for z in range(10):
        pos1=random.randint(0,3)
        Implant+=dictt[pos1]
    return Implant
    
def mutation(implant): #get string that will be implanted with 4 mutation
    kmer=list(implant)
    index=random.randint(0,9)  #random position will mutated in these indexes
    mut_indexes={index}
    while len(mut_indexes) !=4:  #4 mutation
        index=random.randint(0,9) 
        mut_indexes.add(index)
    willMutated=list(mut_indexes)
    for i in range(4): #will be implanted 10-mer with 4 mutation
        newbase=random.randint(0,3)
        while kmer[willMutated[i]] == dictt[newbase]: #get different  base
            newbase=random.randint(0,3)
        kmer[willMutated[i]]=dictt[newbase] #mutated
    
    MutatedString="".join(kmer) 
    return MutatedString
        
#implant 10 mer each 10 Dna strings
def changeLine(line,implanted):
    #get index 0-490
    startIndex=random.randint(0,489)
    #print(startIndex)
    eachString=list(line)
    eachString[startIndex:startIndex+10]=implanted #.lower()
    return "".join(eachString)

def finalStrings(listt,Implant):#  uses mutation(),changeline() functions
    finalStrings=[]
    ImplantwithMutation=[]
    for t in range(10):
        result=mutation(Implant)
        finalStrings.append(changeLine(listt[t],"".join(result)))
        ImplantwithMutation.append(result)
    return finalStrings,ImplantwithMutation

#write to file all strings
def writeToFile(finalStrings,inputfile):
    with open(inputfile, 'w+') as f:  
        for listitem in finalStrings:
            f.write('%s\n' % listitem)
            
def readFile(inputfile):
    with open(inputfile, 'r') as f:  
        f1=f.readlines()
    return f1 #return list of strings from file


#get random 10 motifs with k length 
def motifFinding(inputfile,k): #get file and k-consensus string length
    with open(inputfile, 'r') as f:  
        motifs=[]
        f1=f.readlines()
        for x in f1:
            #get from random index with k length string
            randomStart=random.randint(0,489)
            motifs.append(x[randomStart:randomStart+k]) #get motifs from random positions  with k length
        return motifs    

#motifs=motifFinding(y,10)

def score(motifs):
    #hash each column and sum of all key values -max length =score
    #take column by column and calculate score 
    profile={'A':[],'T':[],'C':[],'G':[]} # A,T,G,C are bases and lists are k length,value in each column
    Score=0
    for i in range(len(motifs[0])):#k length
        columnby=[] #each column list
        mydict={'A':0,'T':0,'C':0,'G':0}
        for j in range(10):    #10 string
            columnby.append(motifs[j][i])
        for x in range(len(columnby)):
            mydict[columnby[x]] +=1 #update mydict for each base at column
        valueList=list(mydict.values())
        colScore=sum(valueList)-max(valueList) #
        Score +=colScore #update general score with column score
        for key in mydict: #update profile column by column
            profile[key].append(mydict[key]/10) 
        
    return Score,profile

def gibbs_score(motifs): 
    #hash each column and sum of all key values -max length =score
    #take column by column and calculate score 
    Score=0
    for i in range(len(motifs[0])):#k length
        columnby=[] #each column list
        mydict={'A':0,'T':0,'C':0,'G':0}
        for j in range(10):    #10 string
            columnby.append(motifs[j][i])
        for x in range(len(columnby)):
            mydict[columnby[x]] +=1 #update mydict for each base at column
        valueList=list(mydict.values())
        colScore=sum(valueList)-max(valueList) #
        Score +=colScore #update general score with column score
            
    return Score

#myscore,myprofile=score(motifs)

def calProb(posMotif,profile):#posMotif string,profile dict
    mult=1
    for i in range(len(posMotif)):
        mult *=profile[posMotif[i]][i]
         
    return mult  

#calProb("TTCGGGGGGG",myprofile)        

#calculate consensus string length strings in each DNA string and return max possibility as a new motif in each row
def newMotifFinding(y,profile): #get inputfile, profile as parameters and calculate new motifs
    newMotifs=[]
    k=len(profile['A']) #get consensus string length
    with open(y, 'r') as f:  
        f1=f.readlines()
        for x in f1:#lines
            prob=-1
            finalMotifInRow=""
            for y in range(0,len(x)-k):#go on the each row string
                possible=x[y:y+k]
                newprob=calProb(possible,profile)
                if newprob >prob:
                    finalMotifInRow=possible
                    prob=newprob
            newMotifs.append(finalMotifInRow)        
    return newMotifs

#newMotifFinding('input.txt',myprofile)    

def calculateBiasRandom(probabilities):# kmer-prob value
    possibleKmers=list(probabilities.keys()) #'AT..CT'
    probs=list(probabilities.values()) #probabilities of each possible motifs
    x=sum(probs)
    for i in range(len(probs)):
        probs[i]=probs[i]/x
    selected=np.random.choice(possibleKmers, 1, p=probs)
    return selected[0] #selected possible motifs


def gibbsProfile(motifs,DNAstrings): #motifs are list of strings
    probabilities={}
    randRow=random.randint(0,9) #1 row will be deleted between 0-9
    #deleted=motifs[randRow] 
    del motifs[randRow] #deleted a random row ,row is 9
    profile={'A':[],'T':[],'C':[],'G':[]} # A,T,G,C are bases and lists are k length,value in each column
    for i in range(len(motifs[0])):#k length
        columnby=[] #each column list
        mydict={'A':0,'T':0,'C':0,'G':0}
        for j in range(9):    #9 string
            columnby.append(motifs[j][i])
        for x in range(len(columnby)):
            mydict[columnby[x]] +=1 #update mydict for each base at column
        for key in mydict: #update profile with +1 and divided 18, column by column
            profile[key].append(mydict[key]+1/18) #/18
    #after profile ,calculate which motif will be selected 
    k=len(motifs[0])    
    for d in range(0,len(DNAstrings[randRow])-k):
        possible=DNAstrings[randRow][d:d+k] #traverse on deleted Dna strings
        prob=calProb(possible,profile)
        if possible not in probabilities:
            probabilities[possible]=prob
            #probabilities.append({possible:prob})
    selectedMotif=calculateBiasRandom(probabilities)
    # add the seleccted motif to current motif list and calculate score
    motifs.append(selectedMotif) #selected motif is added again motif list
    newScore=gibbs_score(motifs) #calculate new score with new motif list
    return newScore,motifs
        
#Flow begins... -------------------
#input file tek olcak implanted tek olacak
DNAstrings=randomBases() #strings
willImplanted=Implant() #implant 
afterMutated,ImplantList=finalStrings(DNAstrings,willImplanted) #strings with (10,4) implant
writeToFile(afterMutated,'inputfile.txt') #write Dna strings to file
stringsInFile=readFile('inputfile.txt')

def randomized(file,k):#file is input file,k-consensus string length
    initialMotif=motifFinding(file,k) #random motifs
    initialScore,initialProfile=score(initialMotif)
    while True:
        newMotifs=newMotifFinding(file,initialProfile)#with profile fin max prob k length in each string
        newScore,newProfile=score(newMotifs)
        if newScore < initialScore:
            initialScore=newScore
            initialProfile=newProfile
            initialMotif=newMotifs
        else:
            return initialScore,initialMotif
        

def gibbs(file,stringsInFile,k):
    initialMotif=motifFinding(file,k) #get random 10 motifs
    finalMotif=initialMotif.copy()
    finalScore=gibbs_score(initialMotif) #initial score
    while True:
        prevfinal=finalScore
        for i in range(50):
            newScore,newMotif=gibbsProfile(initialMotif,stringsInFile)
            if newScore <finalScore:
                finalScore=newScore
                finalMotif=newMotif
            else:
                initialMotif=newMotif
        if finalScore ==prevfinal:
            return finalMotif,finalScore

def consensus(motif,outputfile):
    with open(outputfile, 'w') as f:  
        for listitem in motif:
            f.write('%s\n' % listitem)
            
print(willImplanted)
print("Implant list:\n")
for i in range(len(ImplantList)):
    print(ImplantList[i])

#Randomized
Random_Score1,Random_Motif1=randomized('inputfile.txt',9) 
Random_Score2,Random_Motif2=randomized('inputfile.txt',10) 
Random_Score3,Random_Motif3=randomized('inputfile.txt',11)    
consensus(Random_Motif1,'Randomized_Output9.txt') 
consensus(Random_Motif2,'Randomized_Output10.txt')    
consensus(Random_Motif3,'Randomized_Output11.txt')       
#Gibbs
Gibbs_Motif1,Gibbs_Score1=gibbs('inputfile.txt',stringsInFile,9)
Gibbs_Motif2,Gibbs_Score2=gibbs('inputfile.txt',stringsInFile,10)
Gibbs_Motif3,Gibbs_Score3=gibbs('inputfile.txt',stringsInFile,11)
consensus(Gibbs_Motif1,'Gibbs_Output9.txt') 
consensus(Gibbs_Motif2,'Gibbs_Output10.txt')    
consensus(Gibbs_Motif3,'Gibbs_Output11.txt')       





