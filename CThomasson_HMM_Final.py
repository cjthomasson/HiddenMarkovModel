
'''
Caley Thomasson 

Problem Set #3
Hidden Markov Model
Viterbi Algorithm

'''
import numpy as np
import argparse


def getFrequencyT(seqFile): 
    '''
    gets frequency of amino acids in soluble sequence file in dictionary: solFrequencyDictionary
    '''
    seqFile = open(seqFile, "r")
    seqFile = seqFile.read()
    seqFile = seqFile.rstrip()
    total_length = 0 
    aaDict = {}                                     # aaDictionary format: {amino acid: CountofAminoAcid}
     
    for aa in seqFile:
        if aa != "\n" and aa != " ":
            if aa in aaDict:                          # if aaDictionary.has_key(aa)
                aaDict[aa] += 1                        #then aaDictionary[aa] += 1                                                    
            else:
                aaDict[aa] = 1                        # if aaDictionary[does not have amino acid from aa in for loop listed]
                                                        #add it so count of aa in aaDictionary now = 1 
            total_length += 1                               

    probadd = 0   
    for n in aaDict.keys():                         #for aakey in dictionary
        fraction = aaDict[n] / (total_length)          # do this math: to get %age of amino acid out of whole transmembrane sequence file     
        aaDict[n] = aaDict[n]/(total_length)
        probadd += fraction
    return aaDict

def getFrequencyS(seqFile): 
    '''
    gets frequency of amino acids in soluble sequence file in dictionary: solFrequencyDictionary
    '''
    seqFile = open(seqFile, "r")
    seqFile = seqFile.read()
    seqFile = seqFile.rstrip()
 
    total_length = 0 
    aaDict = {}                                       # aaDictionary format: {amino acid: CountofAminoAcid}
     
    for aa in seqFile:
        if aa != "\n" and aa != " ":
            if aa in aaDict:                          # if aaDictionary.has_key(aa)
                aaDict[aa] += 1                        #then aaDictionary[aa] += 1  
            else:
                aaDict[aa] = 1                        # if aaDictionary[does not have aa]
                                                        #then add it so count of aa in aaDictionary now = 1 
            total_length += 1                               
    probadd = 0   

    for n in aaDict.keys():                         #for aakey in dictionary
       
        fraction = aaDict[n] / (total_length)          # do this math: to get %age of aa out of entire soluble sequence file
        aaDict[n] = aaDict[n]/(total_length)
        probadd += fraction

    return aaDict

def getCatSeqForMatrix(testFile):
    testFile = open(testFile, "r")
    catSeq = ""
    for line in testFile: 
        catSeq += line.strip("\n")
    testFile.close()
    return(catSeq)

def getFirstScore(testFile, p_startsWithS, p_startsWithT): 
    '''
        called in stateProb() to get beginning probabilities
    '''
    testFile = open(testFile, "r")  
    
    for line in testFile: 
        first_aa = line[0]        
    solStartProbForFirstAA = np.log(p_startsWithS)+ np.log(solFrequencyDictionary.get(first_aa))            
    transStartProbForFirstAA = np.log(p_startsWithT) + np.log(transFrequencyDictionary.get(first_aa)) 

    testFile.close()
    return(solStartProbForFirstAA, transStartProbForFirstAA)

def stateProb(testFile):    
    '''
    Calculates scoring matrix and stores max score for each state and the state it was calculated from for trace back
    
    '''
    solStartProbForFirstAA, transStartProbForFirstAA= getFirstScore(testFile, p_startsWithS, p_startsWithT) #start probabilities
    probOfaaInSolState = 0
    probOfaaInTransState = 0    
    aaSequence = getCatSeqForMatrix(testFile)  
    testFile = open(testFile, "r") # open after getFirstScore() is called, it opens and closes file as well. if left open or opened before, file will be at end of file and loop wont work
    rows = 2 # number of states
    columns = len(aaSequence) #number of amino acids in sequence
    aaMatrix = np.zeros((rows+2, columns+1), dtype = np.object)

    for colj in range(1, columns+1): 
        aaMatrix[0][colj] = aaSequence[colj-1] # sets column headers from amino acids in seq
 
    aaMatrix[2][0] = "S" #row label for soluble state
    aaMatrix[3][0] = "T" #row label for transmembrane state 
  
    previousTransScore = transStartProbForFirstAA
    previousSolScore = solStartProbForFirstAA
    aaMatrix[2][1] = solStartProbForFirstAA #S
    aaMatrix[3][1] = transStartProbForFirstAA #T 
    
    for colIndex in range(1, len(aaSequence)): 
        probOfaaInSolState = solFrequencyDictionary.get(aaSequence[colIndex]) 
        probOfaaInTransState = transFrequencyDictionary.get(aaSequence[colIndex])
   
        if((previousSolScore + (np.log(probOfaaInSolState)) + np.log(p_ss)) >  (previousTransScore + np.log(probOfaaInSolState) + np.log(p_ts))):
            maxSolScore = (previousSolScore + (np.log(probOfaaInSolState)) + np.log(p_ss))
            transitionDirectionSol = "S"  #direction the max score comes from to use for trace back

        if((previousSolScore + (np.log(probOfaaInSolState)) + np.log(p_ss)) < (previousTransScore + np.log(probOfaaInSolState) + np.log(p_ts))):
            maxSolScore = previousTransScore + np.log(probOfaaInSolState) + np.log(p_ts)
            transitionDirectionSol = "T"
     
        if((previousSolScore + (np.log(probOfaaInTransState)) + np.log(p_st)) >  (previousTransScore + np.log(probOfaaInTransState) + np.log(p_tt))):
            maxTransScore = previousSolScore + (np.log(probOfaaInTransState)) + np.log(p_st)
            transitionDirectionTrans = "S" 
      
        if((previousSolScore + (np.log(probOfaaInTransState)) + np.log(p_st)) <  (previousTransScore + np.log(probOfaaInTransState) + np.log(p_tt))):
            maxTransScore = previousTransScore + np.log(probOfaaInTransState) + np.log(p_tt)
            transitionDirectionTrans = "T"
          
        previousTransScore = maxTransScore
        previousSolScore = maxSolScore
        aaMatrix[3][colIndex+1] = maxTransScore, transitionDirectionTrans    
        aaMatrix[2][colIndex+1] = maxSolScore, transitionDirectionSol
    
    #begin trace back   
    traceBack = ""
    solEndScore = aaMatrix[2][colIndex+1]
    transEndScore = aaMatrix[3][colIndex+1]   
    stateIndex = 3
    if(solEndScore[0]> transEndScore[0]):
        stateIndex = 2
        traceBack += "S"
    if(transEndScore[0]> solEndScore[0]):
        stateIndex = 3
        traceBack += "T"   
   
    for k in range(colIndex+1, 1, -1):
        traceBackTuple = aaMatrix[stateIndex][k]
        
        if(traceBackTuple[1] == "S"):  
            stateIndex = 2
            traceBack += "S"
 
        if(traceBackTuple[1] == "T"):
            stateIndex = 3
            traceBack += "T"
      
    return (reverse(traceBack))
    testFile.close()
    
def reverse(regularString):
    stringToReverse = regularString
    reverseString=""
    for i in range(len(stringToReverse)-1,-1,-1):
        reverseString += stringToReverse[i]
        
    return(reverseString)      
  
def getTransProb(stateSeqFile, k): 
    '''
Calculates the probability of transition states: SS, ST, TS, TT 
Calculates the beginning probabilities: a sequence starting with "S" (in a soluble state) or "T" (in a trans membrane state)
    '''    
    stateSeqFile = open(stateSeqFile, "r")
      
    totalTransitionCount = 0
    startsWithSCount = 0
    startsWithTCount = 0
    startLinesCounted = 0
    sTransitionCount = 0
    tTransitionCount = 0
   
    for seq in stateSeqFile:
        if(seq[0] == "S"):
            startsWithSCount += 1
        if(seq[0] == "T"):
            startsWithTCount += 1
            
        startLinesCounted +=1
        seq = seq.rstrip()
        k = int(k)
        
        for i in range(0, len(seq)-1): 
            seqPair = seq[i:i+k]              
            
            if(seq[i] == "S"):
                sTransitionCount +=1
            if(seq[i] == "T"):
                tTransitionCount += 1    
            if(seqPair in STATEDICT):
                STATEDICT[seqPair] += 1 
                totalTransitionCount +=1 #transitions
            else:
                STATEDICT[seqPair] = 1
                totalTransitionCount +=1

    p_startsWithS = startsWithSCount/(startLinesCounted)
    p_startsWithT = startsWithTCount/(startLinesCounted)

    p_ss = (STATEDICT.get("SS")/(sTransitionCount))
    p_st = (STATEDICT.get("ST")/(sTransitionCount))
    p_tt = (STATEDICT.get("TT")/(tTransitionCount))
    p_ts = (STATEDICT.get("TS")/(tTransitionCount))

    stateSeqFile.close()
    return(p_startsWithS, p_startsWithT, p_ss, p_st, p_tt, p_ts)

def toString(): 
    print("Frequency of Soluble Sequence File: " )      
    getFrequencyS(solubleSeqFile)
    print("Frequency of Transmembrane Sequence File: ")
    getFrequencyT(transmembraneSeqFile)
    print("Transmission probability of State Sequence File with default transition size = 2: ")
    getTransProb(stateSeqFile, int(sizeOfTransmissionState))
    print("Emission Probability of State Sequence File: ")
  
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 'Counts amino acid frequencies, state probabilities, transition probabilities (including start probabilities)')
    parser.add_argument('-sol', dest ='soluble', required = True, help = "need soluble sequences")
    parser.add_argument('-mem', dest = 'transmembrane', required = True, help = "need transmembrane sequence")
    parser.add_argument('-state', dest = 'state', required = True, help = "need state sequences")
    parser.add_argument('-k', dest = 'k', default = 2, help = "size of Transition state, default = 2" )
    parser.add_argument('-test', dest = 'testFile', required = True, help = "test file transmembrane, soluble seqs")
    args = parser.parse_args()
    
    solubleSeqFile = args.soluble
    transmembraneSeqFile = args.transmembrane
    stateSeqFile = args.state
    sizeOfTransmissionState = args.k
    testFile = args.testFile

STATEDICT = {}     
solFrequencyDictionary = getFrequencyS(solubleSeqFile)
transFrequencyDictionary = getFrequencyT(transmembraneSeqFile)
p_startsWithS, p_startsWithT, p_ss, p_st, p_tt, p_ts = getTransProb(stateSeqFile, sizeOfTransmissionState)
print("State Probability for sequence: ", stateProb(testFile))
