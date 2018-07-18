'''
Created on Apr 25, 2018

@author: Caley Thomasson
Protein Clustering

'''
import argparse
import numpy as np
from array import array

def getAADictandVectors(seqFile): #get Amino Acid Frequency Dictionary and Protein Vectors
   
    seqFile = open(seqFile, "r")
    vecList = [] #vectors of count per protein [[protein1_location, count of A in protein 1, count of B in protein 1...][protein2_location, count of A in protein 2, count of B in protein 2...]..]
    vecFreqList = [] #***Frequency vectors : [[protein1_location , Frequency of A , Frequency of B, ...][protein2_location , Frequency of A , Frequency of B ...]...]
    
    aaDictionary = {'A': 0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':0 , 'I':0 , 'K':0, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'V':0, 'W':0, 'Y':0}
    count = 0
    outerVectorIndex = 0
    for line in seqFile: 
        for aa in line:
            if ">"in line: #get name_location for protein sequence
                addingCurrentVector= [0]*21
                addingCurrentVector[0] = line
                break
            elif ">" not in line:
                count += 1
                if(aa in aaDictionary): #add to total frequency of aa in dictionary
                    aaDictionary[aa] += 1
                    currentCountOfAAInVector = addingCurrentVector[getIndexofAA(aa)]
                    addingCurrentVector[getIndexofAA(aa)] = currentCountOfAAInVector +1
        if(">" in line):
            vecList.insert(outerVectorIndex,addingCurrentVector)
            outerVectorIndex +=1
    
    for i in range(len(vecList)):
        vecList[i][0] = 'PROTEIN' + vecList[i][0]
#vector = frequencies (this divides the count of amino acids currently in the vector 
#by total number of that amino acid in the file to get frequencies for each protein        
    for proteinVector in vecList:
        vecFreqList.append(proteinVector[0])
        for i in range(1, 21):  
            temp_totalCountOfSpecificAA = aaDictionary.get(getAAbyIndex(i))
            vecFreqList.append(proteinVector[i]/temp_totalCountOfSpecificAA)
    return vecFreqList

def getAAbyIndex(index):
    aalist = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']    
    return aalist[index-1]

def getIndexofAA(aa):
    aaDictionary = {'A': 1, 'C':2, 'D':3, 'E':4, 'F':5, 'G':6, 'H':7 , 'I':8 , 'K':9, 'L':10, 'M':11, 'N':12, 'P':13, 'Q':14, 'R':15, 'S':16, 'T':17, 'V':18, 'W':19, 'Y':20}
    return(aaDictionary[aa])
       
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Uses amino acid frequency vectors for clustering')
    parser.add_argument('-aa', dest ='aaseqs', required = True, help = "need aa sequences")
    parser.add_argument('-cof', dest = 'cofractionData', required = True, help = "need cofraction data")
    parser.add_argument('-phylo', dest = 'phylo', required = True, help = "need phylo profiles")
 #   parser.add_argument('-out', dest = 'morpheus', retuired = False, help =  "need output file name ")
    args = parser.parse_args()
    
    seqFile = args.aaseqs
    cofractionData = args.cofractionData
    phyloProfiles = args.phylo


vectorList = getAADictandVectors(seqFile)
with open('morpheus', 'w') as f:
    for line in vectorList:
        line = str(line)
        if 'PROTEIN' in line:
            f.write('\n')
        f.write((str(line)).rstrip()+ '\t')
        
