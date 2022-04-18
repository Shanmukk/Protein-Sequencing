"""
Protein Sequencing Project
Name:
Roll Number:
"""

from dataclasses import replace

from pandas import cut
import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    g = open(filename, "r")
    g = g.read()
    line = g.replace('\n',"")
    return line


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    lst = []
    k = ['UAA', 'UAG', 'UGA']
    i = startIndex
    d = dna.replace("T", "U")
    while i < len(d):
        lst.append(d[i:i+3])
        if d[i:i+3] in k:
            break
        i += 3
    return lst


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    d1 = {}
    import json
    d = open(filename, "r")
    data = json.load(d)
    d = dict((x, k) for k, v in data.items() for x in v) 
    for key,value in d.items():
        s = key.replace("T","U")
        d1[s] = value
    return d1


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    lst = []
    for i in range(len(codons)):
        if i == 0:
            lst.append('Start')
        else:
            lst.append(codonD[codons[i]])
    return lst


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    count = 0
    i = 0
    lst = []
    dna = readFile(dnaFilename)
    while i < len(dna):
        if dna[i:i+3] == 'ATG':
            rna = dnaToRna(dna, i)
            protein = generateProtein(rna, makeCodonDictionary(codonFilename))
            lst.append(protein)
            i += 3*len(rna)
        else:
            i += 1
            count += 1
    return lst


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    lst = []
    for i in proteinList1:
        if i in proteinList2:
            if i not in lst:
                lst.append(i)
    return lst
 


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    lst = []
    for i in proteinList:
        for j in i:
            lst.append(j)
    return lst


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    d ={}
    for i in aaList:
        if i in d:
            d[i] += 1
        else:
            d[i] = 1
    return d


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    d = []
    x,y = aminoAcidDictionary(combineProteins(proteinList1)), aminoAcidDictionary(combineProteins(proteinList2))
    for i,j in x.items():
        x[i] = j/len(combineProteins(proteinList1))
        if i not in y:
            y[i] = 0
    for i,j in y.items():
        y[i] = j/len(combineProteins(proteinList2))
        if i not in x:
            x[i] = 0
        if abs(x[i] - y[i]) > cutoff:
            if i != 'Start' and i != 'Stop':
                d.append([i, x[i], y[i]])
    return d


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    lst = sorted(commonalities)
    print("The following proteins occurred in both DNA Sequences:")
    for i in range(len(lst)):
        for j in range(len(lst[i])):
            if lst[i][j] == "Start" or lst[i][j] == "Stop":
                continue
            else:
                print(lst[i][j], end=' ')
        print()
    print("The following amino acids occurred at very different rates in the two DNA sequences:")
    for i in sorted(differences):
        print(i[0],": ",round(i[1]*100, 2),"% in Seq1,",round(i[2]*100, 2),"% in Seq2")
        print()
    return 


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    lst = []
    x = combineProteins(proteinList1), combineProteins(proteinList2)
    for i in x:
        for j in i:
            if j not in lst:
                lst.append(j)
    return sorted(lst)


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    lst = []
    d = aminoAcidDictionary(combineProteins(proteinList))
    for i in labels:
        if i in d:
            lst.append(d[i]/len(combineProteins(proteinList)))
        else:
            lst.append(0)
    return lst


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    w = 0.35
    plt.bar(xLabels, freqList1, width=-w,  align='edge',label = label1)
    plt.bar(xLabels, freqList2, width= w,  align='edge',label = label2)
    plt.xticks(rotation="vertical")
    plt.legend()
    plt.show()
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    return


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    '''print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    test.week1Tests()
    print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    runWeek1()'''

    ## Uncomment these for Week 2 ##
    
    '''print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()'''
    

    ## Uncomment these for Week 3 ##
    
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    
