# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 10:02:36 2018
OKAY
@author: BETA COMPUTER
"""
import operator
import pylab, scipy, random
import numpy as np
from scipy.stats import shapiro
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.libqsturng import psturng
def getDataBiorad(fileName):
    dataFile = open(fileName, 'r')
    dataFile.readline()
    data = {}
    for line in dataFile:
        print(line)
        gene,sample,exp = line.split(';')
        print(exp)
        if gene not in data:
            data[gene] = [(sample, round(float(exp),3))]
        else:
            data[gene] += [(sample, round(float(exp),3))]
    dataFile.close()
    return data
data = getDataBiorad('Organotypique culture gene study CB1.csv')
sampleList = ['ACAM','ACEA', 'AM2-', 'CTRL']
#for a in data['Bax']:
#    sampleList.append(a[0])

def drawBoxplot(data, gene, sampleList):
    """ Assumes gene is the name of the
    data set of gene
        sampleList is list of sample tested"""
    geneData = data[gene]
    for sample in sampleList:
        draftData = []
        for expression in geneData:
            if sample in expression[0]:
                draftData.append(expression[1])
        print(draftData)
        if sampleList.index(sample) != 0:
            pylab.figure()
        pylab.boxplot(draftData)
        pylab.title(gene +'\n' + sample)
        pylab.ylim(0, 3)
###TEST NORMAL DISTRIBUTION for 1 gene ####     
def testNormality(data, gene, sampleList):
    """ Assumes gene is the name of the data set of 
    gene sampleList is list of sample tested"""
    geneData = data[gene]
    anova = 0
    for sample in sampleList:
        draftData = []
        for expression in geneData:
            if sample in expression[0]:
                draftData.append(expression[1])
        a,b = shapiro(draftData)
        if b > 0.05:
            print(sample + ' of ' + gene + ' is normally distributed')
        else:
            anova += 1
            print(sample + ' of ' + gene + ' is NOT normally distributed')
    if anova != 0:
        print ('''CANT'T do ANOVA test for ''' + gene)
        print()
        return False
    else:
        print ('CAN do ANOVA test for ' + gene)
        print()
        return True
#### WRAP UP FUNCTION TO FIND LIST OF NORMAL DISTRIBUTED GENES######
def findNormalGenes(data, sampleList):
    print('================ NORMAL TEST ==================')
    normalList = []
    for gene in data:
        if testNormality(data,gene, sampleList):
            normalList.append(gene)
    return normalList

normalList = findNormalGenes(data, sampleList)

#### AMONG NORMALLY DISTRIBUTED GENE #####
###FIND GENES THAT HAVE VARIANCE HOMOGENITY####
def testVarEqual(data, gene, sampleList):
    geneData = data[gene]
    draftData = {}
    for sample in sampleList:
        draftData[sample] = []
        for expression in geneData:
            if sample in expression[0]:
                draftData[sample].append(expression[1])
    a,b = scipy.stats.bartlett(draftData['ACEA'],
                         draftData['ACAM'],
                         draftData['AM2-'],
                         draftData['CTRL'])
    if b > 0.05 :
        print('Variances of ' + gene + ' are equal')
        print()
        return True
        
    else:
        print('Variances of ' + gene + ' are NOT equal')
        print()
        return False
            

### WRAP UP FUNCTION, FIND GENES MEETING 2 CONDITION###
### THESE GENES THEN BE ANALYSED BY ANOVA TEST ####
def findANOVAList(data, normalList, sampleList):
    print('========== VARIANCE EQUALITY TEST ===========')
    ANOVAList = []
    for gene in normalList:
        if testVarEqual(data, gene, sampleList):
            ANOVAList.append(gene)
    return ANOVAList

ANOVAList = findANOVAList(data, normalList, sampleList)

### FIND GENES THAT WILL NEED NON PARAMETRIC TEST####
def findNonParaList(data, ANOVAList):
    nonParaList = []
    for gene in data:
        if gene not in ANOVAList:
            nonParaList.append(gene)
    return nonParaList
nonParaList = findNonParaList(data,ANOVAList)

### ANOVA TEST FOR 1 GENE###
def ANOVA(data, gene, sampleList):
    geneData = data[gene]
    draftData = {}
    for sample in sampleList:
        draftData[sample] = []
        for expression in geneData:
            if sample in expression[0]:
                draftData[sample].append(expression[1])
    a,b = scipy.stats.f_oneway(draftData['ACEA'],
                         draftData['ACAM'],
                         draftData['AM2-'],
                         draftData['CTRL'])
    if b < 0.05 :
        print('A significant difference in ' + gene)
        print()
        return True
        
    else:
        print('DIDNOT find significant difference in ' + gene)
        print()
        return False

### KRUSKAL WALLIS TEST FOR 1 GENE ####
def kruskal (data, gene, sampleList):
    geneData = data[gene]
    draftData = {}
    for sample in sampleList:
        draftData[sample] = []
        for expression in geneData:
            if sample in expression[0]:
                draftData[sample].append(expression[1])
    a,b = scipy.stats.kruskal(draftData['ACEA'],
                         draftData['ACAM'],
                         draftData['AM2-'],
                         draftData['CTRL'])
    if b < 0.05 :
        print('A significant difference in ' + gene)
        print()
        return True
        
    else:
        print('DIDNOT find significant difference in ' + gene)
        print()
        return False

### WRAP UP####
##FIND 2 LIST OF GENE THAT RETURN SIGNIFICANT BY ANOVA and KRUSKAL WALLIS##

def snList(data, sampleList, ANOVAList = None, nonParaList = None):
    snListANOVA = []
    snListKruskal= []
    print('============ RUN ANOVA TEST =============')
    for gene in ANOVAList:
        if ANOVA(data, gene, sampleList):
            snListANOVA.append(gene)
    print('============ RUNE KRUSKAL TEST ============')
    for gene in nonParaList:
        if kruskal(data,gene, sampleList):
            snListKruskal.append(gene)
    return snListANOVA, snListKruskal

snANOVA, snKruskal = snList(data,sampleList, ANOVAList, nonParaList )
####################################
### STARTING POSTHOC PHASE#############
#### PREPARING DATA TO DO POSTHOC TEST ####
def setDataPosthoc(data,gene,sampleList):
    geneData = data[gene]
    exp = []
    group = []
    for sample in sampleList:
        for expression in geneData:
            if sample in expression[0]:
                exp.append(expression[1])
                group.append(sample)
    exp = np.array(exp)
    group = np.array(group)
    return exp, group
 
#### DO TUKEYPOSTHOC for 1 NORMAL GENE ####    
def posthocTukey(data, gene, sampleList):
    exp,group = setDataPosthoc(data, gene, sampleList)
    print('Summary post-hoc test Tukey of ' + gene)
    tukeyhsd=pairwise_tukeyhsd(endog =exp, groups = group,alpha=0.05)
    print(tukeyhsd)
    print()
    return tukeyhsd
### WRAP UP TUKEYPOSTHOC ####
def runPosthocTukey(data, snANOVA, sampleList):
    print ('============RUN POSTHOC TUKET TEST ===============')
    for gene in snANOVA:
        posthocTukey(data,gene,sampleList)
        
#runPosthocTukey(data, snANOVA,sampleList)


###PREPARE TO RUN MULTIPLE MANNWHITNEY U TEST ####
samplePairs = []
for a in range(4):
    for b in range(a +1, 4):
        samplePairs.append((sampleList[a],sampleList[b]))

### MWU TEST FOR 1 GENE ####
def MWU(data, gene, sampleList):
    geneData = data[gene]
    draftData = {}
    snDiff = []
    for sample in sampleList:
        draftData[sample] = []
        for expression in geneData:
            if sample in expression[0]:
                draftData[sample].append(expression[1])
    print('Mann Whiney Test with Bonferroni correction.\n' +\
          'GENE: ' + gene + ' with alpha = ' +\
          str(round(0.05/len(samplePairs),4)))
    for sample1, sample2 in samplePairs :
        t,p =scipy.stats.mannwhitneyu(draftData[sample1],
                                 draftData[sample2])
        if p < 0.05/len(samplePairs):
            print(sample1, sample2, 'pvalue =',
                  str(round(p, 4)), 'Reject = True')
            snDiff.append((sample1, sample2,str(round(p,4))))
            
        else:
            print(sample1, sample2, 'pvalue =', 
                  str(round(p, 4)), 'Reject = False')
    return snDiff

### WRAP UP FUNCTION #####
def runMWU(data, nonParaList, sampleList):
    print('============= POSTHOC MANN WHITNEY ==============')
    for gene in nonParaList:
        print()
        MWU(data,gene, sampleList)

#runMWU(data, snKruskal, sampleList)

#####################################
##############CREATING GRAPH##########
#### MUST FIRST PREPARE DATA #####
                
#### NEED TO GEST A LIST of MEANS AND SEMS FOR PLOTTING ####  

def getMeanAndSEM(draftData):
    """ Assumes draftData is a list of expression
    Return Mean and SEM"""
    mean = round(sum(draftData)/len(draftData),4)
    SEM = round(np.std(draftData)/(len(draftData)**0.5),4)
    return mean,SEM

def findMeansAndSEMs(data, gene, sampleList):
    geneData = data[gene]
    means= []
    SEMs = []
    for sample in sampleList:
        draftData = []
        for expression in geneData:
            if sample in expression[0]:
                draftData.append(expression[1])
        means.append(getMeanAndSEM(draftData)[0])
        SEMs.append(getMeanAndSEM(draftData)[1])
    return means, SEMs



#### TUKEY HSD DONT RETURN pVALUE #####
###NEED TO RETRIEVE PVALUE AND PREPARE DATA####
    
def getPValuesTukey(res):
    """ Assume res is the tukeyhsd result 
    from pairwise_tukeyhsd test"""
    rs = res.meandiffs/res.std_pairs
    pvalues = psturng(np.abs(rs), len(res.groupsunique), res.df_total)
    return pvalues

def appendPValues(reject,groupsunique, pvalues):
    i = 0
    snDiff = []
    pair = []
    length = len(groupsunique)
    for a in range(length):
        for b in range(a +1 , length):
            if reject[i] == True:
                snDiff.append((groupsunique[a],groupsunique[b], 
                           str(round(pvalues[i],4))))
            i += 1
    return snDiff

def tukeyPValue(data, gene, sampleList):
    res = posthocTukey(data,gene, sampleList)
    pvalues = getPValuesTukey(res)
    snDiff = appendPValues(res.reject,res.groupsunique, pvalues)
    return snDiff

#### PLOTING SIGNIFICANT DIFFERENCE BAR #####            
def snDiffPlot(dx,rankY,i,j,pValue,xInd,means,SEMs,figure):
    """i, j are indexes of xlabel, dx = abs(i-j)
    """
#    print(k,l,i,j)
    dx = abs(i-j)*1
    x = (xInd[i]+xInd[j])/2
#    y = (1+0.05*(l**1.7))*max([Y[a] for a in range(i,j+1)])
    if dx == 1:
        y = 1.2*max(SEMs[i],SEMs[j]) + max([means[a] for a in range(i,j+1)]) 
    else:
        y = (1+0.1*dx)*max([means[a] for a in range(i,j+1)]) +\
        0.07*(i+j)/2
    print (i,j)
    
    props = {'connectionstyle':'bar,fraction=0.05','arrowstyle':'-',
                 'shrinkA':1*dx,'shrinkB':1*dx,'linewidth':2}
    figure.annotate(pValue, xy=((xInd[i]*0.8+xInd[j]*0.2),y*1.07), zorder=10)
    figure.annotate('', xy=(xInd[i],y), xytext=(xInd[j],y), arrowprops=props)
    
    
#### WRAP UP PLOITING WITH BAR PLOTTING #####   
    
def barPlot(data, gene, sampleList, snDiff):
    means, SEMs = findMeansAndSEMs(data,gene,sampleList)
    dictMeans = {}
    for num,sample in enumerate(sampleList):
        dictMeans[sample] = means[num]
#    print(dictMeans)
    sortedMeans = dict(
            sorted(dictMeans.items(), key=operator.itemgetter(1)))
    for a,b in enumerate(sortedMeans):
        sortedMeans[b] = a +1
#    print(sortedMeans)
    means = np.array(means)
    SEMs= np.array(SEMs)
#    print(means,SEMs)
    xInd = np.arange(len(sampleList))
    bar_kwargs = {'width':0.5,'color':'y','linewidth':2,'zorder':5}
    err_kwargs = {'zorder':0,'fmt':'none','linewidth':2,'ecolor':'k',
                  'capsize':3}
    fig, ax = pylab.subplots()
    ax.p1 = pylab.bar(xInd, means, **bar_kwargs)
    ax.errs = pylab.errorbar(xInd, means, yerr=SEMs, **err_kwargs)
    ### PREPARE TO PLOT SIGNIFICANT DIFF BAR ####
    for pair in snDiff:
        sample1,sample2,pValue = pair
        pValue = 'p = ' + pValue
        if dictMeans[sample1] > dictMeans[sample2]:
            rankY = sortedMeans[sample1]
        else:
            rankY = sortedMeans[sample2]
        i = sampleList.index(sample1)
        j = sampleList.index(sample2)
        dx = abs(i - j)
        snDiffPlot(dx,rankY,i,j,pValue,xInd,means,SEMs,ax)
    pylab.xticks(xInd, sampleList, color='k')
    pylab.title(gene)
    pylab.ylim(ymax= 3)
    return ax

for gene in (snANOVA):
    snDiff = tukeyPValue(data, gene, sampleList)
    print(snDiff)
    barPlot(data, gene, sampleList, snDiff)
for gene in (snKruskal):
    snDiff = MWU (data, gene, sampleList)
    print(snDiff)
    barPlot(data, gene, sampleList, snDiff)        