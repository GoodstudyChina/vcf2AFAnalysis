#!/usr/bin/python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy 
import pdb
import pylab
from optparse import OptionParser
import os
import pandas
import sys
from scipy import interpolate
import pkg_resources
import matplotlib as mpl
import SNVer2SNPv2 
from scipy.stats import kde
import matplotlib.pyplot as plt

# TODO
# scaffoldenden plotten von agp (als punkte)




def extendIntervalDownstream(SNPDataFrame, position, pool, threshold, outlierTolerance):


    end = position

# extend downstream
    extend = True
    errorCount = 0            

    while extend == True:
                    nextSNP = SNPDataFrame[(SNPDataFrame['Position'] > end)].head(1)
                    if (nextSNP.iloc[0]['Frequency' + pool]  <= threshold) and errorCount < outlierTolerance: # next AF okay
                        end = nextSNP.iloc[0,1]
                        errorCount = 0
                        #print 'correct frequency until ' + str(end)
                    elif (nextSNP.iloc[0]['Frequency' + pool]  >= threshold) and errorCount < outlierTolerance: # if frequency to low, go on
                        searchPosition = nextSNP
                        while (searchPosition.iloc[0]['Frequency' + pool]  >= threshold) and errorCount < outlierTolerance: # look as far as outlier tolerance allows
                            #trying to extend interval to next position
                            if (SNPDataFrame[(SNPDataFrame['Position'] > searchPosition.iloc[0,1])].head(1).iloc[0]['Frequency' + pool]  <= threshold): # nextAF to low
                                end = searchPosition.iloc[0,1]
                                break


                            searchPosition = SNPDataFrame[(SNPDataFrame['Position'] > searchPosition.iloc[0,1])].head(1)
                            errorCount += 1
                            #print errorCount
                            
                        
                        
                        if (SNPDataFrame[(SNPDataFrame['Position'] > end)].head(1).iloc[0]['Frequency' + pool]  <= threshold):
                            errorCount = 0
                            #print 'reseting counter'
                        
                        
                    else:
                        extend = False
                        

    return end                        


def extendIntervalUpstream(SNPDataFrame, position, pool, threshold, outlierTolerance):


                start = position

     #extend upstream
                extend = True
                errorCount = 0



                while extend == True:
                    #print "start: " + str(start)
                    previousSNP = SNPDataFrame[(SNPDataFrame['Position'] < start)].tail(1)
                    #print previousSNP
                    #print "AF:"
                    #print previousSNP.iloc[0]['Frequency' + pools[0]]
                    #print "error count: " + str(errorCount)
                    if (previousSNP.iloc[0]['Frequency' + pool]  <= threshold) and errorCount < outlierTolerance:
                        start = previousSNP.iloc[0,1]
                        errorCount = 0
                        #print 'correct frequency until ' + str(start)
                    elif (previousSNP.iloc[0]['Frequency' + pool]  >= threshold) and errorCount < outlierTolerance:
                        searchPosition = previousSNP
                        while len(searchPosition.index)>0 and (searchPosition.iloc[0]['Frequency' + pool]  >= threshold) and errorCount < outlierTolerance: # falls der errorCounter noch nicht zu hoch ist wird weiter verlängert
                            #print 'trying to extend interval to ' + str(searchPosition.iloc[0,1])
                            #print 'AF ' + str(searchPosition.iloc[0]['Frequency' + pools[0]])
                            #print 'nextAF ' + str(SNPDataFrame[(SNPDataFrame['Position'] < searchPosition.iloc[0,1])].tail(1).iloc[0]['Frequency' + pools[0]])
                            if len(SNPDataFrame[(SNPDataFrame['Position'] < searchPosition.iloc[0,1])].index) != 0 and (SNPDataFrame[(SNPDataFrame['Position'] < searchPosition.iloc[0,1])].tail(1).iloc[0]['Frequency' + pool]  <= threshold):
                                start = searchPosition.iloc[0,1]
                                break


                            searchPosition = SNPDataFrame[(SNPDataFrame['Position'] < searchPosition.iloc[0,1])].tail(1)
                            errorCount += 1
                            #print "Error Count"
                            #print errorCount
                            
                        
                        
                        if (SNPDataFrame[(SNPDataFrame['Position'] < start)].tail(1).iloc[0]['Frequency' + pool]  <= threshold):
                            errorCount = 0
                            #print 'reseting counter'
                        
                        
                    else:
                        extend = False
 
                return start




def nDI(SNPDataFrame,threshold, pools, poolNumber, outlierTolerance, filterIntervalsSmallerThan = int(1000)):
    """ Takes a DataFrame with Variants, and returns a List of the longest stretches of Variants equal to or above the threshold. Es werden erst Startpunkte anhand der deltaAF
gesucht. Erweitert werden die Intervalle dann so lange, wie die AF des ersten pools stimmt, dieser also nur ein Allel aufweist."""

    print "searching for intervals"
    # testing if threshold is valid
    if not 0.0 <= threshold and threshold <= 1.0:
        print " Threshold not within range. setting to 0.1"
        threshold = 0.1
    print "Threshold for outliers: " + str(threshold)
    listOfIntervals = []
    filterCount = 0
    listOfFilteredIntervals = []
    #"ersten SNP suchen mit entsprechender AF. So lange weitergehen, bis x SNPs nacheinander (aus dem confidence Interval; binominalverteilung mit coverage um AF), unter den
    #threshold fallen"

    # finding startpoints for interval calculation
    #mean_dAF = (absDelta(SNPDataFrame['Frequency' + pools[0]],SNPDataFrame['Frequency' + pools[1]])).mean()
    sequencingErrorRate = 0.01 # sequencing error rate according to sequencing method
    print "assumed sequencing error rate: " + str(sequencingErrorRate)
    phenotypicDifference = 1 # difference in pools should ideally be 1 = 100%
    print "assumed phenotypic difference between pools: " + str(phenotypicDifference)
    # seeds of intervalls; delta value = 1 - plausible error rate
    print "finding seeds of intervals"

    listOfPutativeStartPositions = SNPDataFrame[(abs(SNPDataFrame['Frequency' + pools[0]] - SNPDataFrame['Frequency' + pools[1]]) >= (phenotypicDifference - 100. * sequencingErrorRate / (SNPDataFrame['Coverage' + pools[0]] + SNPDataFrame['Coverage' + pools[1]])/2)) & (SNPDataFrame['Frequency' + pools[poolNumber]] <= threshold)]['Position'].tolist() 

    pool = pools[poolNumber]
    print "calculating intervals for pool " + pool


    if len(listOfPutativeStartPositions) > 0:
        # initialising
        print "extending seeds"
        start = 0
        end = 0
        length = len(listOfPutativeStartPositions)
        i = 0
        
    # extend the interval as far as possible
        while i < length: # as long as there are still seeds to check
            if start == 0: # initialising new search
                start = listOfPutativeStartPositions[i]
                end = start
                errorCount = 0
            else:
                # extend downstream
                end = extendIntervalDownstream(SNPDataFrame, end,pool, threshold, outlierTolerance)
                #extend upstream
                start = extendIntervalUpstream(SNPDataFrame, start,pool, threshold, outlierTolerance)
                   
                
######## filtering intervals for size
                intervalLength = int(end) - int(start) + 1
                if intervalLength < filterIntervalsSmallerThan:
                       # "storing filtered Interval"
                       listOfFilteredIntervals.append((int(start),int(end),intervalLength))
                       start = 0
                       filterCount+=1
                       i += 1
                else:
                       # "appending interval"
                       listOfIntervals.append((int(start),int(end),intervalLength))
                       while  i < len(listOfPutativeStartPositions) and end >= listOfPutativeStartPositions[i]:
                           # "skip seeds that overlap the interval"
                           i += 1
                       start = 0

    return listOfIntervals, filterCount, listOfFilteredIntervals


### ToDo In einen plot, effizienter
def plotIntervals(chromosome, listOfIntervals, ax, color = 'm'):
    for interval in listOfIntervals:
            start, end, length = interval
            ax.plot([int(start), int(end)],[1.1,1.1],color+'.-')
    return
    


def gff2ExonList(lines):
        '''Takes the lines of a gff file; returns a list of lists, with a list ['genename',start_exon1,end_exon1,start_exon2,end_exon2,...]  for each gene'''
        # initialize
        print "reading exons from gff"
        list_of_genes = []
        exons = None
        inGene = False
        for line in lines:
                line = line.strip().split()
                if line[2] == 'mRNA' and inGene == False:
                        exons =[]
                        exons.append(line[8].split(';')[0])
                        inGene = True
                elif (line[2] == 'exon' or line[2] == 'CDS') and inGene == True:
                        
                        exons.append(line[3])
                        exons.append(line[4])
                elif line[2] == 'mRNA' :
                        list_of_genes.append(exons)
                        exons = []
                        exons.append(line[8].split(';')[0])
                        inGene = True
 #               else:
#                        continue
                
        if exons != None:
                list_of_genes.append(exons)

        print len(list_of_genes)
        print list_of_genes 
               
        return list_of_genes
               
### ToDo in einen plot, effizienter
def plotExonList(list_of_genes, ax):
        print "plotting exons"
        for gene in list_of_genes:
          
            for pos in range(len(gene)):

                if pos == 0:
                    genename = str(gene[0]) #plot genename
                    ax.text(gene[1],-0.05,genename.replace("ID=",""))
                elif pos%2 == 1: # plot exon
                    ax.plot([gene[pos], gene[pos+1]],[-0.1,-0.1],'b')
                    
                elif pos%2 == 0 and pos != len(gene)-1: #plot intron
                    ax.plot([gene[pos], gene[pos+1]],[-0.1,-0.1],'m.-')
            
                else :
                     continue

        return

def plotgff(gff_file, Chr, ax):
    filehandle = open(gff_file,'r')
    lines = [x for x in filehandle.readlines() if Chr in x]
    for line in lines:
            line = line.split()
            ax.plot([int(line[3]), int(line[4])],[-0.1,-0.1],'m.-')

    filehandle.close()
    return



### erweitern für agp file
def makeListOfContigPositions(contigFile, Ns, Chr):
        """Chromosomal References are normally put together by adding up all scaffolds and contigfs in the (presumably)
        right order, adding a defined number of Ns inbetween.  This script calculates, given a list of contigs in fasta format, and depending on the
         number of separating Ns, the start and/prj/gf-gabibeet/data/Assemblies/RefBeet/RefBeet-1.0/RefBeet-1.0.release.contigs.fna end positions of every contig in the reference and returns
         them as a list of tupels [(start, end),...]"""
        f = open(os.path.abspath(contigFile),"r")
        position = int(0)
        contigPos = []
        count =0

        for line in f:
            if line.startswith(">") and Chr in line: # if the line is a fasta header line ... 
                line = f.next() # read in the next line, since it must be the sequence
                start = position
                end = position + len(line)
                contigPos.append((start,end))
                position = (end + Ns)
                count += 1
       
        print  str(count) + ' contigs for '  + Chr +  ' read'
        return contigPos         


def readSNVerSNPfile2dataFrame(path, namesOfPools, chromosome = 'all'):
        """reads in the tabdelimited txt file containing the SNPs and returns a dataframe. This file has to be created out of SNVer output using SNVer2SNP.py. The names of the pools have to be provided
        as a list in the right order (same order as in the file which is to be read in). The returned dataframe object can be large, since the whole data is stored and only later
        filtered for chromosomes, so you might consider prefiltering your input file for the chromosomes you want to analyze. If a list of specified chromosomes is given, Only one chromosome at a time is read in, resulting in less memory
        requirement, but slower processing"""
        # defining how many columns have to be read in, and their datatypes
        names = ['Chromosome', 'Position', 'Ref', 'Alt'] # Standard fields
	
        for name in namesOfPools: # naming the three fields for each pool
                names.append('AltCount' + name)
                names.append('Coverage' + name)
                names.append('Frequency' + name)
	
        formats = ['S30', 'i', 'S14', 'S14'] # Standard fields datatypes

        for name in namesOfPools: # defining the datatypes of the three fileds for each pool
                formats.append('i')
                formats.append('i')
                formats.append('f')

        cols = [0,1,2,3] +  range(5, len(namesOfPools)*3 +5) # which columns to read in from the txtfile. Standard fields + three fields for each pool

        
        print 'Columns of the file to be read in:'
        print cols
        print 'Naming of the columns:'
        print names

	filehandle = open(path,'r')

#        from IPython import embed
#        embed()
	
	if chromosome == 'all':
            dataFrame = pandas.DataFrame(numpy.loadtxt(filehandle, dtype={'names': names,'formats': formats},usecols=cols))
        else: # reading in chunkwise, filtering each chunk for the desired chromosome
            iter_txt = pandas.read_table(filehandle, dtype={'names': names,'formats': formats}, usecols=cols, iterator=True, chunksize=1000, header = 0 , names = names)
            dataFrame = pandas.concat([chunk[chunk['Chromosome'].isin(chromosome)] for chunk in iter_txt])
	filehandle.close()
        print "read: " + path
        
        print "created dataFrame" 
        
        dataFrame.dropna()
        print "dropped NA values"	

        return dataFrame



def readMarkerfile(path, chromosome = 'all'):
	"""read in a Markerfile as dataframe. All but the first and second columns are ignored. the first column is expected to be the Chromosome names, the second to be the positions"""
        names = ['Chromosome', 'Position'] # Standard fields
        formats = ['S30', 'i'] # Standard fields datatypes
        cols = [0,1]
	
	filehandle = open(path,'r')
	
	if chromosome == 'all':
            dataFrame = pandas.DataFrame(numpy.loadtxt(filehandle, dtype={'names': names,'formats': formats},usecols=cols))
        else: # reading in chunkwise, filtering each chunk for the desired chromosome
            iter_txt = pandas.read_table(filehandle, dtype={'names': names,'formats': formats}, usecols=cols, iterator=True, chunksize=1000, header = 0 , names = names)
            dataFrame = pandas.concat([chunk[chunk['Chromosome'].isin(chromosome)] for chunk in iter_txt])
	filehandle.close()
        print "read: " + path

        dataFrame.dropna()
        print "dropped NA values"	

        return dataFrame
    

def checkOrCreateSNPfile(path_VCF):
    if os.path.isfile(path_VCF + ".snp") and not os.stat(path_VCF + ".snp")[6]==0:
            path_SNPs = path_VCF + ".snp"
            print "SNP file already exists"
    elif os.path.isfile(path_VCF + ".snp") and os.stat(path_VCF + ".snp")[6]==0:
            print "existing SNP file is empty. removing"
            os.remove(path_VCF + ".snp")
            print "recreating SNP file"
            path_SNPs = path_VCF + ".snp"
            path_SNPs =  SNVer2SNPv2.parseVCF('GATK', path_VCF, path_SNPs)
            print "created SNP file"
    elif not os.path.isfile(path_VCF + ".snp"):
            print "SNPfile doesnt exist, creating ..."
            path_SNPs = path_VCF + ".snp"
            path_SNPs =  SNVer2SNPv2.parseVCF('GATK', path_VCF, path_SNPs)
            print "created SNP file"
    else:
            print "error: Problem with input "
    
    return path_SNPs


def overlapPositions(df_SNPs1, df_SNPs2):
	"""takes two dataframes with Variations. Returns the Variations common to both dataframes. """
	print "Finding overlap of dataframes"
#	df_SNPs = df_SNPs.join(df_Markers, how='inner', on=['Chromosome', 'Position'])
#	df_SNPs = pandas.concat([df_SNPs, df_Markers], join='inner').groupby(['Chromosome', 'Position'], as_index=False)
	overlap_SNPs = pandas.merge(df_SNPs1, df_SNPs2, how='inner', on=['Chromosome', 'Position'])
	print "Results:"
	print "first file:"
	print df_SNPs1.describe()
	print "second file:"
	print df_SNPs2.describe()
	print "overlap:"
	print overlap_SNPs.describe()
	return overlap_SNPs

def overlapVCF(path1, path2):
	"""return and save the Variants common in two vcf or snp files"""
	vcf1 = readMarkerfile(path1)
	vcf2 = readMarkerfile(path2)
	overlap = overlapPositions(vcf1, vcf2)
	print "writing overlap to:"
	print "overlapVariants.csv"
	overlap.to_csv("overlapVariants.csv", sep='\t', index=False)
	return overlap
	
	
	
def movingAverageOverDataPoints(interval, window_size):
	'''Sliding window over windowsize number or datapoints, for example window_size 20 means 20 datapoints. Die stepSize ist immer 1 und es werden genausoviele Datenpunkte
	zurueckgegeben, wie reingegeben. Es ist also Kein sliding window wo man datenpunkte zusammenfasst, sondern ein smoothing ueber eine angegebene anzahl Datenpunkte.'''
    	window = numpy.ones(int(window_size))/float(window_size)
    	return numpy.convolve(interval, window, 'same')	
		
#### bottleneck ist die size bestimmung der dataframes, das geht doch schneller über index und index.len
def movingAverageOverWindow(dataFrameChr, windowSize, stepSize): # for snp density # langsam, da der dataframe fuer jedes window nach positionen gefiltert wird # bottleneck im ganzen tool. entweder positionen als index nutzen, oder die endposition des letzten filterns merken und als start fuer das naechste nutzen? 
	'''Sliding Window over a region of windowSize, with a stepSize. windowSize 20000 means a chromosomal part of size 20kb. Returns a list of the number of datapoints in each interval, and the middle point of the interval.
the return values are a list of numbers of snps per interval, and a correlating list of positions, which mark the middle points
of each interval.'''
	numbers=[]
	positions=[]

	dataPoints = numpy.arange(1, dataFrameChr['Position'].irow(-1)+windowSize,stepSize) # create a number of datapoints defining the intervals
	for point in range(len(dataPoints))[1:-1]: ## exclude the first and last datapoint, which are start and end
		number = len(dataFrameChr[ (dataFrameChr['Position']>=dataPoints[point]-1) & (dataFrameChr['Position']<dataPoints[point+1])]) # count variations in window
		numbers.append(number)
		positions.append(dataPoints[point])
	return numbers, positions

def SlidingWindow(dataFrameChr, windowSize, stepSize, columnName): # langsam, da der dataframe fuer jedes window nach positionen gefiltert wird # bottleneck im ganzen tool. entweder positionen als index nutzen, oder die endposition des letzten filterns merken und als start fuer das naechste nutzen?
	'''Sliding Window over a region of windowSize, with a stepSize. windowSize 20000 means a chromosomal part of size 20kb. Returns a list of the mean value of the named column in each interval, and the middle point of the interval'''
	numbers=[]
	positions=[]
	dataPoints = numpy.arange(1, dataFrameChr['Position'].irow(-1)+windowSize,stepSize) # create a number of datapoints defining the intervals
	for point in range(len(dataPoints))[1:-1]: ## exclude the first and last datapoint, wich are start and end
		number = dataFrameChr[columnName][ (dataFrameChr['Position']>=dataPoints[point]-1) & (dataFrameChr['Position']<dataPoints[point+1])].mean() # count variations in window
		numbers.append(number)
		positions.append(dataPoints[point])
	return numbers, positions

def smoothTriangle(data,degree,dropVals=False):
        """performs moving triangle smoothing with a variable degree."""
        """note that if dropVals is False, output length will be identical
        to input length, but with copies of data at the flanking regions"""
        triangle=numpy.array(range(degree)+[degree]+range(degree)[::-1])+1
        smoothed=[]
        for i in range(degree,len(data)-degree*2):
                point=data[i:i+len(triangle)]*triangle
                smoothed.append(sum(point)/sum(triangle))
        if dropVals: return smoothed
        smoothed=[smoothed[0]]*(degree+degree/2)+smoothed
        while len(smoothed)<len(data):smoothed.append(smoothed[-1])
        return smoothed


def delta(AFdata1, AFdata2):
        delta = []
        if len(AFdata1)==len(AFdata2): # es muessen gleichviele Datenpunkte sein. Es wird NICHT geprÃ¼ft, ob die Daten an der selben Position sind
                #for i in range(len(AFdata1)):
                #        delta.append(sum(AFdata1[i:i+1])-sum(AFdata2[i:i+1])) # Ist es hier besser, absolutwerte zu nehmen???
		delta = AFdata1 - AFdata2
        
        return delta
		
def absDelta(AFdata1, AFdata2):
        '''Takes to series or columns of dataframes, returns a series with absdelta values'''
        if len(AFdata1)==len(AFdata2): # es muessen gleichviele Datenpunkte sein. Es wird NICHT geprÃ¼ft, ob die Daten an der selben Position sind
            delta = abs(AFdata1 - AFdata2)
        return delta



def splineCurveFit(x,y,smoothness=25,spline=5,nest=-1):
        # a spline is a sufficiently smooth polynomial function that is piecewise-defined, and possesses a high degree of smoothness at the places where the polynomial pieces connect (which are known as knots)
        #  implement a curve fit for the data
        #  Interpolation of an N-D curve
        from numpy import arange, cos, linspace, pi, sin, random
        from scipy.interpolate import splprep, splev

        # spline parameters
        # s=25.0 # smoothness parameter
        # k=5 # spline order
        # nest=-1 # estimate of number of knots needed (-1 = maximal)

        # find the knot points
        tckp,u = splprep([x,y],s=smoothness,k=spline,nest=-1)

        # evaluate spline, including interpolated points
        xnew,ynew = splev(linspace(0,1,400),tckp)
        return xnew,ynew

# Define model function to be used to fit to the data above:
def gauss(x, *p):
        A, mu, sigma = p
        return A*numpy.exp(-(x-mu)**2/(2.*sigma**2))

def calculateFittedGaussianCurve(positions, data):
                # takes a list of positions and corresponding values, returns fitted gaussian curve
                # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
                p0 = [1., 0., 1.]
                from scipy.optimize import curve_fit
                coeff, var_matrix = curve_fit(gauss, positions, data, p0=p0)
                print "Gaussian coefficients: "
                print 'Fitted mean = ', coeff[1]
                print 'Fitted standard deviation = ', coeff[2]
                # Get the fitted curve
                fit = gauss(positions, *coeff)
                return fit, coeff[1], coeff[2]


def filterDataFrameForChromosome(dataFrame,Chrs):
    '''filters the dataframe for all the chromosomes in Chrs, takes a string or list'''
    #	dataFrame = dataFrame[dataFrame['Chromosome']==Chrs]
    if type(Chrs) == list:
        return dataFrame[dataFrame['Chromosome'].isin(Chrs)]
    else:
        return dataFrame[dataFrame['Chromosome'].isin([Chrs])]

def filterDataFrameForCoverage(dataFrame,coverage):
        
        coveragecolumns = []
        for column in dataFrame.columns.tolist():
                if 'Coverage' in column: # select columns that contain the word Coverage
                        coveragecolumns.append(column)

        
        for column in coveragecolumns:
                dataFrame = dataFrame[(dataFrame[column]>coverage)&(dataFrame[column]<(coverage*3))]


        print "SNPs after coverage filtering Cov: " + str(coverage) + " to " + str(coverage*3)
        print dataFrame.describe()
        print "\n \n \n"
	
        return dataFrame
          




def filterDataFrameForNoneInformativeSNPs(dataFrame, noninformativeCutOff):
        '''Filter out SNPs, that have a very similar  AF'''
        print "Filtering noninformative SNPs"
        #print dataFrame.describe()

        Statisticsfile.write("SNPs before filtering: " + " \n")
        Statisticsfile.write(str(dataFrame.describe())+ " \n")        

        frequencycolumns = []
        for column in dataFrame.columns.tolist():
                if 'Frequency' in column: # select columns that contain the word frequency
                        frequencycolumns.append(column)

        dataFrame = dataFrame[abs(dataFrame[frequencycolumns[0]] - dataFrame[frequencycolumns[1]]) >= noninformativeCutOff] 
                
        Statisticsfile.write("SNPs after filtering for uninformative SNPs: " + " \n")
        Statisticsfile.write(str(dataFrame.describe())+ " \n") 		
        print "After filtering"
        print dataFrame.describe()
        
        return dataFrame

def filterZeros(dataFrame):
        '''removes als variants that have 0 allelefrequency in both pools, and are thus noninformative'''
        frequencycolumns = []
        for column in dataFrame.columns.tolist():
                if 'Frequency' in column: # select columns that contain the word frequency
                        frequencycolumns.append(column)

        
        dataFrame[(dataFrame[frequencycolumns[0]] != 0) & (dataFrame[frequencycolumns[1]] != 0.0)]
        return dataFrame	


def filterDataFrameForProportionalSNPs(dataFrame,proportional):
        '''the sum of both allelefrequencies must be inbetween 1 - x < sumAF < 1 + x. x can range from 0 to 1. Thus, the higher x, the more are the pools allowed to differ.'''
        print "Before filtering for proportional SNPs"
        print dataFrame.describe()

        frequencycolumns = []
        for column in dataFrame.columns.tolist():
                if 'Frequency' in column: # select columns that contain the word frequency
                        frequencycolumns.append(column)

        

        dataFrame = dataFrame[(dataFrame[frequencycolumns[0]] + dataFrame[frequencycolumns[1]] > (1 - proportional)) & (dataFrame[frequencycolumns[0]] + dataFrame[frequencycolumns[1]] < (1 + proportional))]
                
		
	
        print "After filtering for proportional SNPs"

        print dataFrame.describe()
        return dataFrame


def readAndFilterVariants(path_SNPs, pools, df_Markers, df_errors, Chrs, noninformative, coverage, proportional, Path2StatisticsFile, low_memory = False):
        
    print "processing"
    print path_SNPs

    if low_memory == True:
        df_SNPs = readSNVerSNPfile2dataFrame(path_SNPs, pools, Chrs) # read in the SNP file as a dataframe, only for the requested chromosomes
    else:
        df_SNPs = readSNVerSNPfile2dataFrame(path_SNPs, pools) # read in the SNP file as a dataframe, for all chromosomes
    stats = df_SNPs.describe()
    print "Unfiltered data: \n"
    print stats
    Statisticsfile = open(Path2StatisticsFile, 'w')  
    Statisticsfile.write("Unfiltered data: \n")
    Statisticsfile.write(str(stats)+ " \n")

    df_SNPs = filterZeros(df_SNPs)

    ## get Chromosomes 
#    df_Chrs = []
#    print "fetching Chromosome Names from DataFrame"
#    searchArray = df_SNPs['Chromosome'].unique()
#    print searchArray.dtype.names
#    for i in range(len(searchArray)):
#            df_Chrs.append(searchArray[i])
#    df_Chrs.sort()
#    print 'Chromosomes found in Data: '  
#    print df_Chrs
    df_Chrs = df_SNPs['Chromosome'].unique().tolist()




    if coverage == 0:
            dataFrame = pandas.DataFrame()
            for Chr in df_Chrs:
                    chromosome = []
                    chromosome.append(Chr)
                    chrSNPs = filterDataFrameForChromosome(df_SNPs,chromosome) # only keep SNPs on Chromosome Chr
                    
                    for name in pools:
                            #from IPython import embed
                            #embed()
        
                            medianCoverage = chrSNPs['Coverage'+name].median()
                            dev = chrSNPs['Coverage'+name].std()
                            print "automatic coverage filtering for " + Chr + " pool: " + name + " with median coverage = " + str(medianCoverage) + " ; filter from " + str(medianCoverage*0.75) + " to " + str(medianCoverage*0.75*2)
                            #chrSNPs = filterDataFrameForCoverage(chrSNPs, medianCoverage-dev) # filter each pool seperately for coverage
                            chrSNPs = filterDataFrameForCoverage(chrSNPs, medianCoverage*0.75) # filter each pool seperately for coverage

                    if dataFrame.empty:
                            dataFrame = chrSNPs
                    else:
                            dataFrame = pandas.concat([dataFrame,chrSNPs])
            dataFrame.describe()
            df_SNPs = dataFrame
            Statisticsfile.write(str(dataFrame.describe())+ " \n")


    elif coverage:
            df_SNPs = filterDataFrameForCoverage(df_SNPs, coverage)


    if type(df_Markers) != bool:
            print "Filtering for Markers"
            df_SNPs = pandas.merge(df_SNPs, df_Markers, how='inner', on=['Chromosome', 'Position'])
            print df_SNPs.describe()


    if type(df_errors) != bool:
        print "Filtering for Errors"
        df_errors['temp'] = 1
        df_SNPs = pandas.merge(df_SNPs, df_errors, how='outer', on=['Chromosome', 'Position'])
        df_SNPs = df_SNPs[df_SNPs['temp'] != 1] # only keep SNPs that were not in df_Errors
        del df_errors
        print df_SNPs.describe()


    if Chrs == "all":
            Chrs = df_Chrs
            print 'Chromosomes found for data Analysis: '  
            print Chrs
    elif frozenset([frozenset(element) for element in  Chrs]) <= frozenset([frozenset(element) for element in df_Chrs]):
            print 'Chromosomes found for data Analysis: '  
            print Chrs
    else:
            print 'Chromosomes you want are not in your data. Will use available Chromosomes instead'
            Chrs = df_Chrs
            print 'Chromosomes found for data Analysis: '  
            print Chrs

    if noninformative:
            df_SNPs = filterDataFrameForNoneInformativeSNPs(df_SNPs, noninformative) # filter out none informative SNPs
                    
    if proportional:
            df_SNPs = filterDataFrameForProportionalSNPs(df_SNPs, proportional) # filter out none informative SNPs
            

    for column in df_SNPs.columns.tolist():
            print column

    Statisticsfile.close()
    return df_SNPs, Chrs

# plotting delta values
def plot_delta(chrSNPs, windowSize, pools, ax, labels,  color = "y", labelname = "delta frequencies" , deltaIndex = "delta" , boost_flag = False):
	# if window size is set, the window size is the length of the Interval
    
            deltaValues = absDelta(chrSNPs['Frequency'+pools[0]],chrSNPs['Frequency'+pools[1]])
            deltaDf = False
            print "calculated delta values"
            if windowSize:
		    print "Length of sliding window for movingAverage: " + str(windowSize)
		    #from IPython import embed
                    #embed()
                    if chrSNPs['Position'].tail(1).iloc[0]>=windowSize :
                            
                            deltaDf = pandas.concat([chrSNPs['Position'], deltaValues], axis=1)
                            deltaDf.columns = ['Position','delta']
                            movingAverageDeltaValues, numbers = SlidingWindow(deltaDf, windowSize, windowSize/2, 'delta')
                            deltaDf.describe()
                            deltaDf.head()


                    else:
                            print "interval longer than chromosome"
                            
                    if len(movingAverageDeltaValues) == len(numbers):
			    labels[deltaIndex], = ax.plot(chrSNPs['Position'], deltaValues, "."+color, label = labelname )
                            print "plotted delta values"
                        
                            labels[deltaIndex + ' rolling mean'], = ax.plot(numbers,list(movingAverageDeltaValues),'b', label= "delta rolling mean")
			    print 'plotting moving average for delta'
                            
			    if boost_flag == True:
				    x = deltaValues
				    max_dAF = max(x)
				    #boost = [1/(max_dAF-y) if y < max_dAF else 1000  for y in x]
				    #movingAverageOverWindow(dataFrameChr, windowSize, stepSize)
				    boost = movingAverageOverDataPoints([1/(max_dAF-y) if y < max_dAF else 1000  for y in x],windowSize) # Algorithm from 1. Schneeberger K, Ossowski S, Lanz C, Juul T, Petersen AH, Nielsen KL, Jørgensen J-E, Weigel D, Andersen SU (2009) SHOREmap: simultaneous mapping and mutation identification by deep sequencing. Nat Methods 6:550–1
				    labels['delta boost'], = ax.plot(chrSNPs['Position'],boost, 'b.')
				    
				    x = chrSNPs['Position']
				    density = kde.gaussian_kde(x)
				    xgrid = numpy.linspace(x.min(), x.max(),  len(x))#chrSNPs['Position'].irow(-1)/windowSize)
				    #labels['delta boost KDE'], = ax.plot(x,numpy.array(boost)/(density(boost)*100000000), 'r-')
				    #labels['delta boost KDE'], = ax.plot(x,1/(density(boost)*100000000), 'r-')
				    
				    
				
                    
		    else:
                            print "no moving average plotted"
                            print "number of moving average delta values: " + str(len(movingAverageDeltaValues))
                            print "number of positions: " + str(numbers)
            else:
                    labels[deltaIndex], = ax.plot(chrSNPs['Position'], deltaValues, "."+color, label = labelname )
                    print "plotted delta values"    

            del deltaDf
            return deltaValues




def plotMovingAverage(chrSNPs, pools, windowSize, ax, labels): # plot moving average for allele frequency
            for pool in pools:

                    (movingAveragePoolValues, positions) = SlidingWindow(chrSNPs, windowSize,windowSize/2 , columnName = 'Frequency'+pool)
                            # movingAverageOverWindow(chrSNPs['Frequency'+pool], windowSize, windowSize/2)
                    labels[pool + ' rolling mean'], = ax.plot(positions, list(movingAveragePoolValues),color = mpl.rcParams['axes.color_cycle'][pools.index(pool)+2],label= pool + " rolling mean")
                        #    labels[pool + ' rolling mean'], = ax.plot(chrSNPs['Position'],list(movingAveragePoolValues),color = mpl.rcParams['axes.color_cycle'][pools.index(pool)+2],label= pool + " rolling mean")


def plotSnpDensity(snpDensityWindowSize, ax, chrSNPs, labels):

# fuegt der grafik ein moving average ueber die SNP density hinzu, sowie ein dazugehoeriges kernel density estimate
    ax2 = ax.twinx()
    ax2.set_ylabel(r"Snps per Window")
    
    if chrSNPs['Position'].irow(-1)>=snpDensityWindowSize :
            print "calculating SNP density"
            (movingAveragePoolValues , positions) = movingAverageOverWindow(chrSNPs, snpDensityWindowSize, snpDensityWindowSize/2)
            
            x = chrSNPs['Position']
            density = kde.gaussian_kde(x,bw_method = 'silverman' )
            xgrid = numpy.linspace(x.min(), x.max(),  len(x))#chrSNPs['Position'].irow(-1)/windowSize)


    else:
            print "not enough values to calculate moving average for SNPs"
    if len(movingAveragePoolValues) == len(positions):
            labels['snp density'], = ax2.plot(positions,movingAveragePoolValues) #color = mpl.rcParams['axes.color_cycle'][pools.index(pool)+2],label= pool + " snp density " + str(windowSize/1000) + "kb")
            labels['snp kde'], = ax2.plot(xgrid, density(xgrid)*1000000000, 'r-')
            n, bons, patches = ax2.hist(x, bins=8, normed=True)
    else:
            print "no snp density plotted"

    return ax2


def save_plot(fig, Path2OutputFile, dpi = 600):
    #mng = plt.get_current_fig_manager()
    fig.set_size_inches(12,6)
    if fig.savefig(Path2OutputFile):
	print "saved file to: \n"
	print Path2OutputFile
