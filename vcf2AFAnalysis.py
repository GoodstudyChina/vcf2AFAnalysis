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
from matplotlib.legend_handler import HandlerLine2D
from AFAnalysis import *

try:
    from IPython import embed
except ImportError:
    print "IPython is not installed"

#from IPython import embed
#embed()


print 'SciPy Version installed: ' + pkg_resources.get_distribution("scipy").version
print 'Pandas Version installed: ' + pkg_resources.get_distribution("pandas").version
print 'matplotlib Version installed: ' + pkg_resources.get_distribution("matplotlib").version




# Parsing the command line options
parser = OptionParser()
parser.add_option('-N', '--noninformative', dest="noninformative", type="float", help='the level of allele frequency for that SNPs must differ to be counted as informative. Useful for vizualization purposes. By default, variants without any informative value (0 in both pools) are filtered.', default=False)
parser.add_option('-v', '--vcf1', dest="vcffile1", type="string", help='vcf file containing the Variations', default=False)
parser.add_option('-V', '--vcf2', dest="vcffile2", type="string", help='vcf file containing the Variations of another experiment. If this option is used, only the delta allele frequencies of the two files will be plotted', default=False)
parser.add_option('-g', '--genes_gff', dest="genes_gff", type="string", help='path to genes_gff file, containing genes to plot. Introns and exons will be plotted. This can be very timeconsuming for large files.')
parser.add_option('-G', '--gff', dest="gff", type="string", help='path to generic file [chr, startpos, endpos] to plot as intervals')
parser.add_option('-w', '--windowSize', dest="windowSize", type="int", help='the genomic length to calculate the rolling mean for. Step size is len/2. The rolling mean of the allele frequency will be calculated for both/all pools by default. if the -d option is set to "only", only the rolling mean of the delta allele frequenvy will be plotted.')
parser.add_option('-c', '--chromosomes', dest="Chrs", type="string", default="all", help='the Chromosomes to plot, delimited by ",", default is "all". Specifying a list of chromosomes, will speed up processing significantly, because only nescessary data is read, and the chromosome determination step is skipped. It is also possible to specify just one chromosome and a range to plot (Chr1:12345-56789)')
parser.add_option('-p', '--pools', dest="pools", type="string", help='the sequenced pools, in the same order as in the SNPs file, delimited by ",". Most functions are only available for two pools.')
parser.add_option('-P', '--proportional', dest="proportional", type="float", help='the level of allelefrequency for which SNPs must be proportional. For visualization.')
parser.add_option('-o', '--output', dest="output", type="string", help='Output filename. If not set, an automatic name will be generated', default=False)
parser.add_option('-x', '--coverage', dest="coverage", type="int", help='Coverage to filter for (coverage > x > coverage*3), "0" for automatic. The automatic filter will filter each chromosome seperately for (medianCoverage*0.75) > X > (medianCoverage*2.25). Automatic filtering is only advised for whole chromosomes.')
parser.add_option('-m', '--mean', dest="mean", help='plot mean', action="store_true")
parser.add_option('-L', '--low_mem', dest="low_mem", help='activate low-memory mode. This will use less memory, and make the script slower, but you will not get statistics about the whole dataset', action="store_true")
parser.add_option('-K', '--contigfile', dest="contigfile", type="string", help='file containing contigpositions, .fna or .agp', default=False)
parser.add_option('-M', '--Markerfile', dest="markerfile", type="string", help='file containing Markerpositions. Chromosomename in the first column, Position in the second.', default=False)
parser.add_option('-E', '--Errorfile', dest="errorfile", type="string", help='file containing Errorpositions. Chromosomename in the first column, Position in the second.', default=False)
parser.add_option('-H', '--Histogram', dest="histogram", help='plot an AF Histogram next to the plot of allele frequencies. If more than one chromosome is plottet, this option is ignored. For F2 progeny from 2 parents, one would expect a normal distribution around 0.5', action="store_true", default=False)
parser.add_option('-d', '--delta', dest="delta", type="string", help='Plot a delta of the allelefrequency between the two pools. Possible: "False","True","only"', default=False)
parser.add_option('-s', '--SNP_density_windowSize', dest="snpDensity", type="int", help='the size of the interval to evalute for the SNP measurement, step size for roling mean is interval/2')
parser.add_option('-i', '--interval', dest="interval", type="int", help="define and plot interval , Try to find intervals that have been selected in both pools, filter out intervals shorter than given int"  ,default = False)
parser.add_option('-b', '--batch', dest="batch", help='plot in batch mode. In this mode, the plot is not shown, but saved to disk directly', action="store_true", default=False)
parser.add_option('-B', '--Boost', dest="boost", help='plot the boosted delta allele frequencies. The algorithm is similar to the one used by SHOREmap', action="store_true", default=False)
parser.add_option('-f', '--filteredCSV', dest="filteredCSV", help='if set, the Variations after filtering will be writen to an external file', action="store_true", default=False)



## TODO:

# Werte nur einmal aus dem array auslesen und zwischenspeichern zb. positions = chrSNPs['Position']
# filter fuer kollabierte Bereiche: Zusammenhanegende Bereiche, deren Coverage zu hoch , aber deren AF zu niedrig ist rausfiltern



### Usage examples:
# -d only -w 1000 , only plot a rolling men for the delta values
# -d only -w 1000 -p green,red -v 360Acc.vcf -V 120Acc.vcf, plot the rolling mean for the delta values of two experiments, that have the same pools, but might differ in other parameters, like coverage or accessions used. 
# if you run into memory problems, consider spliting your input file (.vcf or .snp) in one file per chromosome, or activating low-memory mode
# -c Chr1 -H plot the AFE for Chr1, and a histogram of the AFE side by side. This is useful to see systematic bias or a skew in afe between the pools. In a crossing experiment of two parents, the F2 should both show a normal distribution of AFE aroud 0.5.
# -L -c BmChr3:12000-15000   --contigfile scaffold_10971_BLAST.agp  --windowSize 5000 -i 2500 plot only BmChr3 from base 12000 to 15000. Also plot scaffold and contog borders defined in the .agp file. Also plot a sliding window for the mean allelefrequeencies and search for segregating intervals larger than 2.5 kb 


(options, args) = parser.parse_args()


### Check, if the options are compatible
# -N, -P, -d, -i, -B are only applicable with exactly two pools -p pool1,pool2
if (options.noninformative or options.proportional or options.delta or options.interval or options.boost) and  (len(options.pools.strip().split(','))!= 2):
    print "-N, -P, -d, -i, -B are only applicable with exactly two pools"
    



vcfFileList = []
ax2 = False


if options.vcffile1:
	path_VCF = os.path.abspath(options.vcffile1)
	vcfFileList.append(path_VCF)
else:
	print "No input file given"


                  
if options.vcffile2:
    	path_VCF2 = os.path.abspath(options.vcffile2)
	vcfFileList.append(path_VCF2)


		

windowSize = options.windowSize


chromosome_start = False
chromosome_end = False

Chrs = options.Chrs
if Chrs != "all":
	Chrs = Chrs.strip().split(',')
	# get start and stop if set
	if len(Chrs) == 1 and ':' in Chrs[0] and '-' in Chrs[0]:
            Chromosome, pos = Chrs[0].split(':')
            chromosome_start, chromosome_end = pos.split('-')
            chromosome_start = int(chromosome_start)
            chromosome_end = int(chromosome_end)
            print 'reading only chromosome: '
            print Chromosome
            print 'from: to:'
            print chromosome_start, chromosome_end
            Chrs = []
            Chrs.append(Chromosome)
            


            

print 'Chromosomes to plot:'
print Chrs

pools = options.pools
pools = pools.strip().split(',') 


# Filters
coverage = options.coverage
proportional = options.proportional
noninformative = options.noninformative

# create outputfilename
if not options.output:
        Path2OutputFile = "Pools-" + str(options.pools) + "_Chrs-" + str(options.Chrs) + "_coverage-" + str(coverage) 
        if proportional:
                Path2OutputFile = Path2OutputFile + "_proportional-" + str(proportional)
        if windowSize:
                Path2OutputFile = Path2OutputFile + "_RollingMean-" + str(windowSize)
        if options.markerfile:
                Path2OutputFile = Path2OutputFile + "_Markerfile-"
        if options.noninformative:
                Path2OutputFile = Path2OutputFile + "_noninformative-" + str(options.noninformative)
        if options.windowSize:
                Path2Outputfile = Path2OutputFile + "_windowSize_" +str(options.windowSize)
        if options.errorfile:
                Path2Outputfile = Path2OutputFile + "_errorfile_" 
        if options.histogram:
                Path2Outputfile = Path2OutputFile + "_histogram_" 
        if options.delta:
                Path2Outputfile = Path2OutputFile + "_delta_" +str(options.delta)
        if options.snpDensity:
                Path2Outputfile = Path2OutputFile + "_snpDensity_" +str(options.snpDensity)
        if options.interval:
                Path2Outputfile = Path2OutputFile + "_minIntLen_" +str(options.interval)
            
                
        print "output:"
        print Path2OutputFile
                       
                
	
else:
	Path2OutputFile = os.path.abspath(options.output)

if not Path2OutputFile.endswith('.png'):
	Path2OutputFile = Path2OutputFile + ".png"

Path2StatisticsFile = Path2OutputFile + ".txt"




### read marker and errorfile			

df_Markers = False
if options.markerfile:
	path_Markerfile = os.path.abspath(options.markerfile, Chrs)
	if os.path.isfile(path_Markerfile) and not os.stat(path_Markerfile)[6]==0:
		df_Markers = readMarkerfile(path_Markerfile)
		print "Markerfile read:"
		print df_Markers.describe()
		print df_Markers.head()
                Statisticsfile.write("Using Markerfile:" + "\n")
                Statisticsfile.write(df_Markers.describe())
                Statisticsfile.write( "\n")

	elif os.stat(path_Markerfile)[6]==0:
		print "Markerfile is empty"
	else:
		print "Markerfile can't be read"

df_Errors = False
if options.errorfile:
	path_Errorfile = os.path.abspath(options.errorfile, Chrs)
	if os.path.isfile(path_Errorfile) and not os.stat(path_Errorfile)[6]==0:
		df_Errors = readMarkerfile(path_Errorfile)
		print "Errorfile read:"
		print df_Errors.describe()
		print df_Errors.head()
                Statisticsfile.write("Using Errorfile:" + "\n")
                Statisticsfile.write(df_Errors.describe())
                Statisticsfile.write( "\n")

	elif os.stat(path_Errorfile)[6]==0:
		print "Errorfile is empty"
	else:
		print "Errorfile can't be read"




### check if SNP file exists, else create it
path_SNPs = checkOrCreateSNPfile(path_VCF)



### read in and filter variants

if chromosome_start and chromosome_end:
    df_SNPs, Chrs = readAndFilterVariants(path_SNPs, pools, df_Markers, df_Errors, Chrs, noninformative, coverage, proportional, Path2StatisticsFile, low_memory = True, start = chromosome_start, end = chromosome_end)

elif options.low_mem == True:
    df_SNPs, Chrs = readAndFilterVariants(path_SNPs, pools, df_Markers, df_Errors, Chrs, noninformative, coverage, proportional, Path2StatisticsFile, low_memory = True)
else:
    df_SNPs, Chrs = readAndFilterVariants(path_SNPs, pools, df_Markers, df_Errors, Chrs, noninformative, coverage, proportional, Path2StatisticsFile)


### saving filtered SNPs to disk
if options.filteredCSV:
    df_SNPs.to_csv(Path2OutputFile + "_filtered_SNPS.csv", sep='\t', index=False)

if options.vcffile2:
    path_VCF2 = options.vcffile2
    path_SNPs2 = checkOrCreateSNPfile(path_VCF2)
    df_SNPs2, Chrs2 = readAndFilterVariants(path_SNPs2, pools, df_Markers, df_Errors, Chrs, noninformative, coverage, proportional, Path2StatisticsFile)
    options.delta = 'only'
    if options.filteredCSV:
    	df_SNPs.to_csv(Path2OutputFile + "_filtered_SNPS2.csv", sep='\t', index=False)



# to save mem, df_Markers and df_Errors are deleted after usage
del df_Markers
del df_Errors




### start plotting

print "start plotting \n \n \n"




fig = plt.figure()
# starting boundaries
# left = 0.00
first_SNPs=[]
last_SNPs=[]
labels = {}
listOfAxes=[] # list in which all axis are saved
chrsToPlot = []

Statisticsfile = open(Path2StatisticsFile, 'w')	

for Chr in Chrs: # make one plot for each chromosome

        if not options.vcffile2 :
        

            print 'plotting Chromosome: ' + Chr
            chrSNPs = filterDataFrameForChromosome(df_SNPs,Chr) # only keep SNPs on Chromosome Chr
            if len(chrSNPs)==0:
                    continue
            else:
                    chrsToPlot.append(Chr)
                            

#### initialising the subplot                    
            # calculating the location of the Axes, one per chromosome
            # bottom and height are always the same, the starting point (left) depends on the number of the chromosome, the width on the totl number of chromosomes.
            # With this, all chromosome plots have the same length
            #rect_cones = [left, bottom, width, height]
            rect_cones = [0.+1./len(Chrs)*(len(listOfAxes)), 0.1, 1./len(Chrs), 0.9]
            # creating the axes and adding to a list for later reference
            listOfAxes.append(plt.axes(rect_cones))
            # selecting the latest axis from the list
            ax = listOfAxes[-1]

            
            
            first_SNPs.append(chrSNPs['Position'].irow(0)) # save the length of the Chromosome   
            last_SNPs.append(chrSNPs['Position'].irow(-1)) # save the length of the Chromosome
            print first_SNPs, last_SNPs


            Statisticsfile.write('filtering for Chromosome: ' + Chr + " \n")
            Statisticsfile.write(str(chrSNPs.describe())+ " \n\n")


##### determining color for the plot and plotting the allelefrequency data
            poolcolors = {}
            colors = "bgrcmykw"
            for pool in pools: # for each pool, the trianglesmoothed frequencies should be plotted

                    if 'red' in pool:
                            color = 'r'
                            poolcolors[pool] = color
                    elif 'green' in pool:
                            color = 'g'
                            poolcolors[pool] = color
                    else:
                            color = colors[pools.index(pool)]
                            poolcolors[pool] = color


                    if options.delta == 'only':
                            continue   # don't plot SNPs
                    else:
                            labels[pool ], = ax.plot(chrSNPs['Position'],chrSNPs['Frequency'+pool],"."+color,label=pool , alpha=0.5)
    #                labels[pool + 'smoothed'], = plt.plot(chrSNPs['Position'],smoothTriangle(chrSNPs['Frequency'+pool],3),"-"+color,label=pool + " smoothed")
     #               xnew,ynew = splineCurveFit(chrSNPs['Position'],smoothTriangle(chrSNPs['Frequency'+pool],3))
    #                labels[pool + 'fitted data'],=plt.plot(xnew,ynew,'-.'+color,label='fitted data' + pool)

            



##### modifying subplot labels and ticks
            plt.xlabel(Chr) # name the x axis of the subplot after the chromosome
            
            if len(Chrs) >1 :
                    ax.set_xticklabels([]) # don't plot x axis ticks for any of the plots
            if Chrs.index(Chr) >= 1: # don't plot y axis for subplots
                    ax.set_yticklabels([])


#### plotting delta values
            print "delta: " + str(options.delta)
            if options.delta == 'True' or options.delta == 'only':
                deltaValues = plot_delta(chrSNPs, options.windowSize, pools, ax, labels, boost_flag = options.boost)

		Statisticsfile.write('mean delta: ' + Chr + " \n")
		Statisticsfile.write(str(numpy.mean(deltaValues))+ " \n")
		Statisticsfile.write(str(numpy.std(deltaValues))+ " \n\n")
        	
            elif options.windowSize: # if there is no delta, plot the moving average for the pools
                          plotMovingAverage(chrSNPs, pools, windowSize, ax, labels)


#### plot data, that is not specific to the pools


### plot SNP density
            if options.snpDensity:
                ax2 = plotSnpDensity(options.snpDensity, ax, chrSNPs, labels)				
                

### plot mean
            if options.mean:
                if options.delta:
                # plotting the mean of AF
                        AFmeanIndex = 'AF mean smoothed'
                        chrSNPs["AF_mean"] = chrSNPs[[elem for elem in chrSNPs.columns if 'Frequency' in elem]].mean(axis=1)
                        AFmeanValues = chrSNPs["AF_mean"].tolist()
                        labels[AFmeanIndex], = ax.plot(chrSNPs['Position'],smoothTriangle(AFmeanValues,3),".r",label= "mean frequencies", alpha=0.5)
                        
### plot contig borders
                        
            if options.contigfile:
		# plotting contigs
		if options.contigfile.endswith('.agp'):
			ContigPositions = makeListOfContigPositionsAGP(options.contigfile, Chr)
		elif options.contigfile.endswith('.fna'):
			ContigPositions = makeListOfContigPositions(options.contigfile, Ns , Chr)
		else:
			print 'contigfile must end in .agp or .fna'
			print 'not plotting contigs'
			break
                for contig in ContigPositions:
                    start,end = contig
                    # better use vlines
                    mean_af = calcAFforInterval(chrSNPs, pools[0], Chr, int(start),int(end))
                    if mean_af != 'nan':
                        labels['contigs'], = ax.plot([start,end],[mean_af,mean_af], label="contigs", lw = 2.0, marker = '>')
                    print calcAFforInterval(chrSNPs, pools[0], Chr, int(start),int(end))


### plot genes from gff
            if options.genes_gff:
                print 'plotting genes'
                gff = open(options.genes_gff, 'r')
                listOfGenes = gff2ExonList(gff.readlines())
                plotExonList(listOfGenes, ax)
                gff.close()
                del listOfGenes

### plot undefined gff
            if options.gff:
                    plotgff(options.gff, Chr, ax)


### calculate and plot intervals TODO in funktionen
            if options.interval and len(pools) > 1:

                            print 'starting iterative Interval search'
                            Statisticsfile.write('Interval search for ' + Chr +" \n" )

                            minLengthOfInterval = options.interval
                            for threshold in range(10,30,5):
                                    threshold = threshold /100.
                                    print "Testing threshold: " + str(threshold)
                                    for outlierTolerance in range(1,3):
                                            print "and outlier Tolerance of: " + str(outlierTolerance)
                                            intervals1, filterCount1, filteredIntervals1 = nDI(chrSNPs,  threshold, pools, 0, outlierTolerance, minLengthOfInterval)
                                            intervals2, filterCount2, filteredIntervals2 = nDI(chrSNPs,  threshold, pools, 1, outlierTolerance, minLengthOfInterval)
                                            intervals = intervals1 + intervals2
                                            filterCount = filterCount1 + filterCount2
                                            filteredIntervals = filteredIntervals1 + filteredIntervals2

                                            if len(intervals) >0:
                                                    break
                                    if len(intervals) >0:
                                            break
                            if len(intervals) >0:
                                    print 'Found ' + str(len(intervals)) + ' intervals with threshold ' + str(threshold) + ' and outlier tolerance ' + str(outlierTolerance)
                                    print str(len(intervals1)) + ' intervals for pool ' + pools[0]
                                    print str(len(intervals2)) + ' intervals for pool ' + pools[1]
                                    print 'Length of intervals: ' + str(sum(tup[2] for tup in intervals))
                                    Statisticsfile.write('Found intervals with threshold ' + str(threshold) + ' and outlier tolerance ' + str(outlierTolerance)+ " \n")
                                    Statisticsfile.write(str(len(intervals1)) + ' intervals for pool ' + pools[0])
                                    Statisticsfile.write(str(intervals1)+ " \n")
                                    Statisticsfile.write('Length of intervals: ' + str(sum(tup[2] for tup in intervals1))+ " \n")
                                    Statisticsfile.write('Filtered ' + str(filterCount1) + ' intervals, with length: ' + str(sum(tup[2] for tup in filteredIntervals1))+ " \n")
                                    Statisticsfile.write(str(filteredIntervals1)+ " \n")


                                    if len(intervals1) > 1:
                                        print intervals1
                                        print 'Length of intervals: ' + str(sum(tup[2] for tup in intervals1))
                                        intervals1 = mergeIntervals(intervals1, 0)
                                        print "Smoothed intervals for " + pools[0]
                                        print intervals1
                                        print 'Length of intervals: ' + str(sum(tup[2] for tup in intervals1))
                                        Statisticsfile.write("Smoothed intervals for " + pools[0])
                                        Statisticsfile.write(str(intervals1))
                                        Statisticsfile.write('Length of intervals: ' + str(sum(tup[2] for tup in intervals1)))               
                                    
                                    Statisticsfile.write('Found intervals with threshold ' + str(threshold) + ' and outlier tolerance ' + str(outlierTolerance)+ " \n")
                                    Statisticsfile.write(str(len(intervals2)) + ' intervals for pool ' + pools[1] + "\n")
                                    Statisticsfile.write(str(intervals2)+ " \n")
                                    Statisticsfile.write('Length of intervals: ' + str(sum(tup[2] for tup in intervals2))+ " \n")
                                    Statisticsfile.write('Filtered ' + str(filterCount2) + ' intervals, with length: ' + str(sum(tup[2] for tup in filteredIntervals2))+ " \n")
                                    Statisticsfile.write(str(filteredIntervals2)+ " \n")
                     
                            
                            else:
                                    print 'No intervals could be found'
                                
                            print 'Number of filtered Intervals: ' + str(filterCount1 + filterCount2)


##### assigning color of the interval plots according to the color of the pools
                            if pools[0] in poolcolors:
                                color1 = str(poolcolors[pools[0]])
                            else:
                                color1 = 'r'

                            if pools[1] in poolcolors:
                                color2 = str(poolcolors[pools[1]])
                            else:
                                color2 = 'g'


                            if len(Chrs)> 3:
                                linetype = 'o'
                            else:
                                linetype = '-'
                            if len(intervals1) > 0:
                                labels['intervals ' + pools[0]], = plotIntervals(Chr, intervals1, ax, color1, linetype)
                            #if len(intervals2) > 0:
                                #labels['intervals ' + pools[1]], = plotIntervals(Chr, intervals2, ax, color2, linetype)
                    

#### vergleich von zwei vcf files; es werden die delta values der zwei vcfs geplottet        
        elif options.vcffile1 and options.vcffile2 :

            print 'filtering for Chromosome: ' + Chr
            chrSNPs1 = filterDataFrameForChromosome(df_SNPs,Chr) # only keep SNPs on Chromosome Chr
            chrSNPs2 = filterDataFrameForChromosome(df_SNPs2,Chr) # only keep SNPs on Chromosome Chr
                  
            if len(chrSNPs1)==0 and len(chrSNPs2)==0:
                    continue
            else:
                    chrsToPlot.append(Chr)
                            
                    
            
            #rect_cones = [left, bottom, width, height]
            rect_cones = [0.+1./len(Chrs)*(len(listOfAxes)), 0.1, 1./len(Chrs), 0.9]
            print rect_cones
            listOfAxes.append(plt.axes(rect_cones))
            print listOfAxes
            ax = listOfAxes[-1]
            print ax

            
            if chrSNPs1['Position'].irow(-1) > chrSNPs2['Position'].irow(-1):
                  last_SNPs.append(chrSNPs1['Position'].irow(-1)) # save the length of the Chromosome
            else:
                  last_SNPs.append(chrSNPs2['Position'].irow(-1)) # save the length of the Chromosome
             


            plt.xlabel(Chr) # name the x axis of the subplot after the chromosome
            if len(Chrs) >1 :
                    ax.set_xticklabels([]) # don't plot x axis ticks for any of the plots
            if Chrs.index(Chr) >= 1: # don't plot y axis for subplots
                    ax.set_yticklabels([])

        # plotting delta values TODO::::::
            plot_delta(chrSNPs1, options.windowSize, pools, ax, labels,  "y", "delta frequencies vcf1" ,  "delta vcf1", boost_flag = options.boost)
            plot_delta(chrSNPs2, options.windowSize, pools, ax, labels,  "b", "delta frequencies vcf2" ,  "delta vcf2", boost_flag = options.boost)
                         

                  





## adjusting the subplots according to the length of the Chromosome

start = 0.12 # leave space for the ylabel
plot_end = 0.85 # leave some space to the right
left = start

## iterating through all subplots/axes
for subplots in range(len(chrsToPlot)):

	ax = listOfAxes[subplots] # get one axes instance
	
	# set x and y limits, the plot can not be longer than the last found snp
	ax.axis('tight')
	ax.set_ylim([-0.2, 1.2])
	xmin, xmax, ymin, ymax = ax.axis()
	
	
	# modify the ticks and label for the x-axis
	start, end = ax.get_xlim()
	##### set the size of the ticks
	ax.tick_params(axis='y', labelsize = 13)
	ax.tick_params(axis='x', labelsize = 13)

	if len(chrsToPlot) < 3:
		import matplotlib.ticker as tkr     # has classes for tick-locating and -formatting
		yfmt = tkr.FuncFormatter(numfmt)    # create your custom formatter function


		#"Using standard ticks"
		ax.xaxis.set_ticks(numpy.arange(0, end, 1000000.0))
		ax.ticklabel_format(axis = 'x', style = 'plain')
		ax.xaxis.set_major_formatter(yfmt)

		ax.set_xlabel(chrsToPlot[subplots] + ' [Mb]', fontsize = 15)
		
	else:
                #"adjusting x-ticks and labels"
		ax.xaxis.set_ticks(numpy.arange(0, end, 10000000))
			# rotate the label of the axes
		ax.set_xlabel(chrsToPlot[subplots] ,rotation=45, fontsize = 12)
                #ax.set_xlabel(chrsToPlot[subplots] + ' [10 Mb]',rotation=45, fontsize = 12)
	
	
	# Abmessungen des subplots bestimmen, 
	width = float((last_SNPs[subplots]-xmin))/numpy.sum([a - b for a, b in zip(last_SNPs,first_SNPs)]) * plot_end # breite ist Anteil dieses Chromosomes an der Gesamtlaenge mal 0.9
	
	rect_cones = [left, 0.1, width, 0.8]
	ax.set_position(rect_cones)
	if ax2:
		ax2.set_position(rect_cones)
	left = left + width  # der linke Rand des neuen Chromosomenplots ist das ende des vorherigen
	
	
###### plot the allelefrequency histogram, if only one chromosome is plotted
	if options.histogram and len(chrsToPlot) == 1:
### arranging the plot layout, ax is the selected plot
        # rearrange the SNP plot to half of its former width, so a second plot can be shown
		ax.set_position([left-width+0.02, 0.1, width/2, 0.8])
        # calculate parameters for the histogram plot
		width = width/2
		bottom, height = 0.1, 0.8
		bottom_h = left_h = left+width+0.02
		rect_histy = [left -width + 0.06, bottom, width - 0.02, height]
### create and select a new plot in which to plot the histogram
		listOfAxes.append(plt.axes(rect_histy)) # create
		hist_ax = listOfAxes[-1] # select

### plot the histogram
		if len(pools) == 2:
			chrSNPs[['Frequency'+pools[0],'Frequency'+pools[1]]].plot(kind='hist', orientation='horizontal', bins=100, alpha=0.5, logx=True, ax = hist_ax, ylim=(-0.2, 1.2))
		elif len(pools) == 1:
			chrSNPs['Frequency'+pools[0]].plot(kind='hist', orientation='horizontal', bins=100, logx=True, ax = hist_ax, ylim=(-0.2, 1.2))
		else:
			print "more than two pools, plotting only first two pools"
			chrSNPs[['Frequency'+pools[0],'Frequency'+pools[1]]].plot(kind='hist', orientation='horizontal', bins=100, alpha=0.5, logx=True, ax = hist_ax, ylim=(-0.2, 1.2))

# setting y and y label parameters for histogram.
		hist_ax.tick_params(axis='y', labelleft='off', labelsize = 13)
		hist_ax.set_xlabel('Frequency', fontsize=15)
		hist_ax.tick_params(axis='x', labelsize = 13)


####### labels
labels_list = []
labels_names = []
for key,value in labels.items():
        labels_list.append(key)
        labels_names.append(value)

#print labels_list
#['intervals P2', 'intervals P1', 'delta']
#print labels_names
#[<matplotlib.lines.Line2D object at 0x7fd09dd06f50>, <matplotlib.lines.Line2D object at 0x7fd09dd06210>, <matplotlib.lines.Line2D object at 0x7fd09dec6510>]


###setting xticks and xlabels for all plots

if len(labels_list)>1 or options.histogram: # if there is more than one plot, rotate all x tick labels
	for ax in listOfAxes:
		plt.setp(ax.get_xticklabels(), rotation=45)

if len(labels_list)>3: # if there are more than 3 plots rotate all x labels (Chromosome names)
	for ax in listOfAxes:
		plt.setp(ax.xaxis.get_label(), rotation=45)






##### set a y-label only for the leftmost plot
ax = listOfAxes[0]
if options.delta:
    ax.set_ylabel('delta allele frequency estimate', fontsize=15)
else:
    ax.set_ylabel('allele frequency estimate', fontsize=15)



        
#plt.legend(labels_names, labels_list, loc ='best', ncol = 1, shadow = False, numpoints = 1)
fig.legend(labels_names, labels_list, loc = "upper center", ncol = 6, shadow = False, numpoints = 1, frameon = False)



### clear memory
del chrSNPs
del df_SNPs



### show the plot
if options.batch == False:
    plt.show()

###save plot
fig.set_size_inches(7,5)
fig.savefig(Path2OutputFile, dpi=300)
#save_plot(fig, Path2OutputFile, dpi=300) # speichert den plot so, wie man ihn sieht


Statisticsfile.close()



