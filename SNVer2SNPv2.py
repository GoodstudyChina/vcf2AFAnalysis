#!/usr/bin/python

# script used to convert the output of SNVer, which is in its own vcf like format, to a tab delimited txt format, which  is easier to handle

import sys
import os
import re
import csv
from string import translate,maketrans,punctuation


def calcAF(altBases,allBases): # calculates the allele frequency
    refBases = int(allBases)-int(altBases) # number of reads concordant to the reference

    if refBases == 0:
        allelefrequency = 1
    elif altBases == 0:
        allelefrequency = 0
    else:            
        allelefrequency = float(altBases)/int(allBases)

#    print 'Calc Freq: '
#    print "Alt: " + str(altBases)
#    print "Alt: " + str(float(altBases))
#    print "All: " + str(allBases)
#    print "All: " + str(float(allBases))
    
#    print allelefrequency
    return allelefrequency

def GATKprocessPoolField(poolField):
    # Takes one Poolfield; splits it at the ':' delimiter into in a list of five fields, takes the second field, which is a tuple of two, and gives its two values to the calcAF Funktion
    # define Transitionmatrix, to exchange delimiters by whitespace
    delims = ':'
    T = maketrans(delims, ' '*len(delims))
    poolField = translate(poolField.strip(), T).split()
    if len(poolField) == 5:
	    Bases = poolField[1].split(',')
	    refBases  = int(Bases[0])
	    altBases =  int(Bases[1])
	    allBases = refBases + altBases
	    allelefrequency = calcAF(altBases, allBases)
	    return altBases, allBases, allelefrequency
    else:
	    return 0



def parseVCF(mode, infile,  out):
	

	
	if os.path.isfile(out):
		print "converted variationfile already exists: "
		print out
		return out
	
	filehandle = open(infile, 'r')
	outfile = open(out, 'w+')
	## beginn parsing parse SNVer vcf file
	
	if mode == 'SNVer':
	## parse SNVer vcf file

    
	    for line in filehandle:
        	if 'NA' not in line :
			line = translate(line.strip(), T).split()
           	if not line[0].startswith("#"): # skip header and comments
                    #print line
                	formattedline = line[0:2]
                   	formattedline += line[3:5]
                            
#                    likelihood_allelefrequency = re.search(r'AF=(\d\.*\d*)', line[7])
#                    formattedline.append(likelihood_allelefrequency.group(1))
		  	formattedline.append(0)                    
                    
                   	remainingFields = len(line[9:])
                            
                   	i = 10
                    
                    	while i < len(line):
                          
                            formattedline.append(line[i])
                            formattedline.append(line[i+1])        
                            formattedline.append(calcAF(line[i],line[i+1]))
                            i = i+2

                    	#print('\t'.join(map(str,formattedline)))
			outfile.writelines('\t'.join(map(str,formattedline))+'\n')


## parse GATK vcf file
	elif mode == 'GATK':
		for line in filehandle:
			line = line.strip().split()
			
			if not line[0].startswith('#') and './.' not in line: # skip header and comments and NA values
				formattedline = line[0:2]
				formattedline += line[3:5]
						
			#            likelihood_allelefrequency = re.search(r';AF=(\d\.*\d*)', line[7])
			#            formattedline.append(likelihood_allelefrequency.group(1))
				formattedline.append(0)
				i = 9
				
				printingFlag = True
				
				while i< len(line):
					if ':.:' not in line[i]:
						
						basevalues = GATKprocessPoolField(line[i])
						
						if basevalues != 0:
							altBases, allBases, allelefrequency = basevalues
							formattedline.append(altBases)
							formattedline.append(allBases)
							formattedline.append(allelefrequency)
							i += 1
						else:
							printingFlag = False
							i += 1
					else:
						printingFlag = False
						i += 1
						sys.stderr.write("skipped ") 
						sys.stderr.write('\t'.join(map(str,line))) 
						sys.stderr.write('\n') 
			
				if printingFlag :
					#print('\t'.join(map(str,formattedline)))
					outfile.writelines('\t'.join(map(str,formattedline))+'\n')
		
					
			else:
				sys.stderr.write("skipped ") 
				sys.stderr.write('\t'.join(map(str,line))) 
				sys.stderr.write('\n') 
				
	
	
	filehandle.close()
	outfile.close()
	print "SNP file written" 
        return out    
                



##
##with open('some.csv', 'wb') as f:
##    writer = csv.writer(f)
##    writer.writerows(someiterable)

##
##    if not line[0].startswith("#"): # skip header and comments
##        formatfield = line[7].split(';')
##        if formatfield[0] == 'INDEL':
##            del formatfield[0]
###        print formatfield
##        formattedline = line[0:2]
##        formattedline += line[3:5]
##        allelefrequency = re.search(r'AF1=(\d\.*\d*)', formatfield[2])
##        formattedline.append(allelefrequency.group(1))
##        print formattedline
##

