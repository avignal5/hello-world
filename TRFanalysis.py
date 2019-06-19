#!/usr/bin/env python3
# _*_ coding: Utf-8 _*_
# coding: utf-8

#CHANGED?

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
#from matplotlib import figure as figure
import re
import csv
from collections import defaultdict

#Import Tandem Repeat Finder (TRF) data
dataDict = defaultdict(list)
with open('/Users/avignal/Documents/Stats/2016_PacificBee/TandemRepeatFinder/GCA_003314205.1_INRA_AMelMel_1.0_genomic.fna.2.7.7.80.10.50.1000.dat') as csvFile:
	data=csv.reader(csvFile, delimiter = " ")
	for row in data:
		if row:
			if re.search('Sequence',row[0]):
				chromosome = row[1]
			elif re.search('\d',row[0]): # and re.match('CM', chromosome) :
				dataDict["Chr"].append(chromosome)
				dataDict["Start"].append(int(row[0]))
				dataDict["End"].append(int(row[1]))
				dataDict["Period_Size"].append(int(row[2]))
				dataDict["Copy_Number"].append(float(row[3]))
				dataDict["ConsensusSize"].append(int(row[4]))
				dataDict["Percent_Matches"].append(int(row[5]))
				dataDict["Percent_Indels"].append(int(row[6]))
				dataDict["Score"].append(int(row[7]))				
				dataDict["A"].append(int(row[8]))
				dataDict["C"].append(int(row[9]))
				dataDict["G"].append(int(row[10]))
				dataDict["T"].append(int(row[11]))
				dataDict["Entropy"].append(float(row[12]))
				dataDict["Repseq"].append(row[13])
df = pd.DataFrame.from_dict(dataDict)
df['StartMb'] = df['Start'] / 1000000
df['EndMb'] = df['End'] / 1000000

'''
Edit for repeat name (title and figure name), chromosomes or unknown, for selection
'''

#Figure title
title = "10-1000 bp repeats"

#Repeat name
repeat = "10_1000"

#Select data for the 16 chromosomes, or all data, or only the unknown
df1 = df[df.Chr.str.contains('^CM')]					#Chromosomes only
#df1 = df[df.Chr.str.contains('^QFDB')]					#Unknown only
#df1 = df												#All data

#Select specific repeat types
#df2 = df1												#No selection
#df2 = df1[(df1.Copy_Number > 50) & (df1.Period_Size > 80) & (df1.Period_Size < 110)]
df2 = df1[(df1.Period_Size > 10)]						#Specific search

#Text
t = "Period_Size = 10-1000"

'''
End edit
'''

#List of chromosomes
chromosomes = df1['Chr'].unique().tolist()
chromosomes.reverse()
ybins = len(chromosomes)

chrNames = list(range(1,ybins+1))
chrNames = ["Chr " + str(s) for s in chrNames]
chrNames.reverse()

#Figure settings
space = (1 / (ybins + 1)) * 1 / 6
figHeight = (ybins + 1) * 3
figWidth = (df['EndMb'].max()) * 1.4 # + space * 4

fig = plt.figure(figsize=(figWidth, figHeight))
#plt.axes([0, 10, 0, 10])
fig.suptitle(title, fontsize=96)
fig.text(.5,.2,t, wrap=True, fontsize=24)
bottom = space

#Plot chromosomes
count=0
for chr_i in chromosomes:
    chrom = df.Chr == chr_i
    df_chrom = df[chrom]
    chrom2 = df2.Chr == chr_i
    df2_chrom = df2[chrom2]    
    
    #Define size and draw box
    left = space * 4
    width = (df_chrom['EndMb'].max()) * 1.4 /figWidth * 0.90
    height = (1/ybins) * 4.5/6
    ax = plt.axes([left, bottom, width, height])
    bottom = bottom + height + space
    
    #plot
    ax.scatter(x='StartMb', y='Copy_Number', s=200, color='red', data=df2_chrom) #[chrom])
    ax.set_xlim(0, df_chrom['StartMb'].max())
    ax.set_ylim(0, df2['Copy_Number'].max())
    #ax.set_title(chrNames[count])
    
    #Labels, ticks, grids
    ax.set_ylabel(chrNames[count], fontsize=16)
    major_ticks = np.arange(0,df_chrom['EndMb'].max(),1)
    ax.set_xticks(major_ticks)
    minor_ticks = np.arange(0,df_chrom['EndMb'].max(),0.250000)
    ax.set_xticks(minor_ticks, minor=True)
    ax.grid(which='minor', alpha=0.2)
    ax.grid(which='major', alpha=0.5)
    count = count + 1
plt.savefig('/Users/avignal/Documents/Stats/2016_PacificBee/TandemRepeatFinder/chrsRepeat' + repeat + '.png', format="png")#, dpi=600)
#plt.savefig('/Users/avignal/Documents/Stats/2016_PacificBee/TandemRepeatFinder/chrs.pdf', format="pdf")

#plt.show()

