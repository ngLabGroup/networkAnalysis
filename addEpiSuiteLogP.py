# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 08:58:50 2019

@author: Trevor
"""

#clearn Python import get_ipython
get_ipython().magic('reset -sf')

import os
import networkx as nx
import pandas as pd
import numpy as np
import re
import time

#This reads and writes for Episuite. EpiSuite crashes if you give it too large a dataset, so we send only 10000 compounds at a time and then use this script to build them back into an excel file that can be read in by the network correlations script. 

#read in one of the reference files
df = pd.read_excel(r"C:\ResearchWorkingDirectory\gitRepositories\RawExampleData.xlsx", sheet_name = 'Sheet1')

#set the working directory
os.chdir(r'C:\ResearchWorkingDirectory\gitRepositories')

#EpiSuite can only handle about 10000 inputs at a time. 
#being lazy and not writing a proper loop since this will have to be customized for each PAH network
smilesCodes = list(set(df['From'])  |   set(df['To']) )
smilesCodes1 = smilesCodes[0:10010]
smilesCodes2 = smilesCodes[10000:20010]
smilesCodes3 = smilesCodes[20000:]
#smilesCodes4 = smilesCodes[30000:40010]
#smilesCodes5 = smilesCodes[40000:50010]
#
#smilesCodes6 = smilesCodes[50000:]
#smilesCodes7 = smilesCodes[60000:70010]
#smilesCodes8 = smilesCodes[70000:80010]
#smilesCodes9 = smilesCodes[80000:90010]
#smilesCodes10 = smilesCodes[100000:110010]
#
#smilesCodes11 = smilesCodes[110000:120010]
#smilesCodes12 = smilesCodes[120000:130010]
#smilesCodes13 = smilesCodes[130000:140010]
#smilesCodes14 = smilesCodes[140000:150010]
#smilesCodes15 = smilesCodes[150000:160010]
#
#smilesCodes16 = smilesCodes[160000:170010]
#smilesCodes17 = smilesCodes[170000:180010]
#smilesCodes18 = smilesCodes[180000:190010]
#smilesCodes19 = smilesCodes[190000:200010]
#smilesCodes20 = smilesCodes[200000:210010]
#
#smilesCodes21 = smilesCodes[210000:220010]
#smilesCodes22 = smilesCodes[220000:230010]
#smilesCodes23 = smilesCodes[230000:240010]
#smilesCodes24 = smilesCodes[240000:250010]
#smilesCodes25 = smilesCodes[250000:]
#
#smilesCodes26 = smilesCodes[90000:100010]




#**********************************************
timestr = 'EpiSuiteBatchSmiles'  + time.strftime("%Y%m%d-%H%M%S") +'.txt'

#This lets you write off a series of files to send through EpiSuite
with open(timestr, "w") as output:
    for s in smilesCodes1:
        output.write(str(s)+ '\n') 
output.close()
#*****************************************
######################################
#now run all this through episuite



#once you have the EpiSuite Batch files you can run the lower part
os.chdir('C:\ResearchWorkingDirectory\gitRepositories')
#####################################
#initialize the output dataframe
outData = pd.DataFrame(columns = ['SMILES', 'PredEpiLogP', 'EpirEpiLogP'])

for f in list(range(1,2)):
    #*********************************************
    #Load up the batch mode stuff
    
    filename = str(f).join(['BATCH0' ,'.out'])   
    
    file_object  = open(filename, 'r') 
    
    #This gives you one long text files
    inputData = file_object.read()
    file_object.close()
    inputData.find('(est)')
    
    #********************************************
    for m in re.finditer('(est)', inputData):
        #lets get the calculated logP values        
        logP = float(inputData[(m.start()-6): (m.start()-1) ])
          
        #m.start*() start of the current row: search FROM HERE to the end of the file and fine the next '\n' This gives us the end of the smiles string
        sIndex = inputData[m.start():].find('\n')-1
    #    print(sIndex)
    #    print(m.start())
    
        sVal = inputData[(m.start()+18): (m.start()+sIndex) ]
    
        #check if there are any experimental values
        expVal = inputData[m.start()+7: (m.start()+11) ]
        if expVal != '----':
            #print(inputData[m.start()+5: (m.start()+13) ] )
#            print ('Experimental for', sVal, 'is ',expVal)
#            print ('Predicted for', sVal, 'is ',logP)
            outData = outData.append({'SMILES':sVal, 'PredEpiLogP':logP, 'EpirEpiLogP': float(expVal)}, ignore_index = True)
        else:            
            outData = outData.append({'SMILES':sVal, 'PredEpiLogP':logP}, ignore_index = True)
    #********************************************

#write out the updated data files. 
outData.drop_duplicates(inplace = True)

####Write out the updated excel file
os.chdir('C:\ResearchWorkingDirectory\gitRepositories')

timestr = 'ExampleEpiLogP.xlsx'
print (timestr)
writer = pd.ExcelWriter(timestr)
outData.to_excel(writer,'Sheet1')
writer.save()


#Rescue Feature. This is for fixing the headings if you mess up the file. Typically commented out. 
#**********************************************************
#backup = outData
#outData = pd.read_excel(r"C:\ResearchWorkingDirectory\KOWWIN LogP\FluoreneEpiLogP.xlsx", sheet_name = 'Sheet1')
#
#for i, row in outData.iterrows():
#    toFix = row['SMILES']
#    toFix = toFix.replace(' ','')
#    outData.at[i, 'SMILES'] = toFix
