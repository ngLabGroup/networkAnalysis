# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 14:22:08 2019

@author: Trevor
"""

#******************************************************************

#This script allies the Relevant Pathway filter and writes out the ExamplePathWeight.xlsx

#******************************************************************
get_ipython().magic('reset -sf')

import os
import networkx as nx
import pandas as pd
import numpy as np

os.chdir(r'C:\ResearchWorkingDirectory\gitRepositories\networksPlotting')

from CustomFunctions_PAHteam import FindSourceSink
from CustomFunctions_PAHteam import getIncomingNodes

#need the logP values from the EpiSuiteLogP script
LogP = pd.read_excel(r"C:\ResearchWorkingDirectory\gitRepositories\ExampleEpiLogP.xlsx", sheet_name = 'Sheet1')

#read in the raw DF so that we can build the network
df = pd.read_excel(r"C:\ResearchWorkingDirectory\gitRepositories\RawExampleData.xlsx", sheet_name = 'Sheet1')

df.drop_duplicates(inplace = True)
df.reset_index(inplace = True, drop = True)

G_raw = nx.from_pandas_edgelist(df, source ='From', target = 'To',edge_attr=True, create_using=nx.DiGraph())
    
[SourceN, SinkN] = FindSourceSink(G_raw)

#Select rows based on logP <2. This could be set to whatever you wish
HighLogP = LogP.loc[LogP['PredEpiLogP'] >= 2.0]

allRelNodes = []

for i, row in HighLogP.iterrows():
    s = row['SMILES']
    releventNodes = getIncomingNodes(G_raw, s, SourceN)
    allRelNodes = list(set(releventNodes+allRelNodes))
    
extraSinks = []
pathCompounds = []    
newDF = df.iloc[0:0]
    
#now I have all the high throughput nodes, just need to get the right edges
for a in allRelNodes:
    #first get all the edges that contain a relevent compound of interest
    tempa = df.loc[(df['From'] == a) | (df['To'] == a)]
    newDF = newDF.append(tempa)

    outEdges = list(G_raw.out_edges(a))

    #get all the edges that go To each of the high interest compounds
    for o in outEdges:                 
        tempc = df.loc[(df['From'] == o[0]) & (df['To'] == o[1])]
        #also add the sink node to a reference list for the color scheme below
        extraSinks.append(o[1])
        newDF = newDF.append(tempc)

newDF = newDF.drop_duplicates()
extraSinks = list(set(extraSinks))
pathCompounds = list(set(pathCompounds))

#write off the final result
os.chdir('C:\ResearchWorkingDirectory\gitRepositories')

newDF.to_excel('ExamplePathWeight.xlsx', sheet_name='Sheet1', engine='xlsxwriter' )


