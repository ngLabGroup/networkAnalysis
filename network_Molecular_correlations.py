# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 09:34:49 2018

@author: Trevor
"""
#****************************************************************
#This script computes the correlations between network data and chemistry data
#****************************************************************
#clearn Python import get_ipython
get_ipython().magic('reset -sf')

import networkx as nx
import os
import time
from networkx.algorithms import community
import pandas as pd
import numpy as np
from sklearn.preprocessing import normalize
import rdkit as rd
from rdkit import Chem
from rdkit.Chem import Descriptors
from scipy.stats import pearsonr
from scipy.stats import normaltest
from scipy.stats import spearmanr
from scipy.stats import kendalltau
from scipy.stats import levene
#****************************************************************

#****************************************************************
##Load up another file with other data to add. 

dfQF = pd.read_excel(r"C:\ResearchWorkingDirectory\gitRepositories\ExampleNT_BTW_Nodes.xlsx", sheet_name = 'Sheet1')

epiLogP = pd.read_excel(r"C:\ResearchWorkingDirectory\gitRepositories\ExampleEpiLogP.xlsx", sheet_name = 'Sheet1')

epiLogP.drop_duplicates(inplace = True)
epiLogP.reset_index(inplace = True, drop = True)

smilesList = list(dfQF['SMILES'])
#drop all the extra stuff out of this:
e2 = epiLogP[epiLogP['SMILES'].isin(smilesList)]
e2.reset_index(drop = True, inplace = True)

dfQF = dfQF[dfQF['SMILES'].isin(list(e2['SMILES'])) ]

#Could use this function to add a lot of different stuff
for i, row in dfQF.iterrows():
    
#    dfQF.at[i, 'Betweeness Centrality'] = btw[row['SMILES']]    
#    dfQF.at[i, 'Authority Betweeness'] = a[row['SMILES']]
#    dfQF.at[i, 'Hub Betweeness'] = h[row['SMILES']]
        
    #Chemisry parameters
    m = Chem.MolFromSmiles(row['SMILES'])   
    #general structural components
    dfQF.at[i, 'Fraction CSP3' ] = Chem.Lipinski.FractionCSP3(m)
    dfQF.at[i, 'Molecular Weight'] = rd.Chem.Descriptors.ExactMolWt(m) 
    dfQF.at[i, 'Number of Rings' ] = Chem.Lipinski.NumAromaticRings(m)
 
    #non ring carbons
    patt= Chem.MolFromSmarts("[C]")
    dfQF.at[i, 'Non-Ring Carbons' ] = len(m.GetSubstructMatches(patt))
    
#    #All Carbons
    pat2= Chem.MolFromSmarts("[#6]")
    dfQF.at[i, 'All Carbons' ] = len(m.GetSubstructMatches(pat2)) 
    
#    #H acceptors and Donors
    dfQF.at[i, 'Number of H Acceptors' ] = Chem.Lipinski.NumHAcceptors(m)           
    dfQF.at[i, 'Number of H Donors' ] = Chem.Lipinski.NumHDonors(m)  

#    #Hydroxyl Groups
    dfQF.at[i, 'Alphatic Hydroxyl Groups' ] = rd.Chem.Fragments.fr_Al_OH(m)  
    dfQF.at[i, 'Aromatic Hydroxyl Groups' ] = rd.Chem.Fragments.fr_Ar_OH(m)
    dfQF.at[i, 'Total Hydroxyl Groups' ] = rd.Chem.Fragments.fr_Al_OH(m) + rd.Chem.Fragments.fr_Ar_OH(m)
    
    #Carbonxylic Acids
    dfQF.at[i, 'Carbonyl Groups' ] = rd.Chem.Fragments.fr_C_O(m)  
    dfQF.at[i, 'Carboxylic Acid Groups' ] = rd.Chem.Fragments.fr_COO(m)
    dfQF.at[i, 'Carbonyl O, excluding COOH' ] = rd.Chem.Fragments.fr_C_O_noCOO(m)
    
    #Additional Parameters
    dfQF.at[i, 'Heavy Atom Count' ] = rd.Chem.Lipinski.HeavyAtomCount(m)    
    dfQF.at[i, 'Number Aliphatic Carbocycles' ] = rd.Chem.Lipinski.NumAliphaticCarbocycles(m)    
    dfQF.at[i, 'Number Aliphatic Heterocycles' ] = rd.Chem.Lipinski.NumAliphaticHeterocycles(m)    
    dfQF.at[i, 'Number Aliphatic Rings' ] = rd.Chem.Lipinski.NumAliphaticRings(m)    
    dfQF.at[i, 'Number Aromatic Carbocycles' ] = rd.Chem.Lipinski.NumAromaticCarbocycles(m)   
    dfQF.at[i, 'Number Hetero atoms' ] = rd.Chem.Lipinski.NumHeteroatoms(m)  
    dfQF.at[i, 'Number Rotatable Bonds' ] = rd.Chem.Lipinski.NumRotatableBonds(m)
    dfQF.at[i, 'Number Saturated Heterocycles' ] = rd.Chem.Lipinski.NumSaturatedHeterocycles(m)
    dfQF.at[i, 'Number Saturated Rings' ] = rd.Chem.Lipinski.NumSaturatedRings(m)   
    dfQF.at[i, 'Wildman-Crippen LogP' ] = rd.Chem.Crippen.MolLogP(m)
    smilesString = row['SMILES']
    dfQF.at[i, 'Oxygens' ] = smilesString.count('O')
    dfQF.at[i, 'EpiLogP'] = e2.loc[e2['SMILES'] == row['SMILES'], 'PredEpiLogP'].iloc[0]
    
#Correlation Matrix
#****************************
#create the cross correlation matrix
toCor = list(dfQF)
toCor.remove('SMILES') 
toCor2 = toCor

pearsonDF = pd.DataFrame([], columns = toCor)
kendalltauDF = pd.DataFrame([], columns = toCor)

spearmanDF = pd.DataFrame([], columns = toCor) 
spearmanDFnorm = pd.DataFrame([], columns = toCor) 

iterator = 0

for t in toCor:
    iterator +=1
        
#    correlationDF.loc[t] = [0] * len(toCor)         
    #All Possible Parameters
    
    #print('The toCor variable is ', toCor)
    for t2 in toCor2:
        
        #correlate i with all the other members in the dataset
        y = dfQF[t]
        a = normalize(y[:,np.newaxis], axis=0).ravel()
        
        x = dfQF[t2]
        b = normalize(x[:,np.newaxis], axis=0).ravel()        
        
        [s,p] = normaltest(y)
        if p > 0.01:
            print('The following group is not normalally distributed: ', t)
               
        #pearson
        pearsonDF.loc[t,t2] = np.corrcoef(y,x)[0,1]
        
        #spearmanr      
        [cor, p] = spearmanr(y,x)      
        spearmanDF.loc[t,t2] = cor        
        
        [cor, p] = spearmanr(a,b)      
        spearmanDFnorm.loc[t,t2] = cor
        
        #kendalltau correlations
        [cor, p] = kendalltau(y,x)
        kendalltauDF.loc[t,t2] = cor
        

#Write out two different data frames, one with edgedata, one with node data. 

#########Write out the updated excel file
#os.chdir(r'C:\ResearchWorkingDirectory\DataReferenceFiles\PaperData\Spearmans')
##
#timestr = 'SpearmanFluoreneRaw' + time.strftime("%Y%m%d-%H%M%S") +'.xlsx'
#print (timestr)
#writer = pd.ExcelWriter(timestr)
#spearmanDF.to_excel(writer,'Sheet1')
#writer.save()


#os.chdir('C:\ResearchWorkingDirectory\DataReferenceFiles\PaperData')
##
#timestr = 'DFqfAcen' + time.strftime("%Y%m%d-%H%M%S") +'.xlsx'
#print (timestr)
#writer = pd.ExcelWriter(timestr)
#dfQF.to_excel(writer,'Sheet1')
#writer.save()


   
#code for developing all the correlations
    
#y = df['BetweennessCentrality']
##norm1 = x / np.linalg.norm(x)
#a = normalize(y[:,np.newaxis], axis=0).ravel()
#
#x = df['FromWeight']
##norm1 = x / np.linalg.norm(x)
#b = normalize(x[:,np.newaxis], axis=0).ravel()
#
#randBase = np.random.rand(1,22996)
#
#
#np.corrcoef(a,b)

