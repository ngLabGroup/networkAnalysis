# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 16:40:28 2019

@author: twsle
"""
get_ipython().magic('reset -sf')

import os
import networkx as nx
import pandas as pd
import numpy as np
import rdkit as rd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdFMCS
from rdkit.Chem.Draw import MolToImage
from rdkit.Chem import AllChem
from rdkit import DataStructs

#*******************************
#This function accepts a list of SMILES codes and categorizes based on criteria defined 

def FindSubstructures(testMol, benRingMol, smilesList, df):
#smilesList is the list of smiles codes to substructure matched
#numCarbons is how many different carbons atoms may be different for the structure to be considered similiar
   
    fp1 = AllChem.GetMorganFingerprint(testMol,2,useFeatures=True) 
    
    #first look for a perfect match and return that if available. 
    if testMol in smilesList:
        return 'perfect match'
    
    pat1= rd.Chem.MolFromSmarts("[#6]")
    testCarbons = len(testMol.GetSubstructMatches(pat1)) 
    testRings = rd.Chem.rdMolDescriptors.CalcNumAromaticRings(testMol)    
            
    #first, find all the structures that have the pre-defined backbone structure. This goes in the "ringsMatchDF"
    ringMatchDF = df.loc[df['Rings'] == 1]
    
    #carbs Match gets molecules that have the pre-defined range for a partial match plus 1 (we will narrow this down later)
    carbsMatchDF = ringMatchDF.loc[(ringMatchDF['Carbons'] >= (testCarbons -1)) & (ringMatchDF['Carbons'] <= (testCarbons +1) ) ]
    
    #finally we iterate through the carbsMatchDF and check the maximum common structure and look for the biggest ones
    
    for i,row in carbsMatchDF.iterrows():
        currentSmiles = row['SMILES']
        
        m2 = Chem.MolFromSmiles(row['SMILES'])
        mols = [m2, testMol ]
        res = Chem.rdFMCS.FindMCS(mols,completeRingsOnly=True)
        
        
        ss = rd.Chem.MolFromSmarts(res.smartsString)   
        testS = res.smartsString
        
    
        fp2 = AllChem.GetMorganFingerprint(m2,2,useFeatures=True) 
        carbsMatchDF.at[i,'Carbons'] = testS.count('[#6]') #carbons
        carbsMatchDF.at[i,'Oxygens'] = testS.count('[#8]') #oxygens
        carbsMatchDF.at[i, 'Bonds'] = res.numBonds
        carbsMatchDF.at[i, 'smarts'] = testS
        carbsMatchDF.at[i, 'Dice Sim'] = DataStructs.DiceSimilarity(fp1,fp2)
         
    #sort the dataframe by numCarbMatch, then numOxygenMatch, then
            
#    carbsMatchDF.sort_values(by ='Bonds', ascending = False, inplace = True)
#    carbsMatchDF.sort_values(by ='Oxygens', ascending = False, inplace = True)
    carbsMatchDF.sort_values(by ='Carbons', ascending = False, inplace = True)
#    carbsMatchDF.sort_values(by ='Dice Sim', ascending = False, inplace = True)        

    
    #*****************************
    return carbsMatchDF

#*******************************


#read in a list of smiles
    
dfQF = pd.read_excel(r"C:\ResearchWorkingDirectory\DataReferenceFiles\PaperData\DFqfPhen20190613-224810.xlsx", sheet_name = 'Sheet1')

dfQF.sort_values(by ='Betweeness Centrality', ascending = False, inplace = True)
dfQF.reset_index(inplace = True, drop = True)

tempSmiles = dfQF[['SMILES', 'Betweeness Centrality']]  #.head(1000)
smilesList = list(tempSmiles['SMILES'])


#The user must define a minimum sub-structure for to start from. Usually one or two benzene rings is used. 

#napthalene as a base structure
#benRingMol = Chem.MolFromSmiles('C1C2=C(C=CC=C2)C2=C1C=CC=C2')

#benzene as a base structure
benRingMol = Chem.MolFromSmiles('C1=CC=CC=C1')


#************************************************************************
df = pd.DataFrame(columns=['SMILES','nodeTransferWeights110100','Rings','Carbons'])
for s in smilesList:
    m = rd.Chem.MolFromSmiles(s)

    pat2= rd.Chem.MolFromSmarts("[#6]")
    numCarbs = len(m.GetSubstructMatches(pat2)) 
    
    currentRings = len(m.GetSubstructMatches(benRingMol) )
    
    #supposed to only append if it gets the right base substructure
    #something's wrong here
    ntw = float(dfQF.loc[dfQF['SMILES'] == s, 'nodeTransferWeights110100'])
    
    #only include ones that at least match the substrate:
    if currentRings > 0:
        dfNext = pd.DataFrame(data ={'SMILES': [s], 'nodeTransferWeights110100':[ntw],'Rings':[currentRings], 'Carbons':[numCarbs]})
        df = df.append(dfNext, ignore_index = True)
#************************************************************************
  
    
#run these as individual lines. If you need to change the benRingMol you have to run smilesList loop again
    
testMol = Chem.MolFromSmiles('O=C(O)C=Cc1ccccc1C(=O)O') #,sanitize=False)

output = FindSubstructures(testMol, benRingMol,smilesList,df )
#Sort the output dataframe by Dice Sim. This usually gives the best match. 
#******************************************************


#run the closest match here to get its throughput value
dfQF.loc[dfQF['SMILES'] == 'Oc1cccc(c1O)-c1ccccc1C([O-])=O', 'nodeTransferWeights110100']











