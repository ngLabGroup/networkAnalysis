# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 11:54:29 2018

@author: Trevor
"""
#******************************************************************
#This script applies the nodeThroughput algorithm and write out the ExampleQFNodes.xlsx file
#******************************************************************
#clearn Python import get_ipython
get_ipython().magic('reset -sf')

import time
import os
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
os.chdir(r'C:\ResearchWorkingDirectory\gitRepositories\networksPlotting')

#requires the CustomFunctions_PAHteam file to work. 
from CustomFunctions_PAHteam import FindSourceSink
from CustomFunctions_PAHteam  import hierarchy_posMS

plt.interactive(False)

def assignFlowWeights(G, node):
#Assign Flow Weights accepts a network and a specific Node and returns the same network with the node's Transfer Weight and it's outgoing edge Transfer weights added as attributes, as well as a flag indicating whether the specified node was a source, a sink, or a middle node
    
#Condition 1. Specificed Node is a Source Node
    #a.)If the node as an attribute Transfer Weight, then use this value and assign the outgoing edges as such. This allows multi-source networks to have different intput values on different nodes
    #b.) assume 1 as the weight of the source node if unspecified
    #Return NodeStatus Flag = 1 in either case
    
#Condition 2. Specified Node is a Middle node 
    #a. The Transfer Weights of all incoming edges have been assigned. 
    # Assign the Transfer Weights of the outgoing edges and the specified node. Assign NodeStatus flag = 1 
    #b. some or all Transfer Weights are Missing. This node will have to be re-processed once the rest of the input edges are found. Do not assign anything and return NodeStatus Flag = 0. 
    
#Condition 3. Specified Node is a Terminal Sink (no out degrees)
    #a. The Transfer Weights of all incoming edges have been assigned. 
    # Assign the node's weight and return NodeStatus Flag = 2. 
    #b. The Transfer Weights of all incoming edges have not been assigned, do not assign anything and return NodeStatus flag = 0
    
#*************************************************    
    import networkx as nx
    #print('I am on node ', node)
    inDeg = list(G.in_edges(node))
    #accept a network and a node to start on and return the network with the weights added to the edges. Assume weight of 1 if the node assigned is a source node
    if (inDeg == []):
    #Do this if it is a source node
        #print(node, ' is a source node')
        #still need to update the out edges
        outDeg = list(G.out_edges(node))
        nodeContent = 1/len(SourceN)
        nextNodes = []
        sumOutWeights = 0
        #Get the Cumulative weights of all the outdoing edges
        for oD in outDeg: sumOutWeights += G[oD[0]][oD[1]]['Weight']
        #print('SumOutWeights is ', sumOutWeights)
        #sum the weights of the outgoing edges
        
        my_dictionary = {node:nodeContent}
        nx.set_node_attributes(G, my_dictionary, 'nodeTransferWeight')

        for oD in outDeg:
                #Get the weighted value of the output edge by dividing it's aerobic edge weight by the sum of all out edge weights
                nextEdgeGets = nodeContent* (G[oD[0]][oD[1]]['Weight']/sumOutWeights)
                #print('NextEdgeGets is ', nextEdgeGets)
                my_dictionary = {oD: nextEdgeGets}
                
                #put the data into the out edges
                nx.set_edge_attributes(G,  my_dictionary,'edgeTransferWeight')
                #print('set edge ', my_dictionary)
                nextNodes.append(oD[1])
                #let's print everything out for backup. 
                NodeStatus = 1
        #return if this was a source node
        return(G,nextNodes, NodeStatus) 
    
    else:
    #Do all this if it is not a source node
                
        wwIn = 0
        for i in inDeg:
            #check if all the incoming edgs are assigned yet.
            #print('The in degrees are ', i)
            if 'edgeTransferWeight' in G[i[0]][node]:              
                wwIn += G[i[0]][node]['edgeTransferWeight']
            else:
                #print('could not process node')
                NodeStatus = 0
                #Exit the function if you can't process this node
                return(G, [node], NodeStatus)
        #End of the inDeg for loop
        

        #if you get to here then all of the edgelists were assigned
        outDeg = list(G.out_edges(node))
        #if no OutDeg then it's a sinkNode
                    
        nextNodes = []
        nodeContent = wwIn
        #Still need to assign the value of the node
        sumOutWeights = 0
        #Get the Cumulative weights of all the outdoing edges
        for oD in outDeg: sumOutWeights += G[oD[0]][oD[1]]['Weight']
        
        my_dictionary = {node:nodeContent}
        nx.set_node_attributes(G, my_dictionary, 'nodeTransferWeight')
        
        if len(outDeg) == 0:
            NodeStatus = 2 #terminal Sink Node, successfull processed
            #Can't set any edges
        else:
            NodeStatus = 1 #successfully processed middle node
            
            for oD in outDeg:
                nextEdgeGets = nodeContent* (G[oD[0]][oD[1]]['Weight']/sumOutWeights)
                #print('NextEdgeGets is ', nextEdgeGets)
                my_dictionary = {oD: nextEdgeGets}              
                #put the data into the out edges
                nx.set_edge_attributes(G,  my_dictionary,'edgeTransferWeight')
#                print('set edge ', my_dictionary)
                nextNodes.append(oD[1])
                #let's print everything out for backup.
        return(G, nextNodes, NodeStatus)
#**************************************
os.chdir(r'C:\ResearchWorkingDirectory\gitRepositories\networksPlotting')

#Read in the raw edgelist
df = pd.read_excel(r"C:\ResearchWorkingDirectory\gitRepositories\RawExampleData.xlsx", sheet_name = 'Sheet1')    

#also need the file that matches the rules with the proper weighting
#select whichever weighting scheme you desire
rulesRef = pd.read_excel(r"C:\ResearchWorkingDirectory\gitRepositories\networkAnalysis\RuleData.xls", sheet_name = '110100')
       
rulesRefDict = rulesRef.set_index('Rule').T.to_dict('list')        
        
#Could also use this function to add a lot of different stuff
#*******************************************
for i, row in df.iterrows():
    tempRow = rulesRefDict[row['Rule']]
    df.at[i, 'Weight' ] = (float(tempRow[0])*10)/10
#********************************************
#currently the weighting reference is set up as costs, so need to invert to make it strengths
df['Weight'] = 1/df['Weight']


G = nx.from_pandas_edgelist(df, source ='From', target = 'To', edge_attr=True, create_using=nx.DiGraph())

[SourceN, SinkN] = FindSourceSink(G)

#initialize the "thisLevelNodes" variable
#**********************************
nextLevelNodes = SourceN
allSinks = []
nextNodes = []
nodeContent = 1/len(SourceN)
iterator = 0

#*******************************
while set(allSinks) != set(SinkN):
    thisLevelNodes = list(set(nextLevelNodes))
    nextLevelNodes = []
    for node in thisLevelNodes:
        [G, nextNodes, NodeStatus] = assignFlowWeights(G, node)
        #print('The next NodeStatus is ', NodeStatus)
        
        #if assignFlowWeights returns an empty array then we have a Terminal Sink. Otherwise we have an intermediate node.
        
        if NodeStatus == 0: #send this node around again
            #was not able to process node keep it for the next round
            if node not in nextLevelNodes:
                nextLevelNodes = nextLevelNodes + nextNodes
                #print('failed to process node ', node)
        elif NodeStatus == 1:
            nextLevelNodes = nextLevelNodes + nextNodes
            #print('The nextNodes is', nextNodes)
        elif NodeStatus == 2:
            if node not in allSinks:
                allSinks.append(node)
                #print(node, ' is a Terminal Sink')
    #End the For Loop for this Level Nodes
#*****************************************


#Display code if you need a visual, normally commented out
#***************************************
#nodes = G.nodes()

#pos = hierarchy_posMS(G)
#posIn = pd.read_excel(r"C:\ResearchWorkingDirectory\DataReferenceFiles\PaperData\\PlottingImages\PhenPos_PerOnly.xlsx", sheet_name = 'Sheet1')
#pos = {}
#for pIn in posIn.iterrows():
#    pos[pIn[0]] = list(pIn[1])    
   
#ww = nx.get_edge_attributes(G,'edgeTransferWeight')
#colors = list(nx.get_node_attributes(G,'nodeTransferWeight').values())
#labels = []
    
#bb = dict(G.out_degree() )
#colors = list(bb.values())

#nc = nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=colors, with_labels=False, node_size=100, cmap=plt.get_cmap('jet'), vmin = 0, vmax = 0.5, zorder =0)
#ec = nx.draw_networkx_edges(G, pos, alpha=0.8, width = list(ww.values()))
#ec = nx.draw_networkx_edges(G, pos, alpha=0.8, width = 1)
#plt.colorbar(nc, shrink=0.8)
#plt.show()
#*********************************

#check that the sinks sum to 1
zz2 = nx.get_node_attributes(G,'nodeTransferWeight')

tempZ = 0
for s in SinkN:
    tempZ += zz2[s]
print ('The sum of the sinks is ', tempZ)


#need some write out code here

#quick script to reconstruct the dataframe from the graph object so that I can write it out and pass it to other scripts
df = nx.to_pandas_edgelist(G)

##write out the EdgelistFile so that it can be read in by one of the analyzer programs. 
nodeTransWeights = nx.get_node_attributes(G,'nodeTransferWeight')

#add the betweenness here also - This is a very slow calculation - be patient. It is in this script so that we only have to do it once. 
btw = nx.betweenness_centrality(G)

nodeTransWeights2 = pd.DataFrame(list(nodeTransWeights.items()), columns = ['SMILES','nodeTransferWeights'])

for i, row in nodeTransWeights2.iterrows():
    nodeTransWeights2.at[i, 'Betweenness Centrality'] = btw[row['SMILES']]


##Add the degrees in and out ontol the nodes dataframe
#outDeg = pd.DataFrame.from_dict(dict(G.out_degree()),  orient = 'index', columns = ['OutDeg'])
#outDeg['SMILES'] = outDeg.index
#
#inDeg = pd.DataFrame.from_dict(dict(G.in_degree()),  orient = 'index', columns = ['InDeg'])
#inDeg['SMILES'] = inDeg.index
#
#nodeTransWeights2 = pd.merge(nodeTransWeights2,outDeg, on = 'SMILES' )
#nodeTransWeights2 = pd.merge(nodeTransWeights2,inDeg, on = 'SMILES' )


#Go back though the final dataframe and assign the node weights. Will need two looks one for "From" and one for "To"

for i, row in df.iterrows():
    df.at[i, 'fromNodeTransWeight' ] = nodeTransWeights[row['From']]
    df.at[i, 'ToNodeTransWeight' ] = nodeTransWeights[row['To']]
    df.at[i, 'ToBetweenessCentrality'] = btw[row['To']]
    df.at[i, 'FromBetweenessCentrality'] = btw[row['From']]
    

#Set the output directory
#os.chdir('C:\ResearchWorkingDirectory\DataReferenceFiles\PaperData')
#


#write out a full edgelist
#***************************
#timestr = 'deltaGEdges' + time.strftime("%Y%m%d-%H%M%S") +'.xlsx'
#print (timestr)
#writer = pd.ExcelWriter(timestr)
#df.to_excel(writer,'Sheet1')
#writer.save()


#write out the QF node file. This is the most important file for future calculations
os.chdir(r'C:\ResearchWorkingDirectory\gitRepositories')

#timestr = 'ExampleNT_BTW_Nodes' + time.strftime("%Y%m%d-%H%M%S") +'.xlsx'
timestr = 'ExampleNT_BTW_Nodes' + '.xlsx'

print (timestr)
writer = pd.ExcelWriter(timestr)
nodeTransWeights2.to_excel(writer,'Sheet1')
writer.save()

