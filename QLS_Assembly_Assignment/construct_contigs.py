# -*- coding: utf-8 -*-
"""
Created on Sat Sep 14 22:05:09 2019

@author: Bowei
"""

from Bio import SeqIO
from Bio.Seq import Seq

import numpy as np
import itertools
import re

def createPreInxDict(seqList,w,pre=True):
    '''
    This makes use of the 'Prefix Index' to construct overlap dictionary
    seqList: list of SeqRecord objects that stores all the reads
    w: int; how many chars used as prefix/suffix to define overlap
    pre: Boolean. Default=True. Will do prefix if true and suffix otherwise.
    This will return a dictionary recording the preInx
    '''
    outDict={}
    #Each fake read is ordered so ReadX is just X element in the list
    for seqInx in range(len(seqList)):
        if pre:
            #first w chars
            preInx=str(seqList[seqInx].seq)[0:w]
        else:
            #last w chars
            preInx=str(seqList[seqInx].seq)[-w:]
            
        if preInx in outDict:
            outDict[preInx].append(seqInx)
        else:
            outDict[preInx]=[seqInx]
    return outDict

    
    
def connect2Contigs(seqList,seqIn,preDict,sufDict,w):
    '''
    This function connects 2 sequence together if they overlap
    it will return the coordinates where the cutoff point is.
    This function assumed the overlap exists.
    seqList: list of SeqRecord objects that stores all the reads
    seqIn: The current connected chain
    preDict,sufDict: Prefix dictionary and suffix dictionary created
    w: int; how many chars used as prefix/suffix to define overlap. Passed down
    This will return a list where each element is in the format of [(read#,(startInx,endInx))]
    Indices are following py convention so [startInx,endInx)
    '''
    
    #We are using the fact that all sequence in seqList has the same length
    seqLength=len(seqList[0])
    seqOut=seqIn
    # Left edge of the connected chain
    leftEdge=str(seqList[seqIn[0][0]].seq)[0:w]
    #right edge of the connected chain
    rightEdge=str(seqList[seqIn[-1][0]].seq)[-w:]
        
    if (sufDict.get(leftEdge) != None):
    #avoid key error
        if (len(sufDict[leftEdge])==1):
        # meaning found unique contigs to stitch to the left
            seqOut.insert(0,(sufDict[leftEdge][0],(0,seqLength-w)))
    
    if (preDict.get(rightEdge) != None):    
    #avoid key error
        if (len(preDict[rightEdge])==1):
        # meaning found unique contigs to stitch to the right
            seqOut.append((preDict[rightEdge][0],(w,seqLength)))
        
    return seqOut
    
    
def stitchContigs(seqList,w):
    '''Main function that build up the full contigs
    seqList: fasta file that was read in using the read function
    This outputs a dictionary saying how many contigs it forms and how is it finally connected
    '''
    nodeQueue=set(range(len(seqList))); output_dict={}
    nodeCount=1
    preInxDict=createPreInxDict(seqList,w,pre=True)
    sufInxDict=createPreInxDict(seqList,w,pre=False)

    while nodeQueue:
        #randomly pick a starting point
        startSeq=nodeQueue.pop(); seqIn=[(startSeq,(0,len(seqList[startSeq])))]
        expandContig=connect2Contigs(seqList,seqIn,preInxDict,sufInxDict,w)
        
        while (expandContig != seqIn):
            # we can still expand
            seqIn=expandContig
            expandContig=connect2Contigs(seqList,seqIn,preInxDict,sufInxDict,w)
        
        output_dict[nodeCount]=expandContig
        
        #remove the connected nodes from the two dicts and nodeQueue
        usedNodeList=[k[0] for k in expandContig]
        for key in preInxDict:
            preInxDict[key]=list(set(preInxDict[key])-set(usedNodeList))
        for key in sufInxDict:
            sufInxDict[key]=list(set(sufInxDict[key])-set(usedNodeList))

        nodeQueue=nodeQueue.difference(set(usedNodeList))
        nodeCount+=1

    return output_dict


#minmal example to run
if False:
    #minimal eg
    #CATTCGAATA
    seqList=[SeqIO.FastaIO.SeqRecord(seq=Seq('TTC'),id='Read0'),SeqIO.FastaIO.SeqRecord(seq=Seq('ATT'),id='Read1'),
             SeqIO.FastaIO.SeqRecord(seq=Seq('GAA'),id='Read2'),SeqIO.FastaIO.SeqRecord(seq=Seq('TCG'),id='Read3'),
             SeqIO.FastaIO.SeqRecord(seq=Seq('CAT'),id='Read4'),SeqIO.FastaIO.SeqRecord(seq=Seq('AAT'),id='Read5')]
    #b/c both Read4 and Read5 contains AT, we shouldn't stitch either onto Read1
    w=2
    stitchContigs(seqList,w)

# To note: How to actually get the sequence
# Add a function that acutally writes out the letters
