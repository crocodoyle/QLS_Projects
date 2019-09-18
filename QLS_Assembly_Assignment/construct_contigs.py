# -*- coding: utf-8 -*-
"""
Created on Sat Sep 14 22:05:09 2019

@author: Bowei
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2

import numpy as np
import itertools
import re
import argparse
import time
import matplotlib.pyplot as plt

import os
# os.chdir('D:\\xbw\\my docs\\McGill\\OneDrive - McGill University\\year1\\QLSC600\\module1\\hw1')

from assembly import *


def create_prefix_dict(seq_list, w, pre=True):
    '''
    This makes use of the 'Prefix Index' to construct overlap dictionary
    seq_list: list of SeqRecord objects that stores all the reads
    w: int; how many chars used as prefix/suffix to define overlap
    pre: Boolean. Default=True. Will do prefix if true and suffix otherwise.
    This will return a dictionary recording the preInx
    '''
    out_dict = {}
    #Each fake read is ordered so ReadX is just X element in the list
    for seq_idx, seq in enumerate(seq_list):
        if pre:
            #first w chars
            prefix_idx = str(seq.seq[0:w])
        else:
            #last w chars
            prefix_idx = str(seq.seq[-w:])
            
        if prefix_idx in out_dict:
            out_dict[prefix_idx].append(seq_idx)
        else:
            out_dict[prefix_idx] = [seq_idx]

    return out_dict

    
    
def connect_contig(seq_list, input_contig, pre_dict, suf_dict, w):
    '''
    This function connects 2 sequence together if they overlap
    it will return the coordinates where the cutoff point is.
    This function assumed the overlap exists.
    seq_list: list of SeqRecord objects that stores all the reads
    seq_in: The current connected chain
    pre_dict, suf_dict: Prefix dictionary and suffix dictionary created
    w: int; how many chars used as prefix/suffix to define overlap. Passed down
    This will return a list where each element is in the format of [(read#,(start_idx, end_idx))]
    Indices are following py convention so [start_idx, end_idx)
    '''
    
    #We are using the fact that all sequence in seqList has the same length
    seq_length = len(seq_list[0])
    output_contig = input_contig

    # Left edge of the connected chain
    left_edge = str(seq_list[input_contig[0][0]].seq)[0:w]
    #right edge of the connected chain
    right_edge=str(seq_list[input_contig[-1][0]].seq)[-w:]
        
    if (suf_dict.get(left_edge) != None):
    #avoid key error
        if (len(suf_dict[left_edge]) == 1):
        # meaning found unique contigs to stitch to the left
            output_contig.insert(0,(suf_dict[left_edge][0], (0, seq_length - w)))
    
    if (pre_dict.get(right_edge) != None):    
    #avoid key error
        if (len(pre_dict[right_edge]) == 1):
        # meaning found unique contigs to stitch to the right
            output_contig.append((pre_dict[right_edge][0], (w, seq_length)))
        
    return output_contig
    
    
def assemble_genome(seq_list, w):
    '''Main function that build up the full contigs
    seq_list: fasta file that was read in using the read function
    This outputs a dictionary saying how many contigs it forms and how is it finally connected
    '''
    node_queue = set(range(len(seq_list)))
    output_dict = {}
    n_nodes = 1

    prefix_dict = create_prefix_dict(seq_list, w, pre=True)
    suffix_dict = create_prefix_dict(seq_list, w, pre=False)

    while node_queue:
        #randomly pick a starting point
        start_seq = node_queue.pop()
        seq_in = [(start_seq, (0, len(seq_list[start_seq])))]
        expand_contig = connect_contig(seq_list, seq_in, prefix_dict, suffix_dict, w)
        
        while (expand_contig != seq_in):
            # we can still expand
            seq_in = expand_contig
            expand_contig = connect_contig(seq_list, seq_in, prefix_dict, suffix_dict, w)
        
        output_dict[n_nodes] = expand_contig
        
        #remove the connected nodes from the two dicts and nodeQueue
        used_node_list = [k[0] for k in expand_contig]
        for key in prefix_dict:
            prefix_dict[key] = list(set(prefix_dict[key]) - set(used_node_list))
        for key in suffix_dict:
            suffix_dict[key]=list(set(suffix_dict[key]) - set(used_node_list))

        node_queue = node_queue.difference(set(used_node_list))
        n_nodes += 1

    return output_dict

def outputSequence(seq_key,seq_list,output_dict):
    '''
    Take the key of the output_dict and output the actual sequence that's stitched
    '''
    seq_coord=output_dict[seq_key]  
    outSeq=''
    for piece in seq_coord:
        outSeq += seq_list[piece[0]].seq[piece[1][0]:piece[1][1]]
    
    return outSeq

def frac_covered_in_genome(outSeq,genome):
    """
    Calculated what percentage of genome is covered by the union of the constructed contigs
    """
    foundPos=[0]*len(genome); errorCounter=0
    for contig in outSeq:
        try:
            startInx=str(genome).index(str(contig))
            foundPos[startInx:startInx+len(contig)]=[1]*len(contig)
        except ValueError:
            errorCounter += 1
    return sum(foundPos)/len(genome), errorCounter / len(outSeq)
    
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    #minimal eg
    #CATTCGAATA

    random_genome_filename = 'test1.fa'
    bacteria_genome_filename = 'bacteria_5_3061335.fa'


    read_length = 100
    coverage = 10

    # genome, records, output_filename = generate_reads(random_genome_filename, read_length, coverage)
    # genome_len=len(genome)
    # w = np.log10(genome_len / coverage) / np.log10(4)
    # print('w:', w, int(w))
    # w = int(w)
    # print('w is now used as:', w)
    #
    #
    # print('Assembling contigs from', len(records), 'records')
    # startTime=time.time()
    # contigs = assemble_genome(records, w)
    # print("Took %s seconds to run" % (time.time() - startTime))
    #
    # print('Assembled', len(contigs), 'contigs:')
    # #print(contigs)
    # outSeq=[]
    # for contig in contigs:
    #     outSeq.append(outputSequence(contig,records,contigs))

    coverages = [1, 5, 10]
    read_lengths = [100]

    compare_parameters(coverages, read_lengths)
    
    # To note: How to actually get the sequence
    # Add a function that acutally writes out the letters
    if False:
        f=open('assembled_contigs.fa','w')
        for contig in contigs:
            f.write('>Contig'+str(contig)+'\n')
            f.write(str(outputSequence(contig,records,contigs))+'\n')
        
        f.close()
