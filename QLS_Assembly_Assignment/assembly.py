from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

from construct_contigs import *


random_genome_filename = 'randomGenome.fa'
bacteria_genome_filename = 'bacteria_5_3061335.fa'


def compare_parameters(coverages, read_lengths):

    for read_length in read_lengths:
        n_contigs = []
        n50s = []
        genome_covered = []
        erroneous_contigs = []

        for coverage in coverages:
            fig, axes = plt.subplots(1, len(coverages))

            genome, records, output_filename = generate_reads(random_genome_filename, read_length, coverage)

            w = int(np.log10(len(genome) / coverage) / np.log10(4))

            contigs = assemble_genome(records, w)
            out_seq = []
            for contig in contigs:
               out_seq.append(outputSequence(contig, records, contigs))


            n_contigs.append(len(out_seq))
            print('Number of contigs:', len(out_seq))

            n50s.append(n50(out_seq, len(genome)))

            frac_covered, err_rate = frac_covered_in_genome(out_seq, genome)

            genome_covered.append(frac_covered)
            erroneous_contigs.append(err_rate)

        fig, ax = plt.subplots(1, 4, figsize=(12, 4), dpi=300)

        ax[0].plot(n_contigs, coverage)
        ax[0].set_xlabel('Coverage', fontsize=16)
        ax[0].set_ylabel('# contigs', fontsize=16)

        ax[1].plot(n50s, coverage)
        ax[1].set_xlabel('Coverage', fontsize=16)
        ax[1].set_ylabel('N50', fontsize=16)

        ax[2].plot(genome_covered, coverage)
        ax[2].set_xlabel('Coverage', fontsize=16)
        ax[2].set_ylabel('Fraction of bps in a contig', fontsize=16)

        ax[3].plot(erroneous_contigs, coverage)
        ax[3].set_xlabel('Coverage', fontsize=16)
        ax[3].set_ylabel('Contig Error Rate', fontsize=16)

        plt.plot(n_contigs, coverage)
        plt.show()
        fig.savefig(str(read_length) + '_stats.png', dpi=300)

def n50(contigs, genome_len):
    contig_lengths = []
    for contig in contigs:
        contig_lengths.append(len(str(contig)))

    print('contig lengths:', contig_lengths)

    contig_lengths.sort(reverse=True)

    print('contig lengths', contig_lengths)

    running_length = 0
    for contig_len in contig_lengths:
        running_length += contig_len
        if running_length > genome_len / 2:
            return contig_len





def generate_reads(fasta_filename, read_length, coverage=2):
    '''Reads a fasta file and generates a bunch of fake reads, outputing them to another fasta file'''

    # Read and parse the input fasta file
    parser = SeqIO.parse(fasta_filename, "fasta")

    for record in parser:
        genome = record.seq

    # Generate fake reads


    print('Length of genome:', len(genome))
    print('Read length inputted:', read_length)
    print('Coverage inputted:', coverage)
    num_reads = int((len(genome) * coverage) // read_length)

    print('Number of reads:', num_reads)

    # valid starting locations are from 0 to len(genome) - read_length
    starting_locations = np.random.randint(0, len(genome)-read_length, num_reads)

    records = []
    for i, starting_location in enumerate(starting_locations):
        fake_read = genome[starting_location:starting_location+read_length]

        fake_record = SeqIO.FastaIO.SeqRecord(fake_read, 'Read' + str(i), '', '')

        records.append(fake_record)

    # Generate an output filename from the input filename
    tokens = fasta_filename.split('.')
    output_filename = tokens[0] + '_output_reads.fa'

    # Write to file!
    fasta_writer = SeqIO.FastaIO.FastaWriter(open(output_filename, 'w'))

    fasta_writer.write_file(records)

    return genome, records, output_filename
