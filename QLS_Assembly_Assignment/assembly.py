from Bio import SeqIO
import numpy as np

random_genome_filename = 'randomGenome.fa'
bacteria_genome_filename = 'bacteria_5_3061335.fa'


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

    return len(genome), records, output_filename
