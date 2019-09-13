from Bio import SeqIO
import numpy as np

random_genome_filename = 'randomGenome.fa'
bacteria_genome_filename = 'bacteria_5_3061335.fa'


def fake_reads(genome_filename, read_length, coverage=2):
    '''Reads a fasta file and generates a bunch of fake reads, outputing them to another fasta file'''

    # Read and parse the input fasta file
    parser = SeqIO.parse(genome_filename, "fasta")

    for record in parser:
        genome = record.seq

    # Generate fake reads


    # print('length of genome:', len(genome))

    num_reads = (len(genome) * coverage) // read_length

    # valid starting locations are from 0 to len(genome) - read_length
    starting_locations = np.random.randint(0, len(genome)-read_length, num_reads)

    records = []
    for i, starting_location in enumerate(starting_locations):
        fake_read = genome[starting_location:starting_location+read_length]

        fake_record = SeqIO.FastaIO.SeqRecord(fake_read, 'Read' + str(i), '', '')

        records.append(fake_record)

    # Generate an output filename from the input filename
    tokens = genome_filename.split('.')
    output_filename = tokens[0] + '_fake_reads.fa'

    # Write to file!
    fasta_writer = SeqIO.FastaIO.FastaWriter(open(output_filename, 'w'))

    fasta_writer.write_file(records)

    return output_filename

print('hello world')


fake_reads(random_genome_filename, 5)