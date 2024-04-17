# import pydca modules
from pydca.plmdca import plmdca
from pydca.meanfield_dca import meanfield_dca
from pydca.sequence_backmapper import sequence_backmapper
from pydca.msa_trimmer import msa_trimmer
from pydca.contact_visualizer import contact_visualizer
from pydca.dca_utilities import dca_utilities
import os

# get all simulated files in the simulation directory
sim_dir = "/users/bag/hlq763/hbv_covar3/analysis/sim_seq/archive"
sim_files = [f for f in os.listdir(sim_dir) if f.endswith('.fasta')]

# iterate through all the files and convert them
for sim_seq_file in sim_files:
    # Define the name of the converted file
    sim_seq_path = sim_dir + "/" + sim_seq_file
    converted_seq_path = sim_seq_path.replace('.fasta', '_converted.fasta')
    # Read in the fasta file
    with open(sim_seq_path, 'r') as input_file, open(converted_seq_path, 'w') as output_file:
        # Iterate through each line in the input file
        for line in input_file:
            # Check if the line represents a sequence (as opposed to a header in FASTA format)
            if not line.startswith('>'):
                # Replace 'x' with 'a' and 'y' with 'c' in the sequence
                modified_line = line.replace('x', 'a').replace('y', 'c')
                output_file.write(modified_line)
            else:
                # If the line is a header, write it as is
                output_file.write(line)

# iterate through the converted files and run mfDCA
converted_seq_files = [f for f in os.listdir(sim_dir) if f.endswith('_converted.fasta')]
for i, converted_seq_file in enumerate(converted_seq_files):
    print("Running mfDCA on file " + i + " of ", len(converted_seq_files), "...", flush=True)
    # Define the name of the output file
    converted_seq_path = sim_dir + "/" + converted_seq_file
    output_file = converted_seq_file.replace('_converted.fasta', '_mfdca.csv')
    # Run meanfield DCA
    mfdca_inst = meanfield_dca.MeanFieldDCA(
        converted_seq_path,
        'rna',
        pseudocount = 0.5,
        seqid = 0.8,
    )
    # Save top 100 contacts
    mfdca_FN_APC = mfdca_inst.compute_sorted_FN_APC()
    print("hi")
    for site_pair, score in mfdca_FN_APC[:100]:
        output = ",".join(list(site_pair) + [score])
        print(output)
        print(output, file=open(output_file, "a"))