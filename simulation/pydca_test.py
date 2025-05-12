# import pydca modules
from pydca.meanfield_dca import meanfield_dca
import os
from tqdm_joblib import tqdm_joblib
from tqdm import tqdm
from joblib import Parallel, delayed

# get all simulated files in the simulation directory
path = '/users/bag/hlq763/hbv_covar3/analysis/sim_seq/'
entries = os.listdir(path)
sim_dirs = [os.path.join(path, entry) for entry in entries if os.path.isdir(os.path.join(path, entry))]

# sim_files = [os.path.join(dir, f)for dir in sim_dirs for f in os.listdir(dir)  if f.endswith('rescaled.fasta')]
#

def process_sim_file(sim_seq_file):
    try:
        # Define the names of the converted and output files
        converted_seq_file = sim_seq_file.replace('.fasta', '_converted.fasta')
        output_file = converted_seq_file.replace('_converted.fasta', '_mfdca.csv').replace('simseq', 'simresult')
        # Read the input FASTA file and write the converted file
        with open(sim_seq_file, 'r') as input_file, open(converted_seq_file, 'w') as converted_file:
            for line in input_file:
                if not line.startswith('>'):
                    # Replace 'x' with 'a' and 'y' with 'c' in the sequence
                    modified_line = line.replace('x', 'a').replace('y', 'c')
                    converted_file.write(modified_line)
                else:
                    # Write header lines as-is
                    converted_file.write(line)
        # Compute Mean Field DCA
        mfdca_inst = meanfield_dca.MeanFieldDCA(
            converted_seq_file,
            'rna',
            pseudocount=0.5,
            seqid=0.8,
        )
        mfdca_FN_APC = mfdca_inst.compute_sorted_FN_APC()
        link1_2_rank = [i for i, ((x, y), z) in enumerate(mfdca_FN_APC) if x==0 and y==1][0]
        link2_3_rank = [i for i, ((x, y), z) in enumerate(mfdca_FN_APC) if x==1 and y==2][0]
        # Write the top 0.05% site pairs and scores to the output CSV
        # with open(output_file, "w") as output_file_handle:
        #         output_file_handle.write(str(link1_2_rank + 1) +"\n")
        #         output_file_handle.write(str(link2_3_rank + 1)+ "\n")
        alpha = 0.05
        alpha_index = int(alpha * len(mfdca_FN_APC))
        for site_pair, score in mfdca_FN_APC[:alpha_index]:
            output = ",".join([str(i) for i in list(site_pair)] + [str(score)])
            print(output, file=open(output_file, "a"))
        # delete converted file
        os.remove(converted_seq_file)
    except Exception as e:
        print("Error processing file %s: %s" % (sim_seq_file, e))

def convert_sim_file(sim_seq_file):
    converted_seq_file = sim_seq_file.replace('.fasta', '_converted.fasta')
    with open(sim_seq_file, 'r') as input_file, open(converted_seq_file, 'w') as output_file:
        for line in input_file:
            if not line.startswith('>'):
                modified_line = line.replace('x', 'a').replace('y', 'c')
                output_file.write(modified_line)
            else:
                output_file.write(line)

                


def gather_sim_files(sim_dirs):
    sim_files = []
    for directory in sim_dirs:
        try:
            sim_files.extend([os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('rescaled.fasta') ])
        except Exception as e:
            print("Error accessing directory %s: %s" % (directory, e))
    return sim_files


# Define your simulation directories
# Gather all .fasta files
sim_files = gather_sim_files(sim_dirs)

if not sim_files:
    print("No .fasta files found. Exiting.")

for sim_file in tqdm(sim_files):
    process_sim_file(sim_file)

print("All files have been processed.")

# iterate through all the files and convert them
for sim_seq_file in sim_files:
    print(sim_seq_file, flush=True)
    # Define the name of the converted file
    converted_seq_file = sim_seq_file.replace('.fasta', '_converted.fasta')
    output_file = converted_seq_file.replace('_converted.fasta', '_mfdca.csv').replace('simseq', 'simresult')
    # Read in the fasta file
    with open(sim_seq_file, 'r') as input_file, open(converted_seq_file, 'w') as output_file:
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
    mfdca_inst = meanfield_dca.MeanFieldDCA(
        converted_seq_file,
        'rna',
        pseudocount = 0.5,
        seqid = 0.8,
    )
    mfdca_FN_APC = mfdca_inst.compute_sorted_FN_APC()
    alpha = 0.05
    alpha_index = int(alpha * len(mfdca_FN_APC))
    for site_pair, score in mfdca_FN_APC[:alpha_index]:
        output = ",".join(list(site_pair) + [score])
        print(output, file=open(output_file, "a"))




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
    mfdca_FN_APC = mfdca_inst.compute_sorted_FN_APC()

    # grab the scores for a 

    for site_pair, score in mfdca_FN_APC[:100]:
        output = ",".join(list(site_pair) + [score])
        print(output)
        print(output, file=open(output_file, "a"))

mfdca_inst = meanfield_dca.MeanFieldDCA(
    converted_seq_file,
    'rna',
    pseudocount=0.5,
    seqid=0.8,
)

