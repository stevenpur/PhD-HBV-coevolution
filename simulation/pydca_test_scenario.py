# import pydca modules
from pydca.meanfield_dca import meanfield_dca
import os
from tqdm_joblib import tqdm_joblib
from tqdm import tqdm
from joblib import parallel, delayed
import sys
from dca_threshold import six_sigma_cutoff

# take user parameters from command line

sim_files = [sys.argv[1]]


# define functions
def convert_sim_file(sim_seq_file):
    converted_seq_file = sim_seq_file.replace('.fasta', '_converted.fasta')
    with open(sim_seq_file, 'r') as input_file, open(converted_seq_file, 'w') as output_file:
        for line in input_file:
            if not line.startswith('>'):
                modified_line = line.replace('x', 'a').replace('y', 'c')
                output_file.write(modified_line)
            else:
                output_file.write(line)

# sim_files =  ['/users/bag/hlq763/hbv_covar3/analysis/sim_seq/test/test_comut_n1600u5f16per0.5.fasta']

# sim_files = [os.path.join(dir, f)for dir in sim_dirs for f in os.listdir(dir)  if f.endswith('rescaled.fasta')]
#
for sim_file in sim_files:
    convert_sim_file(sim_file)

# def gather_sim_files(sim_dirs):
#     sim_files = []
#     for directory in sim_dirs:
#         try:
#             sim_files.extend([os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('rescaled.fasta') ])
#         except Exception as e:
#             print("Error accessing directory %s: %s" % (directory, e))
#     return sim_files


# Define your simulation directories
# Gather all .fasta files
# sim_files = gather_sim_files(sim_dirs)

if not sim_files:
    print("No .fasta files found. Exiting.")

print("All files have been processed.")

# iterate through the converted files and run mfDCA
converted_seq_files = [f.replace('.fasta', '_converted.fasta') for f in sim_files]

for i, converted_seq_file in enumerate(converted_seq_files):
    # Define the name of the output file
    # Run meanfield DCA
    for seqid in [1]:
        print(seqid, flush = True)
        output_file = converted_seq_file.replace('_converted.fasta', '_mfdca_'+str(seqid)+'.csv')
        print(seqid)
        mfdca_inst = meanfield_dca.MeanFieldDCA(
            converted_seq_file,
            'rna',
            pseudocount = 0.5,
            seqid = seqid,
        )
        # Save top 100 contacts
        mfdca_FN_APC = mfdca_inst.compute_sorted_FN_APC()

        # gather the scores
        scores = [score for site_pair, score in mfdca_FN_APC]
        for i in range(4):
            scores[i] = 200000 - 0.05*i
        # get the threshold
        threshold, keep_mask = six_sigma_cutoff(scores)

        for i, (site_pair, score) in enumerate(mfdca_FN_APC):
            if keep_mask[i]:
                output = str(site_pair[0]) + ',' + str(site_pair[1]) + ',' + str([score])
                print(output, file=open(output_file, "a"))
        
        for i, (site_pair, score) in enumerate(mfdca_FN_APC):
            output_file2 = output_file.replace('.csv', '_filtered.csv')
            if not keep_mask[i]:
                output = str(site_pair[0]) + ',' + str(site_pair[1]) + ',' + str([score])
                print(output, file=open(output_file2, "a"))


