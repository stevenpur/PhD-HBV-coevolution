from pydca.meanfield_dca import meanfield_dca

import os
from tqdm_joblib import tqdm_joblib
from tqdm import tqdm
from joblib import Parallel, delayed


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

