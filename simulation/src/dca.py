# import pydca modules
from pydca.meanfield_dca import meanfield_dca
import os
from tqdm_joblib import tqdm_joblib
from tqdm import tqdm
from joblib import Parallel, delayed
from dca_threshold import six_sigma_cutoff

def convert_sequence_file(input_file, output_file):
    """Convert sequence file by replacing 'x' with 'a' and 'y' with 'c'."""
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if not line.startswith('>'):
                modified_line = line.replace('x', 'a').replace('y', 'c')
                outfile.write(modified_line)
            else:
                outfile.write(line)

def process_sim_file(sim_seq_file):
    """Process a single simulation file through the DCA pipeline."""
    try:
        # Define file paths
        converted_seq_file = sim_seq_file.replace('.fasta', '_converted.fasta')
        # Create output filename with proper directory structure
        rel_path = os.path.relpath(sim_seq_file, os.path.dirname(os.path.dirname(sim_seq_file)))
        output_file = os.path.join(output_dir, rel_path.replace('.fasta', '_mfdca.csv'))
        
        # Ensure output directory exists
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Convert sequence file
        convert_sequence_file(sim_seq_file, converted_seq_file)
        
        # Compute Mean Field DCA
        mfdca_inst = meanfield_dca.MeanFieldDCA(
            converted_seq_file,
            'rna',
            pseudocount=0.5,
            seqid=0.8,
        )
        
        # Get DCA results
        mfdca_FN_APC = mfdca_inst.compute_sorted_FN_APC()

        # calculate the threshold for true direct coupling
        threshold = six_sigma_cutoff(mfdca_FN_APC)
        # get the site pairs with scores above the threshold
        site_pairs = [site_pair for site_pair, score in mfdca_FN_APC if score > threshold]

        # Write results to output file
        with open(output_file, "w") as outfile:
            for site_pair in site_pairs:
                output = ",".join([str(i) for i in list(site_pair)])
                outfile.write(output + "\n")
        
        # Clean up temporary file
        os.remove(converted_seq_file)
        
    except Exception as e:
        print(f"Error processing file {sim_seq_file}: {e}")

def gather_sim_files(sim_dirs):
    """Gather all simulation files from the given directories."""
    sim_files = []
    for directory in sim_dirs:
        try:
            sim_files.extend([
                os.path.join(directory, f) 
                for f in os.listdir(directory) 
                if f.endswith('rescaled.fasta')
            ])
        except Exception as e:
            print(f"Error accessing directory {directory}: {e}")
    return sim_files

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Process simulation files using DCA.')
    parser.add_argument('--simdir', type=str, required=True,
                      help='Directory containing simulation files')
    parser.add_argument('--output', type=str, required=True,
                      help='Output directory for DCA results')
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_args()
    
    # Get simulation directories
    entries = os.listdir(args.simdir)
    sim_dirs = [
        os.path.join(args.simdir, entry) 
        for entry in entries 
        if os.path.isdir(os.path.join(args.simdir, entry))
    ]
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    
    # Gather simulation files
    sim_files = gather_sim_files(sim_dirs)
    
    if not sim_files:
        print("No .fasta files found. Exiting.")
        return
    
    # Process files in parallel
    with tqdm_joblib(tqdm(total=len(sim_files), desc="Processing Files")):
        Parallel(n_jobs=6)(
            delayed(process_sim_file)(f, args.output) for f in sim_files
        )
    
    print("All files have been processed.")

if __name__ == "__main__":
    main()

