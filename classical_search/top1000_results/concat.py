# use this to concatenate files

import os

def concatenate_files(directory, output_file):
    """
    Concatenates the content of all files in a directory into a single output file.

    Args:
        directory (str): The path to the directory containing the files.
        output_file (str): The path to the output file where the concatenated content will be written.
    """
    with open(output_file, 'w') as outfile:
        for filename in os.listdir(directory):
            file_path = os.path.join(directory, filename)
            if os.path.isfile(file_path):
                with open(file_path, 'r') as infile:
                    outfile.write(infile.read())
                    outfile.write('\n')  # Add a newline between files (optional)

# Example usage:
directory_path = '/Users/raubenb/Desktop/gitrepos/protein-folding-qc/classical_search/top1000_results/top1000_chig_mr_1nn/RMSD_of_the_aligned_structure_collection_20062'  # Replace with the actual path to your directory
output_file_path = '/Users/raubenb/Desktop/gitrepos/protein-folding-qc/classical_search/top1000_results/top1000_chig_mr_1nn/rmsd.tabular'  # Replace with the desired output file path
concatenate_files(directory_path, output_file_path)