# use this to concatenate files

import os
import pandas as pd

def concatenate_files(directory, output_file):
    """
    Concatenates the content of all files in a directory into a single output file.

    Args:
        directory (str): The path to the directory containing the files.
        output_file (str): The path to the output file where the concatenated content will be written.
    """
    with open(output_file + 'rmsd.tabular', 'w') as outfile:
        for filename in os.listdir(directory):
            file_path = os.path.join(directory, filename)
            if os.path.isfile(file_path):
                with open(file_path, 'r') as infile:
                    outfile.write(infile.read())
                    outfile.write('\n')  # Add a newline between files (optional)

def rmsd_stats(rmsd_df, output_file_path, protein_model):

    df = pd.read_csv(rmsd_df)
    protein_model = protein_model
    min_rmsd = df.iloc[:,0].min().round(2)
    mean_rmsd = df.iloc[:, 0].mean().round(2)
    median_rmsd = df.iloc[:, 0].median().round(2)
    std_rmsd = df.iloc[:, 0].std().round(2)
    df_out = pd.DataFrame({'protein_model': [protein_model],
                           'min': [min_rmsd], 
                           'mean': [mean_rmsd],
                           'median': [median_rmsd],
                           'std': [std_rmsd]}
                          )
    df_out.to_csv(output_file_path + 'rmsd_stats.csv', index=False)
    

    # print(f"The mean of the second column is: {mean_rmsd}")

# Example usage:
directory_path = '/Users/raubenb/Desktop/gitrepos/protein-folding-qc/classical_search/all_klvffa_hp_fcc/RMSD of the aligned structure collection 14118'  # Replace with the actual path to your directory
output_file_path = '/Users/raubenb/Desktop/gitrepos/protein-folding-qc/classical_search/all_klvffa_hp_fcc/'  # Replace with the desired output file path
protein_model = 'klvffa_HP'
concatenate_files(directory_path, output_file_path)
rmsd_stats(rmsd_df=output_file_path + 'rmsd.tabular', output_file_path=output_file_path, protein_model=protein_model)