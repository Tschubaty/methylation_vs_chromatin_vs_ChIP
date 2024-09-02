import os
import pandas as pd
import re
import numpy as np
import math

# Base path to your meme-chip results
base_path = r"D:/Users/Daniel Batyrev/Documents/GitHub/meme/single_zip"

# Base path to HOCOMOCOv11_core
HOCOMOCOv11_path = "C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/meme/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme"


def check_motif_contains_CG(motif_id, protein_name, experiment_id, bio_sample):
    if not isinstance(motif_id, str):
        print(f"Motif is empty: {motif_id}")
        return ""

    # Determine the file path based on the motif ID pattern
    if motif_id.isupper() and motif_id.isalpha():
        file_path = os.path.join(base_path, f"{bio_sample}_{protein_name}-human_{experiment_id}.fa",
                                 f"{bio_sample}_{protein_name}-human_{experiment_id}.fa", "meme_out", "meme.txt")
    elif re.match(r"^\d+-[A-Z]+", motif_id):
        file_path = os.path.join(base_path, f"{bio_sample}_{protein_name}-human_{experiment_id}.fa",
                                 f"{bio_sample}_{protein_name}-human_{experiment_id}.fa", "streme_out", "streme.txt")
    else:
        file_path = HOCOMOCOv11_path

    # Read the content of the file
    try:
        with open(file_path, 'r') as file:
            content = file.read()
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return ""

    # Search for the motif section in the file content
    motif_section_pattern = rf"MOTIF {motif_id}.*?letter-probability matrix:(.*?)(?=\nMOTIF|\n\n|\n-|\Z|\nURL)"
    match = re.search(motif_section_pattern, content, re.DOTALL)

    if match:
        matrix_content = match.group(1).strip().split('\n', 1)[1]  # Remove the first line
        rows = matrix_content.strip().split('\n')
        probabilities = np.array([list(map(float, row.split())) for row in rows])

        # Check for CG pattern
        for i in range(len(probabilities) - 1):
            if probabilities[i, 1] == max(probabilities[i, :]) and probabilities[i + 1, 2] == max(
                    probabilities[i + 1, :]):
                return True
        return False
    else:
        print(f"Motif section for {motif_id} not found in {file_path}.")
        return ""


def check_CG_pattern(row):
    return check_motif_contains_CG(row['single Motif ID'], row['Protein'], row['Experiment ID'], row['Bio Sample'])


# Load the data from the specified CSV file
input_file = r"C:\Users\Daniel Batyrev\Documents\GitHub\methylation_vs_chromatin_vs_ChIP\Encode3\meme\meme-chip_results\all_single_experiment_motifs_groups_all.csv"
all_groups_df = pd.read_csv(input_file)

# Add the 'Contains CG' column by applying the check_CG_pattern function to each row
all_groups_df['Contains CG'] = all_groups_df.apply(check_CG_pattern, axis=1)

# Drop the "FIMO Motif source file" column from the DataFrame
if 'FIMO Motif source file' in all_groups_df.columns:
    all_groups_df = all_groups_df.drop(columns=['FIMO Motif source file'])

# Save the updated DataFrame to a new CSV file
output_file = r"C:\Users\Daniel Batyrev\Documents\GitHub\methylation_vs_chromatin_vs_ChIP\Encode3\meme\meme-chip_results\all_single_experiment_motifs_groups_all_with_CG.csv"
all_groups_df.to_csv(output_file, index=False)

print(f"All motifs with CG information saved to: {output_file}")
