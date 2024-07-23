import os
import pandas as pd
import re
import numpy as np
import math

# Base path to your meme-chip results
base_path = r"D:/Users/Daniel Batyrev/Documents/GitHub/successZIP/extracted/"

# Base path to HOCOMOCOv11_core
HOCOMOCOv11_path = "C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/meme/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme"

# with open(HOCOMOCOv11_path, 'r') as file:
#    HOCOMOCOv11_content = file.read()

def check_motif_contains_CG(motif_id, protein_name):
    if not isinstance(motif_id, str):
        print(f"motif is empty: {motif_id}")
        return ""
    if motif_id.isupper() and motif_id.isalpha():
        file_path = os.path.join(base_path, f"pooled_{protein_name}/meme_out/meme.txt")
    elif re.match(r"^\d+-[A-Z]+", motif_id):
        file_path = os.path.join(base_path, f"pooled_{protein_name}/streme_out/streme.txt")
    else:
        file_path = HOCOMOCOv11_path

    with open(file_path, 'r') as file:
        content = file.read()

    # Search for the motif section
    motif_section_pattern = rf"MOTIF {motif_id}.*?letter-probability matrix:(.*?)(?=\nMOTIF|\n\n|\n-|\Z|\nURL)"
    match = re.search(motif_section_pattern, content, re.DOTALL)

    if match:
        #print(f"Matched section for motif {motif_id}:")
        # print(match.group(0))  # Print the entire matched section

        matrix_content = match.group(1).strip().split('\n', 1)[1]  # Remove the first line
        rows = matrix_content.strip().split('\n')
        # print(f"Letter-probability matrix for motif {motif_id}:")
        # print(matrix_content)  # Print the specific capture group (the matrix)
        probabilities = np.array([list(map(float, row.split())) for row in rows])

        for i in range(len(probabilities) - 1):
            if probabilities[i, 1] == max(probabilities[i, :]) and probabilities[i + 1, 2] == max(
                    probabilities[i + 1, :]):
                return True
        return False
    else:
        print(f"Motif section for {motif_id} not found in {file_path}.")
        return ""


def check_CG_pattern(row):
    return check_motif_contains_CG(row['single Motif ID'], row['Protein'])

# Load the data
all_groups_df = pd.read_csv(
    "C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/meme/meme-chip_results"
    "/all_motifs_groups_all.csv")

# Add the 'Contains CG' column
all_groups_df.loc[:, 'Contains CG'] = all_groups_df.apply(check_CG_pattern, axis=1)
# Save the DataFrame to a new CSV file
all_groups_df.to_csv(
    "C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/meme/meme-chip_results"
    "/all_motifs_groups_all_with_CG.csv",
    index=False)

print(
    "All motifs with CG information saved to: C:/Users/Daniel "
    "Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/meme/meme-chip_results"
    "/all_motifs_groups_all_with_CG.csv")
