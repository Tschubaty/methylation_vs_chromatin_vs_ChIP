import pandas as pd
import re
import numpy as np


def check_motif_contains_CG_in_streme(motif_id, motif_source_file):
    try:
        motif_source_file = motif_source_file.replace(
            '/ems/elsc-labs/meshorer-e/daniel.batyrev/Encode3/meme/meme-chip_results/',
            'C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/meme/meme-chip_results/').replace(
            'xml', 'txt')
        with open(motif_source_file, 'r') as file:
            content = file.read()

        # print(f"Content of {motif_source_file} (first 500 characters):")
        # print(content[:500])

        # Search for the motif section
        motif_section_pattern = rf"MOTIF {motif_id}.*?letter-probability matrix:(.*?)(?=\nMOTIF|\Z)"
        match = re.search(motif_section_pattern, content, re.DOTALL)

        if match:
            print(f"Matched section for motif {motif_id}:")
            print(match.group(0))  # Print the entire matched section

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
            print(f"Motif section for {motif_id} not found in {motif_source_file}.")
            return ""

    except Exception as e:
        print(f"Error checking motif {motif_id} in file {motif_source_file}: {e}")
        return ""


def check_CG_pattern(row):
    return check_motif_contains_CG_in_streme(row['FIMO Motif ID'], row['FIMO Motif source file'])


# Load the data
all_groups_df = pd.read_csv(
    "C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/meme/meme-chip_results/all_motifs_groups3.csv")

# Filter the DataFrame to include only rows with FIMO files and the highest concentration value
filtered_df = all_groups_df.dropna(subset=['FIMO gff']).copy()
filtered_df['Distribution'] = pd.to_numeric(filtered_df['Distribution'], errors='coerce')
filtered_df = filtered_df.loc[filtered_df.groupby('Protein')['Distribution'].idxmax()]

# Add the 'Contains CG' column
filtered_df.loc[:, 'Contains CG'] = filtered_df.apply(check_CG_pattern, axis=1)
print(filtered_df.loc[:, 'Contains CG'])
# Save the filtered DataFrame to a new CSV file
filtered_df.to_csv(
    "C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/meme/meme-chip_results/best_motifs_with_CG.csv",
    index=False)

print(
    "Best motifs with CG information saved to: C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/meme/meme-chip_results/best_motifs_with_CG.csv")
