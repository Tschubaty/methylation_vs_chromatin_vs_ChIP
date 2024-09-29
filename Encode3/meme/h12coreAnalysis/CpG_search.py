import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles

# Function to parse the MEME file and identify motifs with CpG sites
def parse_meme_file(meme_file):
    motifs = []
    motif = None
    matrix_started = False

    with open(meme_file, 'r') as file:
        for line in file:
            if line.startswith('MOTIF'):
                if motif:
                    motifs.append(motif)
                motif = {'header': line.strip(), 'matrix': []}
                matrix_started = False
            elif line.startswith('letter-probability matrix'):
                matrix_started = True
            elif matrix_started and line.strip() and not line.startswith(('URL', 'MOTIF')):
                motif['matrix'].append([float(x) for x in line.strip().split()])

        if motif:
            motifs.append(motif)

    return motifs

def find_motifs_with_CpG(motifs):
    motifs_with_CpG = []
    for motif in motifs:
        matrix = motif['matrix']
        for i in range(len(matrix) - 1):
            c_prob = matrix[i][1]  # Probability of 'C'
            g_prob = matrix[i + 1][2]  # Probability of 'G'
            # Define a threshold for considering a CpG site
            if c_prob > 0.5 and g_prob > 0.5:
                motifs_with_CpG.append(motif['header'])
                break
    return motifs_with_CpG

# Step 1: Parse the MEME file and find motifs with CpG
meme_file_path = r'C:\Users\Daniel Batyrev\Documents\GitHub\methylation_vs_chromatin_vs_ChIP\Encode3\meme\H12CORE_meme_format.meme'
motifs = parse_meme_file(meme_file_path)
motifs_with_CpG = find_motifs_with_CpG(motifs)

# Step 2: Adjust the motif names by removing the 'MOTIF ' prefix
adjusted_motifs_with_CpG = {motif.split(' ')[1] for motif in motifs_with_CpG}

# Step 3: Load the provided TSV file (H12CORE motifs)
motif_tsv_path = r'C:\Users\Daniel Batyrev\Documents\GitHub\methylation_vs_chromatin_vs_ChIP\Encode3\meme\H12CORE_motifs.tsv'
motifs_df = pd.read_csv(motif_tsv_path, sep='\t')

# Step 4: Add a new column 'Contains_CpG' based on whether the motif is in the CpG list
motifs_df['Contains_CpG'] = motifs_df['Motif'].apply(lambda x: 'Yes' if x in adjusted_motifs_with_CpG else 'No')

# Step 5: Save the updated DataFrame back to a TSV file
updated_motif_tsv_path = r'C:\Users\Daniel Batyrev\Documents\GitHub\methylation_vs_chromatin_vs_ChIP\Encode3\meme\H12CORE_motifs_with_CpG.tsv'
motifs_df.to_csv(updated_motif_tsv_path, sep='\t', index=False)

# Step 6: Load the UniProt mapping file
idmapping_file_path = r'C:\Users\Daniel Batyrev\Documents\GitHub\methylation_vs_chromatin_vs_ChIP\Encode3\UniProt\idmapping_with_H12CORE_motifs.tsv'
idmapping_df = pd.read_csv(idmapping_file_path, sep='\t')

# Step 7: Append a new column 'Contains_CpG' based on the 'H12CORE_motif' column
idmapping_df['Contains_CpG'] = idmapping_df['H12CORE_motif'].apply(
    lambda x: 'Yes' if x in adjusted_motifs_with_CpG else 'No')

# Step 8: Save the updated UniProt mapping file with CpG information
updated_idmapping_tsv_path = r'C:\Users\Daniel Batyrev\Documents\GitHub\methylation_vs_chromatin_vs_ChIP\Encode3\UniProt\idmapping_with_H12CORE_motifs_with_CpG.tsv'
idmapping_df.to_csv(updated_idmapping_tsv_path, sep='\t', index=False)

# Step 9: Define the sets
# Set 1: All proteins with motifs in H12CORE_motifs
proteins_with_motifs = set(motifs_df['Gene (human)'])

# Subset of Set 1: Proteins with CpG-containing motifs in H12CORE_motifs
proteins_with_CpG_motifs = set(motifs_df[motifs_df['Contains_CpG'] == 'Yes']['Gene (human)'])

# Set 2: All proteins in idmapping
proteins_in_idmapping = set(idmapping_df['From'])

# Subset of Set 2: Proteins in idmapping with H12CORE motifs
proteins_with_idmapping_motifs = set(idmapping_df[idmapping_df['H12CORE_motif'].notnull()]['From'])

# Subset of Set 2: Proteins in idmapping with H12CORE motifs and CpG
proteins_with_idmapping_CpG = set(idmapping_df[idmapping_df['Contains_CpG'] == 'Yes']['From'])

# Step 10: Create the Venn diagram for the two primary sets
plt.figure(figsize=(8, 8))

# Prepare the Venn diagram with two primary sets:
venn = venn2(
    [proteins_with_motifs, proteins_in_idmapping],
    (f"Proteins with H12CORE Motifs:\n {len(proteins_with_motifs)}", f"ChIP Data in min 2 samples:\n {len(proteins_in_idmapping)}")
)

# Optionally, annotate subsets within the diagram
# Example: Number of proteins with CpG motifs within H12CORE motifs
if venn.get_label_by_id('10'):  # Only annotate if the subset exists
    venn.get_label_by_id('10').set_text(f"{len(proteins_with_motifs - proteins_with_idmapping_motifs)}")

if venn.get_label_by_id('01'):  # Only annotate if the subset exists
    venn.get_label_by_id('01').set_text(f"{len(proteins_in_idmapping - proteins_with_idmapping_motifs)}")

if venn.get_label_by_id('11'):  # Only annotate if the subset exists
    venn.get_label_by_id('11').set_text(f"{len(proteins_with_idmapping_motifs)}\n({len(proteins_with_idmapping_CpG)} with CpG)")

# Optionally customize the diagram with circles
venn_circles = venn2_circles(
    [proteins_with_motifs, proteins_in_idmapping],
    linestyle='dashed'
)

# Display the plot
plt.title("Venn Diagram Motif Data Base vs Data available")
plt.show()

# Function to calculate the CpG count and total count for each TF family
def calculate_cpg_counts(df):
    total_count = len(df)
    cpg_count = df['Contains_CpG'].value_counts().get('Yes', 0)
    return pd.Series({'TF_Family_CpG_Count': cpg_count, 'TF_Family_Total_Count': total_count})

# Apply the function to each TF family and create new columns
cpg_counts = motifs_df.groupby('TF family').apply(calculate_cpg_counts).reset_index()

# Merge the CpG counts back into the original DataFrame
motifs_df = motifs_df.merge(cpg_counts, on='TF family', how='left')

# Function to calculate the CpG count and total count for each TF subfamily
def calculate_subfamily_cpg_counts(df):
    total_count = len(df)
    cpg_count = df['Contains_CpG'].value_counts().get('Yes', 0)
    return pd.Series({'TF_Subfamily_CpG_Count': cpg_count, 'TF_Subfamily_Total_Count': total_count})

# Apply the function to each TF subfamily and create new columns
subfamily_cpg_counts = motifs_df.groupby(['TF family', 'TF subfamily']).apply(calculate_subfamily_cpg_counts).reset_index()

# Merge the CpG subfamily counts back into the original DataFrame
motifs_df = motifs_df.merge(subfamily_cpg_counts, on=['TF family', 'TF subfamily'], how='left')

# Save the updated DataFrame back to a TSV file with the new columns
updated_motif_tsv_path = r'C:\Users\Daniel Batyrev\Documents\GitHub\methylation_vs_chromatin_vs_ChIP\Encode3\meme\H12CORE_motifs_with_Family_and_Subfamily_CpG_Counts.tsv'
motifs_df.to_csv(updated_motif_tsv_path, sep='\t', index=False)

# Display the DataFrame to verify the new columns
print(motifs_df.head())

# Aggregate the data to get unique family and counts
family_counts = motifs_df[['TF family', 'TF_Family_CpG_Count', 'TF_Family_Total_Count']].drop_duplicates()

# Sort families by total count for better visualization
family_counts = family_counts.sort_values(by='TF_Family_Total_Count', ascending=False)

# Create the bar plot
plt.figure(figsize=(12, 8))

# Plot bars for CpG-containing motifs
plt.bar(family_counts['TF family'], family_counts['TF_Family_CpG_Count'], color='blue', label='Motifs with CpG')

# Plot bars for motifs without CpG on top of the previous bars
plt.bar(family_counts['TF family'], family_counts['TF_Family_Total_Count'] - family_counts['TF_Family_CpG_Count'],
        bottom=family_counts['TF_Family_CpG_Count'], color='orange', label='Motifs without CpG')

# Add labels and title
plt.xlabel('TF Family')
plt.ylabel('Number of Motifs')
plt.title('Distribution of Motifs with and without CpG by TF Family')
plt.xticks(rotation=90)
plt.legend()

# Display the plot
plt.tight_layout()
plt.show()