import matplotlib.pyplot as plt
import os
import pandas as pd

# Define the base directory containing the results
base_dir = r"D:\Users\Daniel Batyrev\Documents\GitHub\homer\appended_bed"

# Initialize a list to store the best results
best_results = []

# Iterate over each directory in the base directory
for protein_name in os.listdir(base_dir):
    dir_path = os.path.join(base_dir, protein_name)

    # Check if it's a directory
    if os.path.isdir(dir_path):
        # Check if knownResults.txt exists in the directory
        file_path = os.path.join(dir_path, "knownResults.txt")

        if not os.path.isfile(file_path):
            print(f"Error: knownResults.txt not found in {dir_path}. Appending NaN values.")
            # Append NaN values if the file is not found
            best_results.append({
                'Folder': protein_name,
                'Protein': protein_name,
                'Motif Name': None,
                'Short Motif Name': None,
                'Consensus': None,
                'P-value': None,
                'Log P-value': None,
                'q-value (Benjamini)': None,
                '# of Target Sequences with Motif': None,
                'percent of Target Sequences with Motif': None,
                '# of Background Sequences with Motif': None,
                'percent of Background Sequences with Motif': None
            })
            continue

        # Define the hardcoded column names
        column_names = [
            "Motif Name",
            "Consensus",
            "P-value",
            "Log P-value",
            "q-value (Benjamini)",
            "# of Target Sequences with Motif",
            "percent of Target Sequences with Motif",
            "# of Background Sequences with Motif",
            "percent of Background Sequences with Motif"
        ]

        # Read the data into a pandas DataFrame, skipping the first row for headers
        df = pd.read_csv(file_path, sep="\t", header=None, names=column_names, skiprows=1)

        # Remove the '%' sign and convert the percentage columns to float
        df['percent of Target Sequences with Motif'] = df['percent of Target Sequences with Motif'].str.replace('%',
                                                                                                                '').astype(
            float)
        df['percent of Background Sequences with Motif'] = df['percent of Background Sequences with Motif'].str.replace(
            '%', '').astype(float)

        # Convert 'P-value' and 'q-value (Benjamini)' to numeric, coercing errors to NaN
        df['P-value'] = pd.to_numeric(df['P-value'], errors='coerce')
        df['q-value (Benjamini)'] = pd.to_numeric(df['q-value (Benjamini)'], errors='coerce')

        # Drop rows where 'P-value' is NaN
        df_clean = df.dropna(subset=['P-value'])

        # Filter the DataFrame for motifs with P-value < 0.00005
        significant_motifs = df_clean[df_clean['P-value'] < 0.00005]

        # Check if there are any significant motifs
        if not significant_motifs.empty:
            for _, row in significant_motifs.iterrows():
                # Extract short motif name
                short_motif_name = row['Motif Name'].split('/')[0]

                # Add the motif details to the results list
                best_results.append({
                    'Folder': protein_name,
                    'Protein': protein_name,
                    'Motif Name': row['Motif Name'],
                    'Short Motif Name': short_motif_name,
                    'Consensus': row['Consensus'],
                    'P-value': row['P-value'],
                    'Log P-value': row['Log P-value'],
                    'q-value (Benjamini)': row['q-value (Benjamini)'],
                    '# of Target Sequences with Motif': row['# of Target Sequences with Motif'],
                    'percent of Target Sequences with Motif': row['percent of Target Sequences with Motif'],
                    '# of Background Sequences with Motif': row['# of Background Sequences with Motif'],
                    'percent of Background Sequences with Motif': row['percent of Background Sequences with Motif']
                })
        else:
            print(f"Warning: No motifs with P-value < 0.00005 found in {file_path}. Appending NaN values.")
            best_results.append({
                'Folder': protein_name,
                'Protein': protein_name,
                'Motif Name': None,
                'Short Motif Name': None,
                'Consensus': None,
                'P-value': None,
                'Log P-value': None,
                'q-value (Benjamini)': None,
                '# of Target Sequences with Motif': None,
                'percent of Target Sequences with Motif': None,
                '# of Background Sequences with Motif': None,
                'percent of Background Sequences with Motif': None
            })

# Convert the list of best results into a pandas DataFrame
best_results_df = pd.DataFrame(best_results)

# Print the combined dataframe
print(best_results_df)

# Define the output path
output_dir = r"C:\Users\Daniel Batyrev\Documents\GitHub\methylation_vs_chromatin_vs_ChIP\Encode3\homer\parse_results"
os.makedirs(output_dir, exist_ok=True)  # Ensure the directory exists
output_path = os.path.join(output_dir, "best_motifs_results.csv")

# Save the results to a CSV file
best_results_df.to_csv(output_path, index=False)

print(f"Results saved to {output_path}")

# Count the occurrences of each Motif Name
motif_counts = best_results_df['Motif Name'].value_counts().head(50)

# Plot the 50 most occurring Motif Names
plt.figure(figsize=(14, 10))
motif_counts.plot(kind='bar')
plt.title('Top 50 Most Occurring Motif Names')
plt.xlabel('Motif Name')
plt.ylabel('Frequency')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()

# Save the plot as an image file (optional)
plt.savefig(os.path.join(output_dir, 'motif_name_occurrences.png'))

# Display the plot
plt.show()

# Plot the distribution of q-value (Benjamini)
plt.figure(figsize=(12, 8))
plt.hist(best_results_df['q-value (Benjamini)'].dropna(), bins=50, edgecolor='black', alpha=0.7)
plt.title('Distribution of q-value (Benjamini)')
plt.xlabel('q-value (Benjamini)')
plt.ylabel('Frequency')
plt.tight_layout()

# Save the plot as an image file (optional)
plt.savefig(os.path.join(output_dir, 'q_value_benjamini_distribution.png'))

# Display the plot
plt.show()
