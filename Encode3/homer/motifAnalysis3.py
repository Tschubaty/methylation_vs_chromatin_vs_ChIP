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
                'Short Motif Name': None,  # Add this field
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

        # Convert 'P-value' to numeric, coercing errors to NaN
        df['P-value'] = pd.to_numeric(df['P-value'], errors='coerce')

        # Drop rows where 'P-value' is NaN
        df_clean = df.dropna(subset=['P-value'])

        # Check if df_clean is not empty
        if not df_clean.empty:
            # Extract the best motif based on the lowest P-value
            best_motif = df_clean.loc[df_clean['P-value'].idxmin()]

            # Extract short motif name
            short_motif_name = best_motif['Motif Name'].split('/')[0]

            # Add the folder name, protein name, and short motif name to the dataframe
            best_results.append({
                'Folder': protein_name,
                'Protein': protein_name,
                'Motif Name': best_motif['Motif Name'],
                'Short Motif Name': short_motif_name,  # Add short motif name
                'Consensus': best_motif['Consensus'],
                'P-value': best_motif['P-value'],
                'Log P-value': best_motif['Log P-value'],
                'q-value (Benjamini)': best_motif['q-value (Benjamini)'],
                '# of Target Sequences with Motif': best_motif['# of Target Sequences with Motif'],
                'percent of Target Sequences with Motif': best_motif['percent of Target Sequences with Motif'],
                '# of Background Sequences with Motif': best_motif['# of Background Sequences with Motif'],
                'percent of Background Sequences with Motif': best_motif['percent of Background Sequences with Motif']
            })
        else:
            print(f"Warning: No valid P-values found in {file_path}. Appending NaN values.")
            best_results.append({
                'Folder': protein_name,
                'Protein': protein_name,
                'Motif Name': None,
                'Short Motif Name': None,  # Add this field
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

# Count the occurrences of each Short Motif Name
short_motif_counts = best_results_df['Short Motif Name'].value_counts().head(50)

# Plot the 50 most occurring Short Motif Names
plt.figure(figsize=(12, 8))
short_motif_counts.plot(kind='bar')
plt.title('Top 50 Most Occurring Short Motif Names')
plt.xlabel('Short Motif Name')
plt.ylabel('Frequency')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()

# Save the plot as an image file (optional)
plt.savefig(os.path.join(output_dir, 'short_motif_name_occurrences.png'))

# Display the plot
plt.show()
