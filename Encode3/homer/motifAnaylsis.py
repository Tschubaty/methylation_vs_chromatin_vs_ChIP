import matplotlib.pyplot as plt
import os
import pandas as pd

# Define the base directory containing the results
base_dir = r"D:\Users\Daniel Batyrev\Documents\GitHub\homer\appended_bed"

# Initialize a list to store the best results
best_results = []

# Iterate over each subdirectory in the base directory
for root, dirs, files in os.walk(base_dir):
    for file in files:
        if file == "knownResults.txt":
            file_path = os.path.join(root, file)

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
            df['percent of Target Sequences with Motif'] = df['percent of Target Sequences with Motif'].str.replace('%', '').astype(float)
            df['percent of Background Sequences with Motif'] = df['percent of Background Sequences with Motif'].str.replace('%', '').astype(float)

            # Convert 'P-value' to numeric, coercing errors to NaN
            df['P-value'] = pd.to_numeric(df['P-value'], errors='coerce')

            # Drop rows where 'P-value' is NaN
            df_clean = df.dropna(subset=['P-value'])

            # Check if df_clean is not empty
            if not df_clean.empty:
                # Extract the best motif based on the lowest P-value
                best_motif = df_clean.loc[df_clean['P-value'].idxmin()]

                # Extract protein name from the directory path
                protein_name = os.path.basename(root)

                # Add the folder name and protein name to the dataframe
                best_results.append({'Folder': protein_name, 'Protein': protein_name, **best_motif})

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
plt.figure(figsize=(12, 8))
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
