import pandas as pd
import matplotlib.pyplot as plt

# Load the data
df = pd.read_csv("C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/meme/meme-chip_results/all_motifs_groups_all_with_CG.csv")

# Convert 'Distribution' column to numeric, setting errors='coerce' will turn non-numeric values to NaN
df['Distribution'] = pd.to_numeric(df['Distribution'], errors='coerce')

# Replace NaN with 0
df['Distribution'].fillna(0, inplace=True)

# Exclude rows with concentration values equal to 0
df_nonzero = df[df['Distribution'] > 0]

# Plot the distribution of all non-zero concentration values
plt.figure(figsize=(12, 6))
plt.hist(df_nonzero['Distribution'], bins=int(0.1/0.001), range=(0, 0.1), edgecolor='black')
plt.xlabel('Concentration')
plt.ylabel('Frequency')
plt.title('Distribution of All Non-Zero Concentration Values')
plt.xticks(ticks=[x * 0.01 for x in range(0, 11)], rotation=90)
plt.xlim(0, 0.1)  # Focus on the range 0 to 0.1
plt.grid(True)
plt.show()

# Select the highest concentration value per protein
highest_concentration_per_protein = df_nonzero.loc[df_nonzero.groupby('Protein')['Distribution'].idxmax()]

# Plot the distribution of the highest concentration value per protein
plt.figure(figsize=(12, 6))
plt.hist(highest_concentration_per_protein['Distribution'], bins=int(0.1/0.01), range=(0, 0.2), edgecolor='black')
plt.xlabel('Concentration')
plt.ylabel('Frequency')
plt.title('Distribution of Highest Concentration Value per Protein')
plt.xticks(ticks=[x * 0.01 for x in range(0, 21)], rotation=90)
plt.xlim(0, 0.2)
plt.grid(True)
plt.show()

# Filter the DataFrame to include only rows with concentration values >= 0.05
df_filtered = df_nonzero[df_nonzero['Distribution'] >= 0.05]

# Create a set to store unique motif occurrences per protein
unique_motifs_per_protein = set()

# Process the "Known or Similar Motifs" column for filtered data
for _, row in df_filtered.iterrows():
    protein = row['Protein']
    known_similar_motifs = row['Known or Similar Motifs']
    if pd.notna(known_similar_motifs):  # Check if the cell is not NaN
        motifs = known_similar_motifs.split('; ')
        for motif in motifs:
            unique_motifs_per_protein.add((protein, motif))

# Convert the set to a DataFrame
unique_motifs_df = pd.DataFrame(list(unique_motifs_per_protein), columns=['Protein', 'Motif'])

# Count the occurrences of each motif
motif_counts = unique_motifs_df['Motif'].value_counts()

# Filter to include only the top 10 motifs for better visualization
top_motifs = motif_counts.head(10)
other_motifs_count = motif_counts[10:].sum()

# Add the count for other motifs as 'Other'
top_motifs['Other'] = other_motifs_count

# Create a pie chart
plt.figure(figsize=(10, 7))
plt.pie(top_motifs, labels=top_motifs.index, autopct='%1.1f%%', startangle=140)
plt.title('Distribution of Known or Similar Motifs (Concentration >= 0.05)')
plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.show()

print(motif_counts.head(10))

# Count the number of unique "single Motif ID" entries per protein
motif_count_per_protein = df_filtered.groupby('Protein')['single Motif ID'].nunique().reset_index()

# Rename the columns for clarity
motif_count_per_protein.columns = ['Protein', 'Unique Motif ID Count']

# Display the DataFrame
print(motif_count_per_protein)

# Count the number of different groups per protein that have at least one entry with concentration values above 0.05
group_count_per_protein = df_filtered.groupby('Protein')['Group ID'].nunique().reset_index()

# Rename the columns for clarity
group_count_per_protein.columns = ['Protein', 'Group Count']

# Plot the number of different groups per protein with concentration values above 0.05
plt.figure(figsize=(12, 6))
plt.bar(group_count_per_protein['Protein'], group_count_per_protein['Group Count'], color='skyblue')
plt.xlabel('Protein')
plt.ylabel('Number of Different Groups')
plt.title('Number of Different Groups per Protein (Concentration >= 0.05)')
plt.xticks(rotation=90)
plt.show()

# Display the DataFrame
print(group_count_per_protein)
