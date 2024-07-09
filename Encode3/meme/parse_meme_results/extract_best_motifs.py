import pandas as pd
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
import time

# Path to your local ChromeDriver
chromedriver_path = 'C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/meme/chromedriver-win64/chromedriver.exe'

# Path to your local HTML file
html_file_path = 'file:///C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/meme/meme-chip_results/pooled_AGO1/centrimo_out/centrimo.html'

# Initialize the Chrome WebDriver
service = Service(chromedriver_path)
driver = webdriver.Chrome(service=service)

# Load the HTML file
print("Loading HTML file...")
driver.get(html_file_path)

# Allow time for the JavaScript to execute
time.sleep(3)  # Adjust the sleep time if necessary

# Extract the column names
print("Extracting column names...")
header = driver.find_element(By.CSS_SELECTOR, "#motifs thead tr")
columns = [th.text.strip() for th in header.find_elements(By.TAG_NAME, "th")]

print(f"Extracted column names ({len(columns)}): {columns}")

# Extract the data rows
print("Extracting data rows...")
rows = driver.find_elements(By.CSS_SELECTOR, "#motifs tbody tr")

extracted_data = []
for row in rows:
    cols = row.find_elements(By.TAG_NAME, "td")
    if len(cols) > 5:
        data_row = [col.text.strip() for col in cols]
        extracted_data.append(data_row)

print(f"Extracted {len(extracted_data)} rows.")
if extracted_data:
    print(f"Sample data row ({len(extracted_data[0])}): {extracted_data[0]}")


# Close the WebDriver
driver.quit()

# Adjust the columns if necessary
if len(columns) != len(extracted_data[0]):
    columns = columns[:len(extracted_data[0])]
    print(f"Adjusted column names ({len(columns)}): {columns}")

# Convert the extracted data to a DataFrame
extracted_df = pd.DataFrame(extracted_data, columns=columns)

# Remove empty columns
extracted_df.dropna(axis=1, how='all', inplace=True)

# Display the first few rows of the DataFrame
print("Extracted DataFrame:")
print(extracted_df.head())

# Find the best concentration value
best_row = extracted_df.loc[extracted_df['Concentration'].astype(float).idxmax()]

print("Best row with the highest concentration value:")
print(best_row)

print(best_row)

# Load the TSV file
tsv_file_path = "C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/UniProt/idmapping_reviewed_true_AND_model_organ_2024_06_16.tsv"
print(f"Loading TSV file from {tsv_file_path}...")
tsv_df = pd.read_csv(tsv_file_path, sep='\t')

print("Original TSV DataFrame:")
print(tsv_df.head())

# Append the best row's data to the TSV DataFrame
for col in column_names:
    tsv_df[col] = best_row[col]

print("Updated TSV DataFrame with new columns:")
print(tsv_df.head())

# Save the updated DataFrame to a new TSV file
output_tsv_path = "C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/UniProt/updated_idmapping.tsv"
tsv_df.to_csv(output_tsv_path, sep='\t', index=False)

print(f"Updated TSV file saved to: {output_tsv_path}")
