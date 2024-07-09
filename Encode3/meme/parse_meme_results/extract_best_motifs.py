import pandas as pd
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
import os
import time

# Path to your local ChromeDriver
chromedriver_path = 'C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/meme/chromedriver-win64/chromedriver.exe'

# Base path to your meme-chip results
base_path = 'C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/meme/meme-chip_results/'

# Initialize the Chrome WebDriver
service = Service(chromedriver_path)
driver = webdriver.Chrome(service=service)

# Prepare to store the best rows for all proteins
best_rows = []

# Loop over all folders (proteins)
for protein_folder in os.listdir(base_path):
    if protein_folder.startswith('pooled_'):
        protein_name = protein_folder[len('pooled_'):]
        html_file_path = f'file:///{base_path}/{protein_folder}/centrimo_out/centrimo.html'

        # Load the HTML file
        print(f"Loading HTML file for {protein_name}...")
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

        print(f"Extracted {len(extracted_data)} rows for {protein_name}.")
        if extracted_data:
            print(f"Sample data row ({len(extracted_data[0])}): {extracted_data[0]}")

        # Adjust the columns if necessary
        if len(columns) != len(extracted_data[0]):
            columns = columns[:len(extracted_data[0])]
            print(f"Adjusted column names ({len(columns)}): {columns}")

        # Convert the extracted data to a DataFrame
        extracted_df = pd.DataFrame(extracted_data, columns=columns)

        # Replace empty strings with NaN
        extracted_df.replace("", float("NaN"), inplace=True)

        # Remove empty columns
        extracted_df.dropna(axis=1, how='all', inplace=True)

        # Display the first few rows of the DataFrame
        print("Extracted DataFrame after removing empty columns:")
        print(extracted_df.head())

        # Convert E-value column to float and sort by E-value
        extracted_df['E-value'] = extracted_df['E-value'].astype(float)
        sorted_df = extracted_df.sort_values(by='E-value')

        # Limit to top best E-values (e.g., top 10)
        top_best_evalues = 10
        top_df = sorted_df.head(top_best_evalues)

        # Find the best concentration value among the top best E-values
        best_row = top_df.loc[top_df['Concentration'].astype(float).idxmax()].copy()

        print(f"Best row with the highest concentration value for {protein_name}:")
        print(best_row)

        # Add protein name to the best row
        best_row['Protein'] = protein_name
        best_rows.append(best_row)

# Close the WebDriver
driver.quit()

# Create a DataFrame for all best rows and save to a CSV file
best_rows_df = pd.DataFrame(best_rows)
output_csv_path = "C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/meme/meme-chip_results/best_concentration_values.csv"
best_rows_df.to_csv(output_csv_path, index=False)

print(f"Best concentration values saved to: {output_csv_path}")
