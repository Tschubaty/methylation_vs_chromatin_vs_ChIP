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
all_motifs = []

# Function to extract concentration values from the Centrimo file
def extract_concentration_from_centrimo(centrimo_file, alt_id):
    driver.get(centrimo_file)
    time.sleep(3)  # Adjust if necessary
    rows = driver.find_elements(By.CSS_SELECTOR, "#motifs tbody tr")
    for row in rows:
        cols = row.find_elements(By.TAG_NAME, "td")
        if len(cols) > 5:
            if cols[3].text.strip() == alt_id:
                return float(cols[5].text.strip())
    return None

# Loop over all folders (proteins)
for protein_folder in os.listdir(base_path):
    if protein_folder.startswith('pooled_'):
        protein_name = protein_folder[len('pooled_'):]
        html_file_path = f'file:///{base_path}/{protein_folder}/meme-chip.html'

        # Load the HTML file
        print(f"Loading HTML file for {protein_name}...")
        driver.get(html_file_path)

        # Allow time for the JavaScript to execute
        time.sleep(3)  # Adjust the sleep time if necessary

        # Extract the group information
        motif_boxes = driver.find_elements(By.CSS_SELECTOR, ".motifbox")
        print(f"Found {len(motif_boxes)} motif boxes in {protein_name}.")

        for idx, box in enumerate(motif_boxes, start=1):
            print(f"Processing box {idx} of {len(motif_boxes)}")
            try:
                group = box.find_element(By.CSS_SELECTOR, ".motifs")
                rows = group.find_elements(By.TAG_NAME, "tr")
                group_id = f"Group_{idx}"
                group_name = "CentriMo Group ↷" if "CentriMo Group ↷" in box.text else "Other Group"

                for row in rows:
                    cols = row.find_elements(By.TAG_NAME, "td")
                    if len(cols) > 5:
                        motif_id = cols[0].text.strip()
                        motif_found = cols[0].text.strip()
                        program = cols[1].text.strip()
                        e_value = cols[2].text.strip()
                        known_motifs = cols[3].text.strip().replace('\n', '; ')
                        distribution = cols[4].text.strip()
                        spamo_fimo = cols[5].text.strip().replace('\n', '; ')

                        centrimo_link = f"{base_path}/{protein_folder}/centrimo_out/centrimo.html?show=db1+{motif_found}"
                        spamo_fimo_link = f"{base_path}/{protein_folder}/spamo_out_{idx}/spamo.html"

                        motif_data = {
                            "Group Name": group_name,
                            "Group ID": group_id,
                            "Protein": protein_name,
                            "Motif ID": motif_id,
                            "Motif Found": motif_found,
                            "Discovery/Enrichment Program": program,
                            "E-value": e_value,
                            "Known or Similar Motifs": known_motifs,
                            "Distribution": distribution,
                            "SpaMo & FIMO": spamo_fimo_link,
                            "CentriMo Group": "Yes" if "CentriMo Group" in group_name else "No",
                            "Centrimo Link": centrimo_link
                        }

                        all_motifs.append(motif_data)

            except Exception as e:
                print(f"Error processing box: {e}")

# Extract concentration values for each motif
for motif in all_motifs:
    alt_id = motif["Motif Found"]
    centrimo_link = motif["Centrimo Link"]
    concentration_value = extract_concentration_from_centrimo(centrimo_link, alt_id)
    motif["Concentration"] = concentration_value

# Close the WebDriver
driver.quit()

# Create a DataFrame for all motifs and save to a CSV file
motifs_df = pd.DataFrame(all_motifs)
output_csv_path = "C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/meme/meme-chip_results/all_motifs_groups.csv"
motifs_df.to_csv(output_csv_path, index=False)

print(f"All motifs saved to: {output_csv_path}")
