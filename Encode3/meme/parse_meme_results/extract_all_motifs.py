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

# Prepare to store the motif group information for all proteins
all_groups = []
group_counter = 1  # Counter to assign unique IDs to groups

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

        # Check if the HTML file is loaded
        if "404 Not Found" in driver.title:
            print(f"Error: {html_file_path} not found")
            continue

        # Extract the groups and motifs
        print("Extracting groups and motifs...")
        motif_boxes = driver.find_elements(By.CSS_SELECTOR, ".motifbox")
        print(f"Found {len(motif_boxes)} motif boxes in {protein_name}.")

        for box in motif_boxes:
            try:
                # Check the inner HTML of the motif box for debugging
                inner_html = box.get_attribute('innerHTML')
                print(f"Motif box inner HTML: {inner_html[:200]}...")  # Print the first 200 characters

                # Extract group information
                group_table = box.find_element(By.CSS_SELECTOR, "table.motifs")
                group_name = box.find_element(By.CSS_SELECTOR, "span.action a").text.strip()
                group_id = f"Group_{group_counter}"
                group_counter += 1
                print(f"Processing group: {group_name} with ID {group_id}")
                rows = group_table.find_elements(By.CSS_SELECTOR, "tbody tr")
                print(f"Found {len(rows)} rows in group {group_name}.")

                # Extract motif information for each row
                for row in rows:
                    cols = row.find_elements(By.TAG_NAME, "td")
                    print(f"Found {len(cols)} columns in row.")
                    if len(cols) > 5:
                        # Extract the link from the "Discovery/Enrichment Program" column
                        link_element = cols[1].find_element(By.TAG_NAME, "a")
                        link_url = link_element.get_attribute("href") if link_element else ""

                        # Extract all links from the "SpaMo & FIMO" column
                        spamo_links = cols[5].find_elements(By.TAG_NAME, "a")
                        spamo_link_urls = [link.get_attribute("href") for link in spamo_links]
                        spamo_links_combined = "; ".join(spamo_link_urls)

                        motif_data = {
                            "Group Name": group_name,
                            "Group ID": group_id,
                            "Protein": protein_name,
                            "Motif ID": cols[0].text.strip() if len(cols) > 0 else "",  # Add Motif ID
                            "Motif Found": cols[1].text.strip() if len(cols) > 1 else "",
                            "Discovery/Enrichment Program": link_url,
                            "E-value": cols[2].text.strip() if len(cols) > 2 else "",
                            "Known or Similar Motifs": cols[3].text.strip().replace("\n", "; ") if len(cols) > 3 else "",
                            "Distribution": cols[4].text.strip() if len(cols) > 4 else "",
                            "SpaMo & FIMO": spamo_links_combined,
                            "CentriMo Group": "Yes" if "CentriMo Group" in group_name else "No"
                        }
                        all_groups.append(motif_data)
                        print(f"Added motif data: {motif_data}")
            except Exception as e:
                print(f"Error processing box: {e}")
                continue

        print(f"Extracted motifs for {protein_name}.")

# Close the WebDriver
driver.quit()

# Create a DataFrame for all motifs and save to a CSV file
all_groups_df = pd.DataFrame(all_groups)
output_csv_path = "C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/meme/meme-chip_results/all_motifs_groups3.csv"
all_groups_df.to_csv(output_csv_path, index=False)

print(f"All motifs saved to: {output_csv_path}")
