import os
import re
import time
import pandas as pd
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By

# Path to your local ChromeDriver
chromedriver_path = 'C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/meme/chromedriver-win64/chromedriver.exe'

# Base path to your meme-chip results
base_path = r"D:/Users/Daniel Batyrev/Documents/GitHub/meme/single_zip"

# Initialize the Chrome WebDriver
service = Service(chromedriver_path)
driver = webdriver.Chrome(service=service)

# Prepare to store the motif group information for all proteins
all_groups = []
group_counter = 1  # Counter to assign unique IDs to groups

# Loop over all folders (single experiment folders)
for protein_folder in os.listdir(base_path):
    start_group = group_counter
    found_any_concentration_value = False

    # Skip folders that do not follow the expected structure (e.g., ignore irrelevant folders)
    full_protein_path = os.path.join(base_path, protein_folder)
    if not os.path.isdir(full_protein_path):
        continue

    # Extract details from the folder name
    folder_parts = protein_folder.split('_')
    if len(folder_parts) < 3:
        print(f"Skipping folder {protein_folder}: unexpected naming format.")
        continue

    # Extract the bio sample
    bio_sample = folder_parts[0]  # e.g., A549
    # Extract the protein name and remove the '-human' suffix
    protein_name = folder_parts[1].replace('-human', '')  # e.g., PBX3
    # Extract the experiment ID and remove the '.fa' suffix
    experiment_id = folder_parts[2].replace('.fa', '')  # e.g., ENCFF732ZLF
    # Construct the path to the nested directory containing the HTML file
    nested_folder_path = os.path.join(full_protein_path, os.path.basename(full_protein_path))

    html_file_path = f'file:///{nested_folder_path}/meme-chip.html'

    # Load the HTML file
    print(f"Loading HTML file for {protein_name}...")
    driver.get(html_file_path)

    # Allow time for the JavaScript to execute
    time.sleep(3)  # Adjust the sleep time if necessary

    # Click the "Expand All Clusters" button
    try:
        expand_all_button = driver.find_element(By.XPATH,
                                                "//span[@class='action' and text()='Expand All Clusters']")
        expand_all_button.click()
        time.sleep(2)  # Wait for the clusters to expand
        print(f"Clicked 'Expand All Clusters' for {protein_name}.")
    except Exception as e:
        print(f"Error clicking 'Expand All Clusters' for {protein_name}: {e}")

    # Check if the HTML file is loaded
    if "404 Not Found" in driver.title:
        print(f"Error: {html_file_path} not found")
        continue

    # Extract the groups and motifs
    print("Extracting groups and motifs...")
    motif_boxes = driver.find_elements(By.CSS_SELECTOR, ".motifbox")
    print(f"Found {len(motif_boxes)} motif boxes in {protein_name}.")

    for box_index, box in enumerate(motif_boxes):
        try:
            # Extract group information
            group_table = box.find_element(By.CSS_SELECTOR, "table.motifs")

            try:
                group_name_element = box.find_element(By.CSS_SELECTOR, "span.action a")
                group_name = group_name_element.text.strip()
            except Exception as e:
                group_name = "Unknown Group"

            group_id = f"Group_{group_counter}"
            group_counter += 1
            rows = group_table.find_elements(By.CSS_SELECTOR, "tbody tr")

            # Extract motif information for each row
            for row_index, row in enumerate(rows):
                cols = row.find_elements(By.TAG_NAME, "td")
                if len(cols) > 5:
                    # Extract the link from the "Discovery/Enrichment Program" column
                    try:
                        link_element = cols[1].find_element(By.TAG_NAME, "a")
                        link_url = link_element.get_attribute("href") if link_element else ""
                    except Exception as e:
                        link_url = ""

                    # Extract all links from the "SpaMo & FIMO" column
                    spamo_links = cols[5].find_elements(By.TAG_NAME, "a")
                    spamo_link_urls = [link.get_attribute("href") for link in spamo_links]
                    spamo_links_combined = "; ".join(spamo_link_urls)

                    # Find all occurrences of special characters
                    special_char_pattern = r'(\+|=|#motif_)'
                    all_matches = list(re.finditer(special_char_pattern, link_url))

                    # Get the last match
                    if all_matches:
                        last_match = all_matches[-1]
                        start_index = last_match.end()
                        segment = link_url[start_index:]
                    else:
                        segment = None

                    # Now extract the single motif ID using the final pattern
                    if segment:
                        final_pattern = r'^[^\s]+'
                        single_motif_id_match_final = re.search(final_pattern, segment)
                        single_motif_id = single_motif_id_match_final.group(
                            0) if single_motif_id_match_final else None
                    else:
                        single_motif_id = None

                    motif_data = {
                        "Bio Sample": bio_sample,
                        "Protein": protein_name,
                        "Experiment ID": experiment_id,
                        "FIMO Group Name": "",  # Initialize as empty
                        "Group ID": group_id,
                        "FIMO Motif ID": "",  # Will be filled later
                        "single Motif ID": single_motif_id,
                        "Enrichment Program": cols[1].text.strip() if len(cols) > 1 else "",
                        "single Motif source": link_url,
                        "meme-chip E-value": cols[2].text.strip() if len(cols) > 2 else "",
                        "Known or Similar Motifs": cols[3].text.strip().replace("\n", "; ") if len(
                            cols) > 3 else "",
                        "Distribution": "",  # Will be filled later
                        "SpaMo": spamo_links_combined,
                        "FIMO gff": "",  # To be filled later
                        "FIMO Motif source file": ""  # To be filled later
                    }
                    all_groups.append(motif_data)
        except Exception as e:
            continue

    # Extract FIMO commands from the HTML
    print("Extracting FIMO commands...")
    fimo_commands = driver.find_elements(By.XPATH, "//script[contains(., 'fimo')]")
    for command in fimo_commands:
        command_text = command.get_attribute('innerHTML')
        fimo_matches = re.findall(
            r'fimo --parse-genomic-coord --verbosity 1 --oc (.*?) --bgfile .*? --motif (.*?) (.*?) (.*?)',
            command_text)
        for match in fimo_matches:
            fimo_path, motif_id, motif_file, fasta_file = match
            fimo_folder = fimo_path.split('/')[-1]  # Extract the folder name
            for group in all_groups:
                if motif_id in group['single Motif source']:
                    group['FIMO gff'] = f"{fimo_folder}/fimo.gff"
                    group['FIMO Motif ID'] = motif_id
                    group['FIMO Motif source file'] = motif_file
                    group['FIMO Group Name'] = fimo_folder  # Update Group Name with FIMO folder name

    print(f"Extracted motifs for {protein_name}.")

    # Process the centrimo.html to extract concentration values
    centrimo_html_file_path = f'file:///{nested_folder_path}/centrimo_out/centrimo.html'
    driver.get(centrimo_html_file_path)
    time.sleep(3)  # Adjust the sleep time if necessary

    # Extract the concentration values
    centrimo_rows = driver.find_elements(By.CSS_SELECTOR, "#motifs tbody tr")
    centrimo_data = []
    for row in centrimo_rows:
        cols = row.find_elements(By.TAG_NAME, "td")
        if len(cols) > 5:
            data_row = {
                "Motif ID": cols[2].text.strip(),
                "Concentration": cols[5].text.strip()
            }
            centrimo_data.append(data_row)

    # Update the Distribution column with concentration values
    for group in all_groups:
        if group['Protein'] == protein_name:
            for centrimo in centrimo_data:
                if centrimo["Motif ID"] == group["single Motif ID"]:
                    group["Distribution"] = centrimo["Concentration"]
                    found_any_concentration_value = True
    print(f"found_any_concentration_value in {protein_name}: {found_any_concentration_value}")

# Close the WebDriver
driver.quit()

# Post-process the data to copy only specified columns by group name
all_groups_df = pd.DataFrame(all_groups)
for group_id in all_groups_df['Group ID'].unique():
    group_data = all_groups_df[all_groups_df['Group ID'] == group_id]
    if not group_data.empty:
        all_groups_df.loc[
            all_groups_df['Group ID'] == group_id, ['FIMO Group Name', 'SpaMo', 'FIMO gff', 'FIMO Motif source file',
                                                    'FIMO Motif ID']] = group_data.iloc[0][
            ['FIMO Group Name', 'SpaMo', 'FIMO gff', 'FIMO Motif source file', 'FIMO Motif ID']].values

# Ensure that "FIMO Group Name" is empty if "FIMO gff" is empty
all_groups_df.loc[all_groups_df['FIMO gff'] == "", 'FIMO Group Name'] = ""

# Create a DataFrame for all motifs and save to a CSV file
output_csv_path = "C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/meme/meme-chip_results/all_single_experiment_motifs_groups_all.csv"
all_groups_df.to_csv(output_csv_path, index=False)

print(f"All motifs saved to: {output_csv_path}")
