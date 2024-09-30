import zipfile
import os

# Define the directory containing the zip files
zip_dir = r"D:\Users\Daniel Batyrev\Documents\GitHub"

# List all files in the directory
files = os.listdir(zip_dir)

# Iterate over all files
for file in files:
    # Check if the file is a zip file
    if file.endswith(".zip"):
        # Define the full path to the zip file
        zip_path = os.path.join(zip_dir, file)
        
        # Define the extraction directory
        extract_dir = os.path.join(zip_dir, file.replace(".zip", ""))
        
        # Create the extraction directory if it doesn't exist
        os.makedirs(extract_dir, exist_ok=True)
        
        # Extract the zip file
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)
        
        print(f"Extracted {file} to {extract_dir}")

print("Extraction complete.")
