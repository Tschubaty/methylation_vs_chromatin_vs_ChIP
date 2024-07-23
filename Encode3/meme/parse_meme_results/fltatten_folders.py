import os
import shutil

# Base path to your extracted folders
base_path = 'D:/Users/Daniel Batyrev/Documents/GitHub/successZIP/extracted'

# Loop over all folders in the base path
for folder_name in os.listdir(base_path):
    if folder_name.startswith('pooled_'):
        inner_folder_path = os.path.join(base_path, folder_name, folder_name)
        if os.path.exists(inner_folder_path):
            print(f"Processing {inner_folder_path}...")
            # Move all contents of the inner folder to the outer folder
            for item in os.listdir(inner_folder_path):
                src = os.path.join(inner_folder_path, item)
                dest = os.path.join(base_path, folder_name)
                try:
                    if os.path.isdir(src):
                        shutil.move(src, dest)
                    else:
                        shutil.move(src, dest)
                except Exception as e:
                    print(f"Error moving {src} to {dest}: {e}")
            # Remove the now-empty inner folder
            try:
                os.rmdir(inner_folder_path)
            except Exception as e:
                print(f"Error removing {inner_folder_path}: {e}")

print("Flattening complete.")
