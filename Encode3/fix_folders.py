import os
import shutil

base_path = r"D:\Users\Daniel Batyrev\Documents\GitHub"
keyword = "_CTCF-human_"

# Find all directories that contain the specified keyword
dirs_to_fix = [d for d in os.listdir(base_path) if keyword in d]

for dir_name in dirs_to_fix:
    inner_path = os.path.join(base_path, dir_name, dir_name)
    if os.path.exists(inner_path):
        # Move all files from the inner directory to the outer directory
        for filename in os.listdir(inner_path):
            shutil.move(os.path.join(inner_path, filename), os.path.join(base_path, dir_name))
        
        # Remove the now empty inner directory
        os.rmdir(inner_path)

print("Folder structure fixed.")
