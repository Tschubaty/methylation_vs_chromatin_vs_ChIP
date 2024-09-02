import os

# Base directory path
base_dir = r'D:\Users\Daniel Batyrev\Documents\GitHub\meme\single_zip'

# Walk through all subdirectories in the base directory
for root, dirs, files in os.walk(base_dir):
    # Process only directories that start with "fimo_out_"
    if os.path.basename(root).startswith('fimo_out_'):
        super_folder = os.path.basename(os.path.dirname(root))
        print(f"Processing super folder: {super_folder}, subfolder: {os.path.basename(root)}")
        
        input_file_path = os.path.join(root, 'fimo.tsv')
        output_file_path = os.path.join(root, 'fimo_cpg_sites.bed')

        # Check if the input file exists
        if not os.path.isfile(input_file_path):
            print(f"  - Skipping {root} as it does not contain fimo.tsv")
            continue  # Skip if fimo.tsv is not found

        # Open the input and output files
        with open(input_file_path, 'r') as infile, open(output_file_path, 'w') as outfile:
            for line in infile:
                # Skip header lines and comments
                if line.startswith('#') or line.startswith('motif_id'):
                    continue

                # Split the line into columns
                columns = line.strip().split('\t')

                # Ensure the line has the correct number of columns
                if len(columns) != 10:
                    continue  # Skip lines that do not have exactly 10 columns

                # Unpack the columns
                motif_id, motif_alt_id, sequence_name, start, stop, strand, score, p_value, q_value, matched_sequence = columns

                # Convert start to integer
                start = int(start)

                # Find CpG sites in the matched sequence
                for i in range(len(matched_sequence) - 1):
                    if matched_sequence[i:i+2].upper() == 'CG':
                        if strand == '+':
                            cpg_start = start + i - 1  # BED format is 0-based, subtract 1 from start
                            cpg_end = cpg_start + 2  # End is 1-based, so add 2 for CpG
                        elif strand == '-':
                            cpg_end = start + (len(matched_sequence) - 1 - i) + 1  # End position
                            cpg_start = cpg_end - 2  # Start position is 2 bases before the end

                        # Write to BED file
                        outfile.write(f'{sequence_name}\t{cpg_start}\t{cpg_end}\t{motif_id}\t{score}\t{strand}\n')

        print(f"  - Finished processing {root}")

print("Processing complete.")
