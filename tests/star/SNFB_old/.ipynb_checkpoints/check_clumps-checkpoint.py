import glob
import os

def concatenate_clump_files():
    # Find all files matching the pattern
    pattern = "output_00002/clump_00002.txt0000*"
    files = glob.glob(pattern)
    
    if not files:
        print(f"No files found matching pattern: {pattern}")
        return
    
    # Sort files to ensure consistent ordering
    files.sort()
    print(f"Found {len(files)} files: {files}")
    
    output_file = "clump_all.txt"
    header_written = False
    
    with open(output_file, 'w') as outfile:
        for i, filename in enumerate(files):
            print(f"Processing {filename}...")
            
            with open(filename, 'r') as infile:
                lines = infile.readlines()
                
                # Skip empty files
                if not lines:
                    continue
                
                # Write header only from the first file
                if not header_written and lines:
                    outfile.write(lines[0])  # Write header
                    header_written = True
                
                # Write data lines (skip header for all files)
                for line in lines[1:]:
                    outfile.write(line)
    
    print(f"Successfully concatenated {len(files)} files into {output_file}")
    
    # Show summary
    with open(output_file, 'r') as f:
        total_lines = sum(1 for _ in f)
    print(f"Clump concatenated output file contains {total_lines} lines (including 1 header)")

if __name__ == "__main__":
    concatenate_clump_files()


