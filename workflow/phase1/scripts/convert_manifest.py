import pandas as pd
import sys

if len(sys.argv) < 2:
    print("Usage: python convert_manifest.py <path_to_csv>")
    sys.exit(1)

input_csv = sys.argv[1]
df = pd.read_csv(input_csv)

# Combine run_id and bioproject as required by your current Smk function
# Then append the scientificname/hostname for your reference
with open('download_manifest.txt', 'w') as f:
    for index, row in df.iterrows():
        # Keep the format 'run_id:bioproject' at the start
        # Add the hostname/scientificname at the end
        f.write(f"{row['run_id']}:{row['bioproject']},{row['scientificname']}\n")

print(f"Successfully created download_manifest.txt with host info from {input_csv}")