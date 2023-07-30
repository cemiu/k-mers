#!/bin/bash

mdir=$(dirname $(realpath "$0"))
cd "$mdir"

for file in "bin/fasta_to_sqlite" "bin/post_process_kmers" "bin/extract_pdb_coordinates"; do
    if ! [ -e "$file" ]; then
        echo "Could not locate binaries. Starting compilation."
        scripts/buildcpp.sh
        break
    fi
done

# check for pdb files
if ! find pdb/ -name '*.ent.gz' -print -quit | grep -q .; then
    echo "Error: No pdb files found in the pdb directory or its subdirectories"
    exit 1
fi

# # List of uniprot files.
# # SwissProt only: 92 MB .gz file to 250 MB database
# files=("uniprotkb/uniprot_sprot.fasta.gz")

# # SwissProt+TrEMBL: 62 GB .gz files to ~191 GB database
# # files=("uniprotkb/uniprot_sprot.fasta.gz", "uniprotkb/uniprot_trembl.fasta.gz")

options=("Process all PDBs" "Process PDBs matching UniProt" "Exit")
echo "Select processing type:"
select process_option in "${options[@]}"; do [[ "$process_option" ]] && break; done
[ "$process_option" == "Exit" ] && { echo "Exiting."; exit 0; }

uniprot_files=("uniprotkb/uniprot_sprot.fasta.gz" "uniprotkb/uniprot_trembl.fasta.gz")

if [[ "$process_option" == "Process PDBs matching UniProt" ]]; then
    if [[ ! -f "$dbfile" ]]; then
        echo "Database file not found."
        # checking for uniprot files
        found_files=0
        for filepath in "${uniprot_files[@]}"; do
            if [[ -f "$mdir/$filepath" ]]; then
                found_files=$((found_files + 1))
            fi
        done
        if (( found_files == 0 )); then
            echo "Uniprot files are missing. These are necessary to match PDBs with UniProt."
            echo "You can download the required files from https://www.uniprot.org/downloads and place them in the 'uniprotkb' directory."
            echo "The required files are uniprot_sprot.fasta.gz (SwissProt only), which expands from a 92 MB .gz file to a 250 MB database,"
            echo "or both uniprot_sprot.fasta.gz and uniprot_trembl.fasta.gz (SwissProt+TrEMBL), which expand from 62 GB .gz files to a ~191 GB database."
            exit 1
        else
            echo "Creating database."
            # Check if gzip is installed
            command -v gzip >/dev/null 2>&1 || { echo >&2 "gzip required but it's not installed. Aborting."; exit 1; }
            # Only process the files found
            for filepath in "${uniprot_files[@]}"; do
                if [[ -f "$mdir/$filepath" ]]; then
                    # logic for creating database here...
                    echo "Processing $filepath..."
                fi
            done

            echo "Creating index."
            sqlite3 "$dbfile" "CREATE INDEX IF NOT EXISTS idx_id ON sequences (id);"
            echo "Done creating database."
        fi
    fi
fi

k=12
while : ; do
    read -p "Enter k-mer length (default is 12): " input_k
    if [[ -z "$input_k" ]] || ( [[ ! "$input_k" =~ ^[0-9]+$ ]] || [[ "$input_k" -le 100 ]] ); then
        echo "Invalid input. Please enter a valid integer greater than 100 or press enter to use the default value."
    else
        k="$input_k"
        break
    fi
done

echo "Generating k-mers from PDBs."

# run python pipeline
if [[ "$process_option" == "Process all PDBs" ]]; then
    PYTHONPATH="${PYTHONPATH}:$mdir" python3 kmers/pipeline.py --handle_all_pdbs true
    echo "Extracting most frequent k-mers of length k=$k"
    ./bin/post_process_kmers -a -k "$k" > kmers.txt
else
    PYTHONPATH="${PYTHONPATH}:$mdir" python3 kmers/pipeline.py --handle_all_pdbs false
    echo "Extracting most frequent k-mers of length k=$k"
    ./bin/post_process_kmers -k "$k" > kmers.txt
fi
echo "Done generating k-mers."
echo "Finished. Results in \`kmers.txt\`"
echo "Top k-mers"
cat kmers.txt | head
