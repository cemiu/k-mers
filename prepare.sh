#!/bin/bash

mdir=$(dirname $(realpath "$0"))
cd "$mdir"

# List of uniprot files.
# SwissProt only: 92 MB .gz file to 250 MB database
files=("uniprotkb/uniprot_sprot.fasta.gz")

# SwissProt+TrEMBL: 62 GB .gz files to ~191 GB database
# files=("uniprotkb/uniprot_sprot.fasta.gz", "uniprotkb/uniprot_trembl.fasta.gz")


# Check if sqlite3 is installed
command -v sqlite3 >/dev/null 2>&1 || { echo >&2 "sqlite3 required but it's not installed (brew install sqlite3 [on macos]). Aborting."; exit 1; }

# Check if g++ is installed
command -v g++ >/dev/null 2>&1 || { echo >&2 "g++ required but it's not installed. Aborting."; exit 1; }

# Check if gunzip is installed
command -v gunzip >/dev/null 2>&1 || { echo >&2 "gunzip required but it's not installed. Aborting."; exit 1; }


# Compile C++ program
echo "Compiling C++ programs..."
arch -x86_64 g++ -std=c++17 -o "$mdir"/bin/fasta_to_sqlite "$mdir"/cpp_scripts/fasta_to_sqlite/*.cpp -lsqlite3
if [ $? -ne 0 ]; then
    echo "Compilation failed. Please check the source code."
    exit 1
fi

g++ -std=c++17 -o "$mdir"/bin/post_process_kmers "$mdir"/cpp_scripts/post_process_kmers/*.cpp
if [ $? -ne 0 ]; then
    echo "Compilation failed. Please check the source code."
    exit 1
fi

g++ -std=c++17 -o "$mdir"/bin/extract_pdb_coordinates "$mdir"/cpp_scripts/extract_pdb_coordinates/*.cpp
if [ $? -ne 0 ]; then
    echo "Compilation failed. Please check the source code."
    exit 1
fi

echo "Compiled all C++ programs."

# Check if the database file exists
dbfile="$mdir/uniprotkb/uniprot_sequences.db"
if [[ -f "$dbfile" ]]; then
    read -p "The file $dbfile already exists. Do you want to overwrite it? (yes/no): " response
    if [[ "$response" == "no" ]]; then
        echo "Aborting."
        exit 0
    elif [[ "$response" == "yes" ]]; then
        rm "$dbfile"
    else
        echo "Invalid response. Please type 'yes' or 'no'."
        exit 1
    fi
fi

for filepath in "${files[@]}"; do
    # Check if file exists
    if [[ ! -f "$mdir/$filepath" ]]; then
        echo "The required file $filepath is missing."
        echo "Please make sure that uniprot_trembl.fasta.gz and uniprot_sprot.fasta.gz are in the 'uniprotkb' directory."
        echo "The file can be downloaded from https://www.uniprot.org/downloads"
        exit 1
    fi
done

# Create database
for filepath in "${files[@]}"; do
    echo "Processing $filepath..."
    gunzip -c "$mdir/$filepath" | ./bin/fasta_to_sqlite
    echo "Done processing $filepath"
done

# Create indexes on the database
echo "Creating indexes..."
sqlite3 "$dbfile" "CREATE INDEX IF NOT EXISTS idx_id ON sequences (id);"

# this index is too large (doubles the DB size)
# sqlite3 $dbfile "CREATE INDEX IF NOT EXISTS idx_sequence ON sequences (sequence);"
echo "Indexes created."

echo "Extracting Residue Coordinates and Generating k-mers..."
PYTHONPATH="${PYTHONPATH}:$mdir" python3 kmers/pipeline.py --handle_all_pdbs true
echo "Done generating k-mers"

k=12

echo "Extracting most frequent k-mers of length k=$k"
./bin/post_process_kmers -a -k "$k > kmers.txt"
echo "Finished. Results in \`kmers.txt\`"
echo "Top k-mers"
cat kmers.txt | head
