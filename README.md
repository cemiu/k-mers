# k-mers
Description and purpose of the project will come at a later point.

## Installation

### Prerequisites
`g++`: through gcc installation

`gzip / gunzip`: usually pre-installed

`python3`: installed and part of path

optional: `sqlite3`: Debian/Ubuntu: `apt install sqlite3`; macOS: `brew install sqlite3`

### Preparation

Clone the repo
```
git clone https://github.com/cemiu/kmers.git && cd kmers
```
Two folders need to be populated:
- `pdb` (required)
- `uniprotkb` (optional)
### pdb
The `pdb` folder has to contain experimental PDB files in the .ent.gz format.

Instructions for downloading can be found here:

https://www.rcsb.org/docs/programmatic-access/file-download-services

https://files.wwpdb.org/pub/pdb/data/structures/divided/pdb/

Alternatively a mirror can be found here: https://pycom.brunel.ac.uk/misc/pdb_2023-07-28.tar (42 GB)

Once downloaded they have to be placed in the `pdb` folder **without** being decompressed.
It does not matter whether they are divided; e.g. `pdb/file.ent.gz` or `pdb/<folder>/file.ent.gz`.

### uniprotkb
**Optionally**, the `uniprotkb` folder can be populated with the uniprotkb fasta files.
This is only needed, if the PDBs should be associated to a Protein. If this is not required, **skip this step**.

Files:
- `uniprot_sprot.fasta.gz` (400 MB after processing)
  - Use only Swiss-Prot, has the majority of PDB coverage
- Optionally, `uniprot_trembl.fasta.gz` (250 GB after processing)
  - Use TrEMBL to match more PDBs; might be useful for max. coverage

The latter might result in (slightly) more PDBs which can be associated to a Protein. The difference is expected to be trivial.

Place the files in the `uniprotkb` folder without decompressing them.

Once `run.sh` has been executed and the database `uniprotkb/uniprot_sequences.db` has been created, the files can be deleted.

### Running

To run the script, execute:
```
./run.sh
```

This will
- Compile C++ binaries, if they don't exist
- Ask whether to process the uniprotkb files
  - If yes, create the database `uniprotkb/uniprot_sequences.db`, if it doesn't exist
- Ask for k-mer size (default: k=12)
- Extracts 3d k-mer from the PDBs (into `pdb_output` folder)
- Extracts k-mer of length k into `kmer.txt`, along with frequency
