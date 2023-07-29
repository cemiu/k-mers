# k-mers
Description and purpose of the project will come at a later point.

## Installation

### Prerequisites
`g++`: through gcc installation
`sqlite3`: Debian/Ubuntu: `apt install sqlite3`; macOS: `brew install sqlite3`
`gzip / gunzip`: usually pre-installed
`python3`: installed and part of path

### Preperation
Clone the repo
```
git clone https://github.com/cemiu/kmers.git && cd kmers
```
Two folders need to be populated.
#### pdb
The pdb folder has to contain all experimental PDB file in the .ent.gz format.
Instructions for downloading can be found here:
https://www.rcsb.org/docs/programmatic-access/file-download-services
https://files.wwpdb.org/pub/pdb/data/structures/divided/pdb/

Alternatively a mirror can be found here: https://pycom.brunel.ac.uk/misc/pdb_2023-07-28.tar (42 GB)

Once downloaded they have to be placed in the `pdb` folder **without** being uncompressed. It does not matter whether they are in `pdb/file.ent.gz` or `pdb/<folder>/file.ent.gz`.

#### uniprotkb
The project requires `uniprot_sprot.fasta.gz` (400 MB after processing)
Optionally, `uniprot_trembl.fasta.gz` can be used, to match more PDBs (250 GB after processing).

The latter might result in (slightly) more PDBs which can be associated to a Protein. The difference is expected to be trivial.

Place the files in the `uniprotkb` folder without uncompressing them.
By default, only Swiss-Prot is used. To also use TrEMBL, uncomment line 11 in `run.sh`.

### Running

To run the script, execute:
```
./run.sh
```

This will
- Compile C++ binaries
- Process Swiss-Prot / TrEMBL into a database 
  - (`uniprot_sprot.fasta.gz` / `uniprot_trembl.fasta.gz`) can be deleted afterwards
- Extract 3d k-mer from the PDBs
- TODO: process the k-mers
