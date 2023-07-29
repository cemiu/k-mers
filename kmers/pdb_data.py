
"""
PDBData takes a byte stream (produces by extract_pdb_coordinates) and parses it into a PDBData object.

The format of the byte stream is as follows:
    1. First line is resolution:
        'resolut: {resolution}' (float)
    2. Second line is uniprot_id. Can be comma separated if multiple uniprot_ids are found:
        'uniprot: {uniprot_id}' (str)
    3. Third line is matched sequence. This is the sequence that is found in the PDB file (SEQRES)
        'matched: {sequence}' (str)
    4. Fourth line is parsed sequence. These are the residues that are found in the PDB file (ATOM)
        'parsed:  {parsed_sequence}' (str)
    5. Lines 5 to n are other sequences that are found in the PDB file (SEQRES)
        'other:   {other_sequence}' (str)
    6. Lines n+1 is a blank line
    7. Lines n+2 to m are the coordinates of the PDB file (ATOM)
        '{residue_name [e.g. A/E/M]} {x} {y} {z}' (char, float, float, float)
"""


class PDBData:
    """Takes in a byte stream"""
    def __init__(self, pdb_byte_stream):
        self._pdb_id = None
        self._resolution = None
        self._parsed_sequence = None
        self._first_residue_number = None
        self._residue_list = []
        self._coordinates = []

        self._input_uniprot_ids = []
        self._input_matched_sequence = None
        self._input_other_sequences = []

        self._parse(pdb_byte_stream)

    def _parse(self, pdb_byte_stream):
        lines = pdb_byte_stream.decode('utf-8').split("\n")

        # Parse pdb id
        self._pdb_id = lines[0].split(':  ')[1]

        # Parse the resolution
        self._resolution = float(lines[1].split(': ')[1])

        # Parse the uniprot ids
        self._input_uniprot_ids = lines[2].split(': ')[1].split(',')

        # Parse the matched sequence
        self._input_matched_sequence = lines[3].split(': ')[1]
        # if self._input_matched_sequence == 'N/A':
        #     print('N/A')

        # Parse the parsed sequence
        self._parsed_sequence = lines[4].split(':  ')[1]

        # Parse the other sequences
        idx = 5
        while lines[idx].startswith('other:   '):
            self._input_other_sequences.append(lines[idx].split(': ')[1])
            idx += 1

        # parse the first residue number (initres)
        self._first_residue_number = int(lines[idx].split(': ')[1])

        # Parse the coordinates
        idx += 2  # skip the initres and blank line
        while idx < len(lines) and lines[idx].strip():  # Ensure we are not at the end or at a blank line
            parts = lines[idx].split()
            residue_name = parts[0]
            x, y, z = map(float, parts[1:])
            self._residue_list.append(residue_name)
            self._coordinates.append((x, y, z))
            idx += 1

    @property
    def pdb_id(self):
        return self._pdb_id

    @property
    def resolution(self):
        return self._resolution

    @property
    def first_residue_number(self):
        return self._first_residue_number

    @property
    def residue_sequence_parsed(self):
        """
        The parsed sequence is the sequence extracted from the PDB ATOM records.
        :return:
        """
        return self._parsed_sequence

    @property
    def residue_sequence_matched(self):
        """
        The matched sequence is the sequence extracted from the PDB SEQRES records,
        of which the parsed sequence is a substring of.
        :return:
        """
        return self._input_matched_sequence

    @property
    def uniprot_ids(self):
        return self._input_uniprot_ids

    @property
    def residue_list(self):
        return self._residue_list

    @property
    def coordinates(self):
        return self._coordinates
