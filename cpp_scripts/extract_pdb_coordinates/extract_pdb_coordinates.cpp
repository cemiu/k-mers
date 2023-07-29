#include <iostream>
#include <unordered_map>
#include <string>
#include <vector>
#include <unordered_set>
#include <sstream>

#include "AtomDataParser.h"
#include "Constants.h"
#include "Utils.h"


// bool invalidSequence = false;
bool invalidAA = false;
bool hasUnknownResidues = false;
std::stringstream parsedSequence;
std::vector<std::string> errorOutput;

int firstCAResidue = 0;

ResidueConfirmation validateAtomSequence(int &prevCAResiduePosition, const int &resSeq) {
    if (prevCAResiduePosition + 1 != resSeq) {
        if (prevCAResiduePosition == -1) { // initial value
            prevCAResiduePosition = resSeq;
            firstCAResidue = resSeq;
            return RESIDUE_VALID;
        }

        if (prevCAResiduePosition == resSeq)
            return RESIDUE_DUPLICATE;

        std::stringstream errorString;
        errorString << "missing residues; prev=" << prevCAResiduePosition << ", next=" << resSeq;
        errorOutput.push_back(errorString.str());

        prevCAResiduePosition = resSeq;
        return RESIDUE_OUT_OF_SEQUENCE;
    }
    prevCAResiduePosition = resSeq;
    return RESIDUE_VALID;
}

//////////////////////////
//////////////////////////
///// PROCESS ROWS ///////
//////////////////////////
//////////////////////////

void processAtom(std::string &line, std::vector<std::string> &output, int &prevCAResiduePosition, bool &isSequenceValid) {
    AtomData data;
    parseAtomData(line, data, 0);
    if (!data.isValidAtom) {
        return;
    }

    switch(validateAtomSequence(prevCAResiduePosition, data.resSeq)) {
    case RESIDUE_VALID:
        break; // continue
    case RESIDUE_DUPLICATE:
        return; // skip to next
    case RESIDUE_OUT_OF_SEQUENCE:
        isSequenceValid = false;
        break;
    }

    char aminoAcid;
    try {
        aminoAcid = aminoAcidLookup.at(data.resName);
    } catch (std::out_of_range) { // should never throw if pdb is valid & not unknown
        if (data.resName == "UNK") {
            hasUnknownResidues = true;
            return;
        }

        throw std::runtime_error("Unexpected atom type: " + data.resName);
    }

    // Selenocysteine, Pyrrolysine, GLX, ASX, too rare, skip
    if (aminoAcid == 'U' || aminoAcid == 'O' || aminoAcid == 'Z' || aminoAcid == 'B') {
        invalidAA = true;
    }

    // construct output string
    std::stringstream ss;
    ss << aminoAcid << ' ' << data.x << ' '  << data.y << ' ' << data.z;
    // std::cout << aminoAcid << ' ' << data.x << ' '  << data.y << ' ' << data.z << std::endl;

    // construct sequence string
    parsedSequence << aminoAcid;

    output.push_back(ss.str());
}

// Checks whether the parsed input, so far, produced a valid, sequential
// list of residues with coordinates.
// Returns PDBParsingCode.SUCCESS if successful, and a specific error code
// otherwise.
PDBParsingCode isPDBInvalid(float &resolution, bool &isSequenceValid, int &prevCAResiduePosition, bool &anyCAAtomsPresent) {
    bool isResolutionValid = resolution < 2.5;
    if (!isResolutionValid) { // resolution too low
        return RESOLUTION_TOO_LOW;
    }

    if (resolution == -1) // no valid resolution remark returned
        return RESOLUTION_NOT_SPECIFIED;

    if (!isSequenceValid) // missing non-terminal residues
        return MISSING_NON_TERMINAL_RESIDUES;

    if (prevCAResiduePosition == -1) { // no single CA atom found
        if (anyCAAtomsPresent) // if any model had, but last one didn't
            return MISSING_NON_TERMINAL_RESIDUES;
        return NO_ALPHA_CARBON_ATOMS_FOUND;
    }

    return SUCCESS;
}

PDBParsingCode isPDBInvalid(float &resolution, bool &isSequenceValid, int &prevCAResiduePosition, bool &anyCAAtomsPresent, bool &isNotProtein) {
    if (isNotProtein)
        return IS_NOT_PROTEIN;
    // if (invalidSequence)
    //     return INVALID_SEQUENCE;
    if (hasUnknownResidues)
        return HAS_UNKNOWN_RESIDUE;
    if (invalidAA)
        return EXCLUDE_RARE_AMINO_ACIDS;
    return isPDBInvalid(resolution, isSequenceValid, prevCAResiduePosition, anyCAAtomsPresent);
}

void resetPDBOutput(std::vector<std::string> &output, bool &isSequenceValid, int &prevCAResiduePosition, bool &anyCAAtomsPresent) {
    if (!anyCAAtomsPresent && output.size()) {
        anyCAAtomsPresent = true;
    }

    output.clear();
    isSequenceValid = true;
    hasUnknownResidues = false;
    invalidAA = false;
    prevCAResiduePosition = -1;
    firstCAResidue = 0;
    parsedSequence.str("");
}

bool hasMatched = true;

std::vector<std::string> processSequences(std::unordered_map<char, std::stringstream> &sequenceStreams) {
    std::unordered_set<std::string> uniqueSequences;
    std::vector<std::string> sequences;
    std::string matchedSequence = "N/A";

    if (sequenceStreams.size() == 0)
        return sequences;

    for (const auto & [_chainId, stream] : sequenceStreams) {
        auto sequence = stream.str();
        if (sequence.size() != 0) {
            auto sequencePosition = sequence.find(parsedSequence.str());
            if (parsedSequence.str().size() > 0 && sequencePosition != std::string::npos) {
                matchedSequence = sequence;
            } else {
                uniqueSequences.insert(sequence);
            }
        }
    }

    // line 4: matched sequence (parsed contained within matched)
    if (matchedSequence != "")
        std::cout << "matched: " << matchedSequence << std::endl;
    // line 5: sequence parsed from ATOM records
    std::cout << "parsed:  " << parsedSequence.str() << std::endl;

    sequences.push_back(matchedSequence);

    // line 6+: all other parsed sequences
    for (std::string seq: uniqueSequences) {
        std::cout << "other:   " << seq << std::endl;
        sequences.push_back(seq);
    }

    if (matchedSequence == "N/A")
        hasMatched = false;

    return sequences;
}

// Small script for parsing PDB files.
// Takes in a stream of a PDB file as input.
int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

    std::vector<std::string> output;
    std::unordered_set<std::string> uniprotIds;
    bool isSequenceValid = true;
    bool processed_atom = false;
    bool anyCAAtomsPresent = false;
    bool isNotProtein = false;
    auto resolution = -1.0f;
    int prevCAResiduePosition = -1;

    std::unordered_map<char, std::stringstream> sequenceStreams;
    
    std::string line;
    while (getline(std::cin, line)) {
        std::string param = line.substr(0, 6);

        if (param == "HEADER") {
            auto headerType = processHeader(line);
            if (headerType != PROTEIN) {
                isNotProtein = true;
                // break;
            }
        } else if (param == "REMARK") {
            processRemark(line, resolution);
        } else if (param == "DBREF ") {
            processDBRef(line, uniprotIds);
        } else if (param == "DBREF1") {
            processDBRef1(line, uniprotIds);
        } else if (param == "SEQRES") {
            processSequence(line, sequenceStreams);
        } else if (param == "ATOM  ") { // HETATM residues are skipped
            processAtom(line, output, prevCAResiduePosition, isSequenceValid);
        } else if (param == "TER   ") { // end of one chain
            auto pdbValidity = isPDBInvalid(resolution, isSequenceValid, prevCAResiduePosition, anyCAAtomsPresent);
            // std::cout << code_name[pdbValidity] << std::endl;
            if (pdbValidity == SUCCESS)
                break; // terminate parser, output PDB

            // if at first you don't succeed, try, try again (parse next model)
            resetPDBOutput(output, isSequenceValid, prevCAResiduePosition, anyCAAtomsPresent);
        } // else ignore line, until end is reached
    }

    // line 1 -- pdb id
    // "pdb_id:  201L" (printed in Utils.cpp)

    // line 2 -- resolution
    std::cout << "resolut: " << resolution << std::endl;

    // line 3 -- uniprot IDs
    std::string allUniprotIds = concatenateString(uniprotIds);
    std::cout << "uniprot: " << allUniprotIds << std::endl;

    // line 4 -- matched sequence (atom record substring of reqres)
    // line 5 -- parsed sequence (atom records)
    // line 6-n -- other sequences (reqres sequence)
    std::vector<std::string> sequences = processSequences(sequenceStreams);
    

    auto pdbValidity = isPDBInvalid(resolution, isSequenceValid, prevCAResiduePosition, anyCAAtomsPresent, isNotProtein);
    if (pdbValidity != SUCCESS) {
        std::cerr << code_name[pdbValidity] << std::endl;

        for (auto error : errorOutput)
            std::cout << error << std::endl;
        return pdbValidity;
    }

    if (uniprotIds.size() == 0) {
        std::cerr << code_name[NO_UNIPROT_ID] << std::endl;
        return NO_UNIPROT_ID;
    }

    // if (!hasMatched) {
    //     std::cerr << "NO_MATCH" << std::endl;
    //     return 103;
    //     // std::cerr << "NO MATCHED SEQUENCE" << std::endl;
    // }

    // line n+1: sequence number of initial residue (starts with 1)
    std::cout << "initres: " << firstCAResidue << std::endl;

    // line n+2: empty line
    std::cout << std::endl;

    // lines n+3 to end: coordinates in format <residue> <x> <y> <z>
    for (std::string pos : output) {
        std::cout << pos << std::endl;
    }

    return 0;
}

