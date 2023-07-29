#ifndef ATOMDATAPARSER_H
#define ATOMDATAPARSER_H

#include <string>

struct AtomData
{
    bool isValidAtom; // whether the atom is valid (CA) 
    std::string resName; // residue name (AA)
    int resSeq; // residue sequence number
    float x, y, z;
};

void parseAtomData(const std::string& str, AtomData& data, size_t offset = 0);

#endif // ATOMDATAPARSER_H
