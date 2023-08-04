#ifndef ATOMDATAPARSER_H
#define ATOMDATAPARSER_H

#include <string>
#include <sstream>

struct AtomData
{
    bool isValidAtom; // whether the atom is valid (CA)
    std::string resName; // residue name (AA)
    int resSeq; // residue sequence number
    float x, y, z;

    AtomData(const std::string& str, size_t offset = 0)
    {
        std::string atom_name = str.substr(-offset + 12, 4);
    
        isValidAtom = atom_name == " CA "; // whether the atom is ca
        if (!isValidAtom)
            return;

        resName = str.substr(-offset + 17, 3);
        resSeq = std::stoi(str.substr(-offset + 22, 4));
        x = std::stof(str.substr(-offset + 30, 8));
        y = std::stof(str.substr(-offset + 38, 8));
        z = std::stof(str.substr(-offset + 46, 8));
    }
};

#endif // ATOMDATAPARSER_H
