#include "AtomDataParser.h"
#include <string>
#include <sstream>
// #include <iostream>

void parseAtomData(const std::string& str, AtomData& data, size_t offset)
{
    std::string atom_name = str.substr(-offset + 12, 4);
    std::stringstream atom_ss(atom_name);
    atom_ss >> atom_name;
    
    data.isValidAtom = atom_name == "CA"; // whether the atom is ca

    data.resName = str.substr(-offset + 17, 3);
    data.resSeq = std::stoi(str.substr(-offset + 22, 4));
    data.x = std::stof(str.substr(-offset + 30, 8));
    data.y = std::stof(str.substr(-offset + 38, 8));
    data.z = std::stof(str.substr(-offset + 46, 8));
}
