#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>
#include <unordered_set>

#include "Constants.h"

std::string concatenateString(const std::vector<std::string>& strings);
std::string concatenateString(const std::unordered_set<std::string>& strings);

PDBType processHeader(const std::string &line);
float extractResolution(const std::string &line);
void processRemark(const std::string &line, float &resolution);
void processDBRef(const std::string &line, std::unordered_set<std::string> &uniprotIds);
void processDBRef1(const std::string &line, std::unordered_set<std::string> &uniprotIds);
void processSequence(const std::string &line, std::unordered_map<char, std::stringstream> &sequenceStreams);

#endif // UTILS_H
