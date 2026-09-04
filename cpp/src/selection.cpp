#include "ensembleql/selection.hpp"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>

namespace ensembleql {
namespace {
std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    return value.substr(first, value.find_last_not_of(" \t\r\n") - first + 1);
}
std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}
} // namespace

Selection::Selection(std::string expression) : expression_(trim(std::move(expression))) {
    if (expression_.empty()) throw std::invalid_argument("Selection cannot be empty");
}

std::vector<std::size_t> Selection::resolve(const Topology& topology) const {
    std::istringstream input(expression_);
    std::string keyword, value, extra;
    input >> keyword;
    keyword = lower(keyword);
    std::vector<std::size_t> result;
    if (keyword == "protein") {
        if (input >> extra) throw std::runtime_error("Malformed selection '" + expression_ + "'");
        for (const Atom& atom : topology.atoms()) if (atom.residue_name != "HOH" && atom.residue_name != "WAT") result.push_back(atom.index);
    } else {
        if (!(input >> value) || (input >> extra)) throw std::runtime_error("Malformed selection '" + expression_ + "'");
        if (keyword == "resid") {
            int first{}, last{};
            const auto colon = value.find(':');
            try {
                first = std::stoi(value.substr(0, colon));
                last = colon == std::string::npos ? first : std::stoi(value.substr(colon + 1));
            } catch (const std::exception&) { throw std::runtime_error("Unknown residue selection: resid " + value); }
            if (first > last) throw std::runtime_error("Residue range start exceeds end: " + value);
            for (const Atom& atom : topology.atoms()) if (atom.residue_index >= first && atom.residue_index <= last) result.push_back(atom.index);
        } else if (keyword == "name") {
            for (const Atom& atom : topology.atoms()) if (atom.name == value) result.push_back(atom.index);
        } else if (keyword == "resname") {
            for (const Atom& atom : topology.atoms()) if (atom.residue_name == value) result.push_back(atom.index);
        } else if (keyword == "chain") {
            if (value.size() != 1) throw std::runtime_error("Chain selection expects one character: " + value);
            for (const Atom& atom : topology.atoms()) if (atom.chain == value.front()) result.push_back(atom.index);
        } else {
            throw std::runtime_error("Unknown selection keyword: " + keyword);
        }
    }
    if (result.empty()) throw std::runtime_error("Selection '" + expression_ + "' matched zero atoms");
    return result;
}

} // namespace ensembleql
