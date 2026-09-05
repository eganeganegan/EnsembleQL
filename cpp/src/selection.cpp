#include "ensembleql/selection.hpp"

#include <algorithm>
#include <cctype>
#include <functional>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace ensembleql {
namespace {
using AtomPredicate = std::function<bool(const Atom&)>;

std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    return value.substr(first, value.find_last_not_of(" \t\r\n") - first + 1);
}

std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

std::string upper(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::toupper(c));
    });
    return value;
}

std::vector<std::string> tokenize(const std::string& expression) {
    std::vector<std::string> tokens;
    std::string current;
    auto flush = [&]() {
        if (!current.empty()) {
            tokens.push_back(current);
            current.clear();
        }
    };
    for (const char character : expression) {
        if (std::isspace(static_cast<unsigned char>(character))) {
            flush();
        } else if (character == '(' || character == ')') {
            flush();
            tokens.emplace_back(1, character);
        } else {
            current.push_back(character);
        }
    }
    flush();
    return tokens;
}

class SelectionParser {
public:
    explicit SelectionParser(const std::string& expression)
        : expression_(expression), tokens_(tokenize(expression)) {}

    AtomPredicate parse() {
        if (tokens_.empty()) fail("selection cannot be empty");
        auto predicate = parse_or();
        if (!at_end()) fail("unexpected token '" + peek() + "'");
        return predicate;
    }

private:
    AtomPredicate parse_or() {
        auto left = parse_and();
        while (match("or")) {
            auto right = parse_and();
            left = [left = std::move(left), right = std::move(right)](const Atom& atom) {
                return left(atom) || right(atom);
            };
        }
        return left;
    }

    AtomPredicate parse_and() {
        auto left = parse_not();
        while (match("and")) {
            auto right = parse_not();
            left = [left = std::move(left), right = std::move(right)](const Atom& atom) {
                return left(atom) && right(atom);
            };
        }
        return left;
    }

    AtomPredicate parse_not() {
        if (match("not")) {
            auto operand = parse_not();
            return [operand = std::move(operand)](const Atom& atom) { return !operand(atom); };
        }
        return parse_primary();
    }

    AtomPredicate parse_primary() {
        if (match("(")) {
            auto predicate = parse_or();
            if (!match(")")) fail("expected ')'");
            return predicate;
        }
        if (at_end()) fail("expected selection predicate");
        const std::string keyword = lower(advance());
        if (keyword == "protein") {
            return [](const Atom& atom) { return atom.residue_name != "HOH" && atom.residue_name != "WAT"; };
        }
        if (keyword == "hydrogen") {
            return [](const Atom& atom) {
                if (!atom.element.empty()) return upper(atom.element) == "H";
                const auto first = std::find_if(atom.name.begin(), atom.name.end(), [](unsigned char c) {
                    return std::isalpha(c);
                });
                return first != atom.name.end() && std::toupper(static_cast<unsigned char>(*first)) == 'H';
            };
        }
        if (keyword != "resid" && keyword != "name" && keyword != "resname" && keyword != "chain") {
            fail("unknown selection keyword '" + keyword + "'");
        }
        if (at_end() || peek() == "(" || peek() == ")" || lower(peek()) == "and" || lower(peek()) == "or") {
            fail("expected a value after '" + keyword + "'");
        }
        const std::string value = advance();
        if (keyword == "resid") {
            int first{}, last{};
            const auto colon = value.find(':');
            try {
                std::size_t consumed{};
                const std::string start = value.substr(0, colon);
                first = std::stoi(start, &consumed);
                if (consumed != start.size()) throw std::invalid_argument("trailing text");
                if (colon == std::string::npos) {
                    last = first;
                } else {
                    const std::string end = value.substr(colon + 1);
                    consumed = 0;
                    last = std::stoi(end, &consumed);
                    if (end.empty() || consumed != end.size()) throw std::invalid_argument("trailing text");
                }
            } catch (const std::exception&) {
                throw std::runtime_error("Unknown residue selection: resid " + value);
            }
            if (first > last) throw std::runtime_error("Residue range start exceeds end: " + value);
            return [first, last](const Atom& atom) { return atom.residue_index >= first && atom.residue_index <= last; };
        }
        if (keyword == "name") return [value](const Atom& atom) { return atom.name == value; };
        if (keyword == "resname") return [value](const Atom& atom) { return atom.residue_name == value; };
        if (value.size() != 1) throw std::runtime_error("Chain selection expects one character: " + value);
        return [chain = value.front()](const Atom& atom) { return atom.chain == chain; };
    }

    bool match(std::string_view expected) {
        if (at_end() || lower(peek()) != expected) return false;
        ++position_;
        return true;
    }
    bool at_end() const noexcept { return position_ >= tokens_.size(); }
    const std::string& peek() const { return tokens_[position_]; }
    std::string advance() { return tokens_[position_++]; }
    [[noreturn]] void fail(const std::string& reason) const {
        throw std::runtime_error("Malformed selection '" + expression_ + "': " + reason);
    }

    std::string expression_;
    std::vector<std::string> tokens_;
    std::size_t position_{};
};
} // namespace

Selection::Selection(std::string expression) : expression_(trim(std::move(expression))) {
    if (expression_.empty()) throw std::invalid_argument("Selection cannot be empty");
}

std::vector<std::size_t> Selection::resolve(const Topology& topology) const {
    const auto predicate = SelectionParser(expression_).parse();
    std::vector<std::size_t> result;
    for (const Atom& atom : topology.atoms()) {
        if (predicate(atom)) result.push_back(atom.index);
    }
    if (result.empty()) throw std::runtime_error("Selection '" + expression_ + "' matched zero atoms");
    return result;
}

} // namespace ensembleql
