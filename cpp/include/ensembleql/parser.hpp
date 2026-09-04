#pragma once

#include "ensembleql/ast.hpp"

#include <stdexcept>
#include <string>

namespace ensembleql {

class QueryError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class Parser {
public:
    ast::Query parse(const std::string& text) const;
};

} // namespace ensembleql
