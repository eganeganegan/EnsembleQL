#include "ensembleql/parser.hpp"

#include <algorithm>
#include <cctype>
#include <optional>
#include <regex>
#include <sstream>
#include <vector>

namespace ensembleql {
namespace {
std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    return value.substr(first, value.find_last_not_of(" \t\r\n") - first + 1);
}
std::string upper(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return value;
}
std::vector<std::string> split_args(const std::string& text) {
    std::vector<std::string> result;
    std::size_t start = 0;
    int depth = 0;
    for (std::size_t i = 0; i < text.size(); ++i) {
        if (text[i] == '(') ++depth;
        else if (text[i] == ')') --depth;
        else if (text[i] == ',' && depth == 0) { result.push_back(trim(text.substr(start, i - start))); start = i + 1; }
    }
    result.push_back(trim(text.substr(start)));
    return result;
}

struct ContactOptions {
    Distance cutoff{0.45};
    ContactMode mode{ContactMode::Atom};
};

ContactOptions contact_options(const std::vector<std::string>& args, const std::string& observable,
                               bool allow_mode) {
    if (args.size() < 2) throw QueryError(observable + " requires two selections");
    ContactOptions options;
    bool cutoff_seen = false;
    bool mode_seen = false;
    for (std::size_t index = 2; index < args.size(); ++index) {
        std::string value = trim(args[index]);
        const auto equals = value.find('=');
        std::string key;
        if (equals == std::string::npos) {
            key = "CUTOFF"; // Preserve the positional third-argument syntax.
        } else {
            key = upper(trim(value.substr(0, equals)));
            value = trim(value.substr(equals + 1));
        }
        if (key == "CUTOFF") {
            if (cutoff_seen) throw QueryError("Duplicate " + observable + " cutoff option");
            options.cutoff = parse_distance(value);
            if (options.cutoff.nm < 0.0) throw QueryError("Contact cutoff must be non-negative");
            cutoff_seen = true;
        } else if (key == "MODE" && allow_mode) {
            if (mode_seen) throw QueryError("Duplicate " + observable + " mode option");
            const std::string mode = upper(value);
            if (mode == "ATOM") options.mode = ContactMode::Atom;
            else if (mode == "RESIDUE") options.mode = ContactMode::Residue;
            else throw QueryError("CONTACT_COUNT mode must be 'atom' or 'residue'");
            mode_seen = true;
        } else {
            throw QueryError("Unknown " + observable + " option: " +
                             (equals == std::string::npos ? value : trim(args[index].substr(0, equals))));
        }
    }
    return options;
}

bool rg_mass_weighted(const std::vector<std::string>& args) {
    if (args.empty() || args.size() > 2) {
        throw QueryError("RG requires one selection and an optional mass_weighted option");
    }
    if (args.size() == 1) return false;
    const auto equals = args[1].find('=');
    if (equals == std::string::npos || upper(trim(args[1].substr(0, equals))) != "MASS_WEIGHTED") {
        throw QueryError("Unknown RG option: " + trim(args[1].substr(0, equals)));
    }
    const std::string value = upper(trim(args[1].substr(equals + 1)));
    if (value == "TRUE") return true;
    if (value == "FALSE") return false;
    throw QueryError("RG mass_weighted option must be true or false");
}

class Cursor {
public:
    explicit Cursor(std::string text) : text_(std::move(text)) {}
    bool eof() { skip(); return position_ >= text_.size(); }
    bool consume_word(const std::string& expected) {
        skip();
        if (upper(text_.substr(position_, expected.size())) != expected) return false;
        const auto end = position_ + expected.size();
        if (end < text_.size() && (std::isalnum(static_cast<unsigned char>(text_[end])) || text_[end] == '_')) return false;
        position_ = end; return true;
    }
    void expect_word(const std::string& expected) {
        if (!consume_word(expected)) throw QueryError("Expected '" + expected + "' near: " + remaining());
    }
    bool consume(char value) { skip(); if (position_ < text_.size() && text_[position_] == value) { ++position_; return true; } return false; }
    std::string identifier() {
        skip(); const auto start = position_;
        while (position_ < text_.size() && (std::isalnum(static_cast<unsigned char>(text_[position_])) || text_[position_] == '_')) ++position_;
        if (start == position_) throw QueryError("Expected expression near: " + remaining());
        return upper(text_.substr(start, position_ - start));
    }
    std::string parenthesized_content() {
        skip(); if (!consume('(')) throw QueryError("Expected '(' near: " + remaining());
        const auto start = position_; int depth = 1;
        while (position_ < text_.size() && depth) {
            if (text_[position_] == '(') ++depth;
            else if (text_[position_] == ')') --depth;
            ++position_;
        }
        if (depth != 0) throw QueryError("Unclosed '(' in query");
        return text_.substr(start, position_ - start - 1);
    }
    std::optional<ast::ComparisonOp> comparison() {
        skip();
        const std::pair<const char*, ast::ComparisonOp> ops[] = {{">=", ast::ComparisonOp::GreaterEqual}, {"<=", ast::ComparisonOp::LessEqual}, {"!=", ast::ComparisonOp::NotEqual}, {"==", ast::ComparisonOp::Equal}, {">", ast::ComparisonOp::Greater}, {"<", ast::ComparisonOp::Less}};
        for (const auto& [token, op] : ops) { const std::string value(token); if (text_.substr(position_, value.size()) == value) { position_ += value.size(); return op; } }
        return std::nullopt;
    }
    std::string quantity(bool unit_required) {
        skip();
        const std::regex pattern(unit_required ? R"(^([+-]?[0-9]+(?:\.[0-9]+)?\s*[A-Za-z]+))" : R"(^([+-]?[0-9]+(?:\.[0-9]+)?(?:\s*[A-Za-z]+)?))");
        std::smatch match; const std::string rest = text_.substr(position_);
        if (!std::regex_search(rest, match, pattern)) throw QueryError("Expected numeric quantity near: " + remaining());
        position_ += match[1].length(); return trim(match[1].str());
    }
    std::string remaining() const { return trim(text_.substr(position_)); }
private:
    void skip() { while (position_ < text_.size() && std::isspace(static_cast<unsigned char>(text_[position_]))) ++position_; }
    std::string text_; std::size_t position_{};
};

ast::ExprPtr parse_expr(Cursor& cursor);

ast::ExprPtr parse_primary(Cursor& cursor) {
    ast::ExprPtr result;
    if (cursor.consume('(')) {
        result = parse_expr(cursor);
        if (!cursor.consume(')')) throw QueryError("Expected ')' after grouped expression");
    } else {
        const std::string name = cursor.identifier();
        const auto args = split_args(cursor.parenthesized_content());
        if (name == "CONTACT") {
            const auto options = contact_options(args, name, false);
            result = std::make_shared<ast::ContactExpr>(args[0], args[1], options.cutoff);
        } else if (name == "DISTANCE") {
            if (args.size() != 2) throw QueryError("DISTANCE requires exactly two selections");
            result = std::make_shared<ast::DistanceExpr>(args[0], args[1]);
        } else if (name == "CONTACT_COUNT") {
            const auto options = contact_options(args, name, true);
            result = std::make_shared<ast::ContactCountExpr>(args[0], args[1], options.cutoff, options.mode);
        } else if (name == "RG") {
            result = std::make_shared<ast::RgExpr>(args[0], rg_mass_weighted(args));
        } else {
            throw QueryError("Unknown observable: " + name);
        }
        if (auto op = cursor.comparison()) {
            if (name == "CONTACT") throw QueryError("CONTACT is already boolean and cannot be compared");
            const bool dimensionless = name == "CONTACT_COUNT";
            const std::string quantity = cursor.quantity(!dimensionless);
            double threshold{};
            if (dimensionless) {
                if (quantity.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz") != std::string::npos) throw QueryError("CONTACT_COUNT threshold must be unitless");
                threshold = std::stod(quantity);
            } else threshold = parse_distance(quantity).nm;
            result = std::make_shared<ast::ComparisonExpr>(result, *op, threshold, !dimensionless);
        }
    }
    if (cursor.consume_word("FOR")) {
        const auto op = cursor.comparison();
        if (!op || *op != ast::ComparisonOp::GreaterEqual) throw QueryError("FOR currently requires '>=' followed by a time quantity");
        result = std::make_shared<ast::ForExpr>(result, parse_duration(cursor.quantity(true)));
    }
    return result;
}

ast::ExprPtr parse_and(Cursor& cursor) {
    auto left = parse_primary(cursor);
    while (cursor.consume_word("AND")) left = std::make_shared<ast::AndExpr>(left, parse_primary(cursor));
    return left;
}

ast::ExprPtr parse_or(Cursor& cursor) {
    auto left = parse_and(cursor);
    while (cursor.consume_word("OR")) left = std::make_shared<ast::OrExpr>(left, parse_and(cursor));
    return left;
}

ast::ExprPtr parse_expr(Cursor& cursor) {
    auto left = parse_or(cursor);
    if (cursor.consume_word("FOLLOWED_BY")) {
        auto right = parse_or(cursor);
        std::optional<Duration> within;
        if (cursor.consume_word("WITHIN")) within = parse_duration(cursor.quantity(true));
        return std::make_shared<ast::FollowedByExpr>(left, right, within);
    }
    if (cursor.consume_word("OVERLAPS")) return std::make_shared<ast::OverlapsExpr>(left, parse_or(cursor));
    if (cursor.consume_word("BEFORE")) return std::make_shared<ast::BeforeExpr>(left, parse_or(cursor));
    if (cursor.consume_word("AFTER")) return std::make_shared<ast::AfterExpr>(left, parse_or(cursor));
    return left;
}
} // namespace

ast::Query Parser::parse(const std::string& text) const {
    try {
        Cursor cursor(text);
        cursor.expect_word("FIND");
        if (cursor.consume_word("EVENT")) cursor.consume_word("WHERE");
        auto root = parse_expr(cursor);
        cursor.consume(';');
        if (!cursor.eof()) {
            const auto rest = cursor.remaining();
            std::istringstream words(rest); std::string first; words >> first;
            if (!first.empty() && std::all_of(first.begin(), first.end(), [](unsigned char c) { return std::isalpha(c) || c == '_'; }))
                throw QueryError("Unknown temporal operator: " + first);
            throw QueryError("Unexpected query text: " + rest);
        }
        return {std::move(root)};
    } catch (const QueryError&) {
        throw;
    } catch (const std::exception& error) {
        throw QueryError(error.what());
    }
}

} // namespace ensembleql
