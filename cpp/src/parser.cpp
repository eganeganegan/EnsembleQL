#include "ensembleql/parser.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
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
                               bool allow_mode, Distance default_cutoff = {0.45}) {
    if (args.size() < 2) throw QueryError(observable + " requires two selections");
    ContactOptions options;
    options.cutoff = default_cutoff;
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

struct DistanceAngleOptions {
    Distance distance;
    Angle angle;
};

DistanceAngleOptions distance_angle_options(const std::vector<std::string>& args,
                                            const std::string& observable,
                                            Distance default_distance,
                                            Angle default_angle,
                                            const std::string& angle_key) {
    if (args.size() < 2) throw QueryError(observable + " requires two selections");
    DistanceAngleOptions result{default_distance, default_angle};
    bool distance_seen = false, angle_seen = false;
    for (std::size_t index = 2; index < args.size(); ++index) {
        const auto equals = args[index].find('=');
        if (equals == std::string::npos) throw QueryError("Expected named " + observable + " option");
        const std::string key = upper(trim(args[index].substr(0, equals)));
        const std::string value = trim(args[index].substr(equals + 1));
        if (key == "DISTANCE" || key == "CUTOFF") {
            if (distance_seen) throw QueryError("Duplicate " + observable + " distance option");
            result.distance = parse_distance(value); distance_seen = true;
        } else if (key == upper(angle_key) || key == "ANGLE") {
            if (angle_seen) throw QueryError("Duplicate " + observable + " angle option");
            result.angle = parse_angle(value); angle_seen = true;
        } else {
            throw QueryError("Unknown " + observable + " option: " + trim(args[index].substr(0, equals)));
        }
    }
    if (result.distance.nm < 0.0 || result.angle.degrees < 0.0 || result.angle.degrees > 180.0) {
        throw QueryError("Invalid " + observable + " distance or angle criterion");
    }
    return result;
}

bool boolean_option(const std::vector<std::string>& args, std::size_t index,
                    const std::string& observable, const std::string& expected_key) {
    const auto equals = args[index].find('=');
    if (equals == std::string::npos || upper(trim(args[index].substr(0, equals))) != expected_key) {
        throw QueryError("Unknown " + observable + " option: " + trim(args[index].substr(0, equals)));
    }
    const std::string value = upper(trim(args[index].substr(equals + 1)));
    if (value == "TRUE") return true;
    if (value == "FALSE") return false;
    throw QueryError(observable + " " + trim(args[index].substr(0, equals)) + " option must be true or false");
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
        const std::regex pattern(unit_required ? R"(^([+-]?[0-9]+(?:\.[0-9]+)?\s*[A-Za-z][A-Za-z0-9]*))" : R"(^([+-]?[0-9]+(?:\.[0-9]+)?))");
        std::smatch match; const std::string rest = text_.substr(position_);
        if (!std::regex_search(rest, match, pattern)) throw QueryError("Expected numeric quantity near: " + remaining());
        position_ += match[1].length(); return trim(match[1].str());
    }
    std::string integer() {
        skip();
        const auto start = position_;
        while (position_ < text_.size() &&
               std::isdigit(static_cast<unsigned char>(text_[position_]))) ++position_;
        if (start == position_) throw QueryError("Expected positive integer near: " + remaining());
        return text_.substr(start, position_ - start);
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
        } else if (name == "HBOND") {
            const auto options = distance_angle_options(args, name, {0.35}, {150.0}, "MIN_ANGLE");
            result = std::make_shared<ast::HydrogenBondExpr>(args[0], args[1], options.distance, options.angle);
        } else if (name == "DIHEDRAL") {
            if (args.size() != 4) throw QueryError("DIHEDRAL requires exactly four single-atom selections");
            result = std::make_shared<ast::DihedralExpr>(args[0], args[1], args[2], args[3]);
        } else if (name == "HELIX") {
            if (args.empty() || args.size() > 2) throw QueryError("HELIX requires a selection and optional minimum_fraction");
            double fraction = 0.5;
            if (args.size() == 2) {
                const auto equals = args[1].find('=');
                if (equals == std::string::npos || upper(trim(args[1].substr(0, equals))) != "MINIMUM_FRACTION") {
                    throw QueryError("Unknown HELIX option: " + trim(args[1].substr(0, equals)));
                }
                const std::string value = trim(args[1].substr(equals + 1));
                std::size_t consumed{};
                fraction = std::stod(value, &consumed);
                if (consumed != value.size()) throw QueryError("HELIX minimum_fraction must be numeric");
                if (!std::isfinite(fraction) || fraction < 0.0 || fraction > 1.0) {
                    throw QueryError("HELIX minimum_fraction must be between 0 and 1");
                }
            }
            result = std::make_shared<ast::HelixExpr>(args[0], fraction);
        } else if (name == "SASA") {
            if (args.empty() || args.size() > 3) throw QueryError("SASA requires a selection and optional probe/points options");
            Distance probe{0.14}; std::size_t points = 96;
            bool probe_seen = false, points_seen = false;
            for (std::size_t index = 1; index < args.size(); ++index) {
                const auto equals = args[index].find('=');
                if (equals == std::string::npos) throw QueryError("Expected named SASA option");
                const std::string key = upper(trim(args[index].substr(0, equals)));
                const std::string value = trim(args[index].substr(equals + 1));
                if (key == "PROBE") {
                    if (probe_seen) throw QueryError("Duplicate SASA probe option");
                    probe = parse_distance(value); probe_seen = true;
                }
                else if (key == "POINTS") {
                    if (points_seen) throw QueryError("Duplicate SASA points option");
                    std::size_t consumed{};
                    const auto parsed = std::stoull(value, &consumed);
                    if (consumed != value.size()) throw QueryError("SASA points must be a positive integer");
                    points = static_cast<std::size_t>(parsed);
                    points_seen = true;
                } else throw QueryError("Unknown SASA option: " + key);
            }
            if (probe.nm < 0.0 || points < 6) throw QueryError("SASA requires a non-negative probe and at least six points");
            result = std::make_shared<ast::SasaExpr>(args[0], probe, points);
        } else if (name == "RMSD") {
            if (args.empty() || args.size() > 2) throw QueryError("RMSD requires a selection and optional align option");
            result = std::make_shared<ast::RmsdExpr>(args[0], args.size() == 1 ? true : boolean_option(args, 1, name, "ALIGN"));
        } else if (name == "COORDINATION_NUMBER") {
            const auto options = contact_options(args, name, false, {0.35});
            result = std::make_shared<ast::CoordinationNumberExpr>(args[0], args[1], options.cutoff);
        } else if (name == "SALT_BRIDGE") {
            const auto options = contact_options(args, name, false, {0.4});
            result = std::make_shared<ast::SaltBridgeExpr>(args[0], args[1], options.cutoff);
        } else if (name == "AROMATIC_STACKING") {
            const auto options = distance_angle_options(args, name, {0.55}, {30.0}, "MAX_ANGLE");
            if (options.angle.degrees > 90.0) throw QueryError("AROMATIC_STACKING maximum angle must not exceed 90 degrees");
            result = std::make_shared<ast::AromaticStackingExpr>(args[0], args[1], options.distance, options.angle);
        } else if (name == "SURFACE_DISTANCE") {
            if (args.size() != 2) throw QueryError("SURFACE_DISTANCE requires exactly two selections");
            result = std::make_shared<ast::SurfaceDistanceExpr>(args[0], args[1]);
        } else if (name == "ORIENTATION") {
            if (args.size() != 3) throw QueryError("ORIENTATION requires two selections and axis=x|y|z");
            const auto equals = args[2].find('=');
            if (equals == std::string::npos || upper(trim(args[2].substr(0, equals))) != "AXIS") {
                throw QueryError("ORIENTATION third argument must be axis=x|y|z");
            }
            const std::string axis = upper(trim(args[2].substr(equals + 1)));
            if (axis.size() != 1 || (axis[0] != 'X' && axis[0] != 'Y' && axis[0] != 'Z')) {
                throw QueryError("ORIENTATION axis must be x, y, or z");
            }
            result = std::make_shared<ast::OrientationExpr>(args[0], args[1], axis[0]);
        } else {
            throw QueryError("Unknown observable: " + name);
        }
        if (auto op = cursor.comparison()) {
            if (name == "CONTACT" || name == "HBOND" || name == "HELIX" ||
                name == "SALT_BRIDGE" || name == "AROMATIC_STACKING") {
                throw QueryError(name + " is already boolean and cannot be compared");
            }
            ast::ValueDimension dimension = ast::ValueDimension::Distance;
            if (name == "CONTACT_COUNT" || name == "COORDINATION_NUMBER") dimension = ast::ValueDimension::Unitless;
            else if (name == "DIHEDRAL" || name == "ORIENTATION") dimension = ast::ValueDimension::Angle;
            else if (name == "SASA") dimension = ast::ValueDimension::Area;
            const std::string quantity = cursor.quantity(dimension != ast::ValueDimension::Unitless);
            double threshold{};
            if (dimension == ast::ValueDimension::Unitless) {
                threshold = std::stod(quantity);
            } else if (dimension == ast::ValueDimension::Distance) threshold = parse_distance(quantity).nm;
            else if (dimension == ast::ValueDimension::Angle) threshold = parse_angle(quantity).degrees;
            else threshold = parse_area(quantity).nm2;
            result = std::make_shared<ast::ComparisonExpr>(result, *op, threshold, dimension);
        }
    }
    if (cursor.consume_word("FOR")) {
        const auto op = cursor.comparison();
        if (!op || *op != ast::ComparisonOp::GreaterEqual) throw QueryError("FOR currently requires '>=' followed by a time quantity");
        result = std::make_shared<ast::ForExpr>(result, parse_duration(cursor.quantity(true)));
    }
    if (cursor.consume_word("REPEATS")) {
        const auto op = cursor.comparison();
        if (!op || *op != ast::ComparisonOp::GreaterEqual) {
            throw QueryError("REPEATS requires '>=' followed by a positive integer count");
        }
        const std::string quantity = cursor.integer();
        std::size_t count{};
        try {
            std::size_t consumed{};
            count = static_cast<std::size_t>(std::stoull(quantity, &consumed));
            if (consumed != quantity.size() || count == 0) throw std::invalid_argument("count");
        } catch (const std::exception&) {
            throw QueryError("REPEATS count must be a unitless positive integer");
        }
        std::optional<Duration> within;
        if (cursor.consume_word("WITHIN")) within = parse_duration(cursor.quantity(true));
        result = std::make_shared<ast::RepeatsExpr>(result, count, within);
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
    while (true) {
        if (cursor.consume_word("FOLLOWED_BY")) {
            auto right = parse_or(cursor);
            std::optional<Duration> within;
            if (cursor.consume_word("WITHIN")) within = parse_duration(cursor.quantity(true));
            left = std::make_shared<ast::FollowedByExpr>(left, right, within);
        } else if (cursor.consume_word("IMMEDIATELY_FOLLOWED_BY")) {
            left = std::make_shared<ast::ImmediatelyFollowedByExpr>(left, parse_or(cursor));
        } else if (cursor.consume_word("OVERLAPS")) {
            left = std::make_shared<ast::OverlapsExpr>(left, parse_or(cursor));
        } else if (cursor.consume_word("DURING")) {
            left = std::make_shared<ast::DuringExpr>(left, parse_or(cursor));
        } else if (cursor.consume_word("UNTIL")) {
            left = std::make_shared<ast::UntilExpr>(left, parse_or(cursor));
        } else if (cursor.consume_word("BEFORE")) {
            left = std::make_shared<ast::BeforeExpr>(left, parse_or(cursor));
        } else if (cursor.consume_word("PRECEDES")) {
            left = std::make_shared<ast::PrecedesExpr>(left, parse_or(cursor));
        } else if (cursor.consume_word("AFTER")) {
            left = std::make_shared<ast::AfterExpr>(left, parse_or(cursor));
        } else {
            break;
        }
    }
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
