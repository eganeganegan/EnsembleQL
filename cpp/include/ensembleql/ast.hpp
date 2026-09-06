#pragma once

#include "ensembleql/topology.hpp"
#include "ensembleql/units.hpp"

#include <memory>
#include <optional>
#include <string>

namespace ensembleql::ast {

enum class Kind { Contact, Distance, ContactCount, Rg, Comparison, For, FollowedBy, Overlaps, Before, After, And, Or };

struct Expr {
    explicit Expr(Kind node_kind) : kind(node_kind) {}
    virtual ~Expr() = default;
    Kind kind;
};

using ExprPtr = std::shared_ptr<Expr>;

struct ContactExpr final : Expr {
    ContactExpr(std::string first, std::string second, Distance distance = {0.45})
        : Expr(Kind::Contact), a(std::move(first)), b(std::move(second)), cutoff(distance) {}
    std::string a, b;
    Distance cutoff;
};

struct DistanceExpr final : Expr {
    DistanceExpr(std::string first, std::string second)
        : Expr(Kind::Distance), a(std::move(first)), b(std::move(second)) {}
    std::string a, b;
};

struct ContactCountExpr final : Expr {
    ContactCountExpr(std::string first, std::string second, Distance distance = {0.45},
                     ContactMode contact_mode = ContactMode::Atom)
        : Expr(Kind::ContactCount), a(std::move(first)), b(std::move(second)),
          cutoff(distance), mode(contact_mode) {}
    std::string a, b;
    Distance cutoff;
    ContactMode mode;
};

struct RgExpr final : Expr {
    explicit RgExpr(std::string value, bool weighted = false)
        : Expr(Kind::Rg), selection(std::move(value)), mass_weighted(weighted) {}
    std::string selection;
    bool mass_weighted;
};

enum class ComparisonOp { Less, LessEqual, Greater, GreaterEqual, Equal, NotEqual };
struct ComparisonExpr final : Expr {
    ComparisonExpr(ExprPtr value, ComparisonOp operation, double threshold_value, bool is_distance)
        : Expr(Kind::Comparison), operand(std::move(value)), op(operation), threshold(threshold_value), distance_threshold(is_distance) {}
    ExprPtr operand;
    ComparisonOp op;
    double threshold;
    bool distance_threshold;
};

struct ForExpr final : Expr {
    ForExpr(ExprPtr value, Duration minimum) : Expr(Kind::For), operand(std::move(value)), duration(minimum) {}
    ExprPtr operand;
    Duration duration;
};

struct BinaryExpr : Expr {
    BinaryExpr(Kind node_kind, ExprPtr lhs, ExprPtr rhs)
        : Expr(node_kind), left(std::move(lhs)), right(std::move(rhs)) {}
    ExprPtr left, right;
};

struct FollowedByExpr final : BinaryExpr {
    FollowedByExpr(ExprPtr lhs, ExprPtr rhs, std::optional<Duration> limit)
        : BinaryExpr(Kind::FollowedBy, std::move(lhs), std::move(rhs)), within(limit) {}
    std::optional<Duration> within;
};
struct OverlapsExpr final : BinaryExpr { OverlapsExpr(ExprPtr a, ExprPtr b) : BinaryExpr(Kind::Overlaps, std::move(a), std::move(b)) {} };
struct BeforeExpr final : BinaryExpr { BeforeExpr(ExprPtr a, ExprPtr b) : BinaryExpr(Kind::Before, std::move(a), std::move(b)) {} };
struct AfterExpr final : BinaryExpr { AfterExpr(ExprPtr a, ExprPtr b) : BinaryExpr(Kind::After, std::move(a), std::move(b)) {} };
struct AndExpr final : BinaryExpr { AndExpr(ExprPtr a, ExprPtr b) : BinaryExpr(Kind::And, std::move(a), std::move(b)) {} };
struct OrExpr final : BinaryExpr { OrExpr(ExprPtr a, ExprPtr b) : BinaryExpr(Kind::Or, std::move(a), std::move(b)) {} };

struct Query { ExprPtr root; };

} // namespace ensembleql::ast
