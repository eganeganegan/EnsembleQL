#pragma once

#include "ensembleql/topology.hpp"
#include "ensembleql/units.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>

namespace ensembleql::ast {

enum class Kind {
    Contact, Distance, ContactCount, Rg, HydrogenBond, Dihedral, Helix, Sasa,
    Rmsd, CoordinationNumber, SaltBridge, AromaticStacking, SurfaceDistance,
    Orientation, Comparison, For, Repeats,
    FollowedBy, ImmediatelyFollowedBy, Overlaps, During, Until, Before,
    Precedes, After, And, Or
};

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

struct HydrogenBondExpr final : Expr {
    HydrogenBondExpr(std::string donor_selection, std::string acceptor_selection,
                     Distance distance = {0.35}, Angle angle = {150.0})
        : Expr(Kind::HydrogenBond), donors(std::move(donor_selection)),
          acceptors(std::move(acceptor_selection)), distance_cutoff(distance),
          minimum_angle(angle) {}
    std::string donors, acceptors;
    Distance distance_cutoff;
    Angle minimum_angle;
};

struct DihedralExpr final : Expr {
    DihedralExpr(std::string a, std::string b, std::string c, std::string d)
        : Expr(Kind::Dihedral), selections{std::move(a), std::move(b),
                                          std::move(c), std::move(d)} {}
    std::array<std::string, 4> selections;
};

struct HelixExpr final : Expr {
    HelixExpr(std::string value, double fraction = 0.5)
        : Expr(Kind::Helix), selection(std::move(value)), minimum_fraction(fraction) {}
    std::string selection;
    double minimum_fraction;
};

struct SasaExpr final : Expr {
    SasaExpr(std::string value, Distance probe_radius = {0.14}, std::size_t point_count = 96)
        : Expr(Kind::Sasa), selection(std::move(value)), probe(probe_radius),
          points(point_count) {}
    std::string selection;
    Distance probe;
    std::size_t points;
};

struct RmsdExpr final : Expr {
    RmsdExpr(std::string value, bool should_align = true)
        : Expr(Kind::Rmsd), selection(std::move(value)), align(should_align) {}
    std::string selection;
    bool align;
};

struct CoordinationNumberExpr final : Expr {
    CoordinationNumberExpr(std::string center_selection, std::string neighbor_selection,
                           Distance distance = {0.35})
        : Expr(Kind::CoordinationNumber), centers(std::move(center_selection)),
          neighbors(std::move(neighbor_selection)), cutoff(distance) {}
    std::string centers, neighbors;
    Distance cutoff;
};

struct SaltBridgeExpr final : Expr {
    SaltBridgeExpr(std::string positive_selection, std::string negative_selection,
                   Distance distance = {0.4})
        : Expr(Kind::SaltBridge), positive(std::move(positive_selection)),
          negative(std::move(negative_selection)), cutoff(distance) {}
    std::string positive, negative;
    Distance cutoff;
};

struct AromaticStackingExpr final : Expr {
    AromaticStackingExpr(std::string first, std::string second,
                         Distance distance = {0.55}, Angle angle = {30.0})
        : Expr(Kind::AromaticStacking), first_ring(std::move(first)),
          second_ring(std::move(second)), centroid_cutoff(distance),
          maximum_angle(angle) {}
    std::string first_ring, second_ring;
    Distance centroid_cutoff;
    Angle maximum_angle;
};

struct SurfaceDistanceExpr final : Expr {
    SurfaceDistanceExpr(std::string molecule_selection, std::string surface_selection)
        : Expr(Kind::SurfaceDistance), molecule(std::move(molecule_selection)),
          surface(std::move(surface_selection)) {}
    std::string molecule, surface;
};

struct OrientationExpr final : Expr {
    OrientationExpr(std::string origin_selection, std::string target_selection,
                    char reference_axis = 'z')
        : Expr(Kind::Orientation), origin(std::move(origin_selection)),
          target(std::move(target_selection)), axis(reference_axis) {}
    std::string origin, target;
    char axis;
};

enum class ComparisonOp { Less, LessEqual, Greater, GreaterEqual, Equal, NotEqual };
enum class ValueDimension { Unitless, Distance, Angle, Area };
struct ComparisonExpr final : Expr {
    ComparisonExpr(ExprPtr value, ComparisonOp operation, double threshold_value,
                   ValueDimension value_dimension)
        : Expr(Kind::Comparison), operand(std::move(value)), op(operation),
          threshold(threshold_value), dimension(value_dimension) {}
    ExprPtr operand;
    ComparisonOp op;
    double threshold;
    ValueDimension dimension;
};

struct ForExpr final : Expr {
    ForExpr(ExprPtr value, Duration minimum) : Expr(Kind::For), operand(std::move(value)), duration(minimum) {}
    ExprPtr operand;
    Duration duration;
};

struct RepeatsExpr final : Expr {
    RepeatsExpr(ExprPtr value, std::size_t minimum_count,
                std::optional<Duration> limit = std::nullopt)
        : Expr(Kind::Repeats), operand(std::move(value)), count(minimum_count), within(limit) {}
    ExprPtr operand;
    std::size_t count;
    std::optional<Duration> within;
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
struct ImmediatelyFollowedByExpr final : BinaryExpr {
    ImmediatelyFollowedByExpr(ExprPtr a, ExprPtr b)
        : BinaryExpr(Kind::ImmediatelyFollowedBy, std::move(a), std::move(b)) {}
};
struct OverlapsExpr final : BinaryExpr { OverlapsExpr(ExprPtr a, ExprPtr b) : BinaryExpr(Kind::Overlaps, std::move(a), std::move(b)) {} };
struct DuringExpr final : BinaryExpr { DuringExpr(ExprPtr a, ExprPtr b) : BinaryExpr(Kind::During, std::move(a), std::move(b)) {} };
struct UntilExpr final : BinaryExpr { UntilExpr(ExprPtr a, ExprPtr b) : BinaryExpr(Kind::Until, std::move(a), std::move(b)) {} };
struct BeforeExpr final : BinaryExpr { BeforeExpr(ExprPtr a, ExprPtr b) : BinaryExpr(Kind::Before, std::move(a), std::move(b)) {} };
struct PrecedesExpr final : BinaryExpr { PrecedesExpr(ExprPtr a, ExprPtr b) : BinaryExpr(Kind::Precedes, std::move(a), std::move(b)) {} };
struct AfterExpr final : BinaryExpr { AfterExpr(ExprPtr a, ExprPtr b) : BinaryExpr(Kind::After, std::move(a), std::move(b)) {} };
struct AndExpr final : BinaryExpr { AndExpr(ExprPtr a, ExprPtr b) : BinaryExpr(Kind::And, std::move(a), std::move(b)) {} };
struct OrExpr final : BinaryExpr { OrExpr(ExprPtr a, ExprPtr b) : BinaryExpr(Kind::Or, std::move(a), std::move(b)) {} };

struct Query { ExprPtr root; };

} // namespace ensembleql::ast
