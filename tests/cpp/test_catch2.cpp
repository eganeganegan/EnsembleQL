#include "ensembleql/event.hpp"
#include "ensembleql/geometry.hpp"
#include "ensembleql/parser.hpp"
#include "ensembleql/temporal.hpp"
#include "ensembleql/units.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace ensembleql;

TEST_CASE("units are dimensional and normalized") {
    REQUIRE(parse_distance("10A").nm == 1.0);
    REQUIRE(parse_duration("5ns").ps == 5000.0);
    REQUIRE_THROWS_WITH(parse_distance("2ns"), "Distance expected, received time unit 'ns'");
}

TEST_CASE("geometry uses inclusive contact boundaries") {
    Frame frame{0.0, {{0.0, 0.0, 0.0}, {0.3, 0.4, 0.0}}, std::nullopt};
    REQUIRE(distance(frame.coordinates[0], frame.coordinates[1]) == 0.5);
    REQUIRE(contacts(frame, {0}, {1}, 0.5).size() == 1);
    REQUIRE(contacts(frame, {0}, {1}, 0.499).empty());
}

TEST_CASE("temporal joins include exact boundaries") {
    const Event left{0, 2, "left", {}, {}}, right{2, 4, "right", {}, {}};
    REQUIRE(temporal::before(left, right));
    REQUIRE(temporal::overlaps(left, right));
    REQUIRE(temporal::followed_by({left}, {right}, 0.0).size() == 1);
}

TEST_CASE("parser builds temporal AST independently of execution") {
    const auto query = Parser().parse(
        "FIND CONTACT(resid 17, resid 42) FOLLOWED_BY "
        "CONTACT(resid 17, resid 53) WITHIN 3ns;");
    REQUIRE(query.root->kind == ast::Kind::FollowedBy);
    REQUIRE_THROWS_AS(Parser().parse("FIND CONTACT(resid 17);"), QueryError);
}
