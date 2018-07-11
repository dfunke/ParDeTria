#include "Geometry.h"

#define GTEST_HAS_TR1_TUPLE 0

#include <gtest/gtest.h>
#include <bitset>

#include "utils/Generator.h"
#include "Partitioner.h"

TEST(Geometry2D, InSimplex) {

    dPoints<2, float> points;
    points.resize(7);

    const tIdType A = 0;
    const tIdType B = 1;
    const tIdType C = 2;

    // A (0,0)
    //points[A].id = A;
    points[A].coords = {{0, 0}};

    // B (2,0,0)
    //points[B].id = B;
    points[B].coords = {{1, 0}};

    // C (1,2,0)
    //points[C].id = C;
    points[C].coords = {{0, 1,}};

    dSimplex<2, float> s({{A, B, C}});

    EXPECT_TRUE(s.inSimplex(dPoint<2, float>({{.25, .25}}), points));
    EXPECT_TRUE(s.inSimplex(dPoint<2, float>({{.5, .5}}), points));
    EXPECT_FALSE(s.inSimplex(dPoint<2, float>({{1, 1}}), points));
}

TEST(Geometry2D, ClosestVertex) {

    dPoints<2, float> points;
    points.resize(7);

    const tIdType A = 0;
    const tIdType B = 1;
    const tIdType C = 2;

    points[A].coords = {{0, 0}};
    points[B].coords = {{1, 0}};
    points[C].coords = {{0, 1,}};

    dSimplex<2, float> s({{A, B, C}});

    EXPECT_EQ(s.findNearestVertex(dPoint<2, float>({{.25, .25}}), points), A);
    EXPECT_EQ(s.findNearestVertex(dPoint<2, float>({{.75, .25}}), points), B);
    EXPECT_EQ(s.findNearestVertex(dPoint<2, float>({{.25, .85}}), points), C);
    EXPECT_EQ(s.findNearestVertex(dPoint<2, float>({{.95, .75}}), points), B);
}

TEST(Geometry2D, NearestSimplex) {

    dPoints<2, float> points;
    points.resize(7);

    dSimplices<2, float> simplices(1, 8, 8);

    const tIdType A = 0;
    const tIdType B = 1;
    const tIdType C = 2;
    const tIdType D = 3;

    points[A].coords = {{0, 0}};
    points[B].coords = {{1, 0}};
    points[C].coords = {{0, 1,}};
    points[D].coords = {{1, 1,}};

    simplices.unsafe_at(1) = dSimplex<2,float>(1, {{A, B, C}}, {{2, dSimplex<2, float>::cINF, dSimplex<2, float>::cINF}});
    simplices.unsafe_at(2) = dSimplex<2,float>(2, {{B, C, D}}, {{dSimplex<2, float>::cINF, dSimplex<2, float>::cINF, 1}});

    EXPECT_EQ(simplices.findNearestSimplex(dPoint<2, float>({{.25, .25}}), points), 1);
    EXPECT_EQ(simplices.findNearestSimplex(dPoint<2, float>({{.75, .75}}), points), 2);
    EXPECT_EQ(simplices.findNearestSimplex(dPoint<2, float>({{1.25, 1.25}}), points), 2);
}