#include "Geometry.h"

#define GTEST_HAS_TR1_TUPLE 0

#include <gtest/gtest.h>

TEST(Geometry3D, HashingTest) {

    dSimplex<3, double> s;
    s.vertices = {{2553, 3528, 4184, 6413}};
    s.fingerprint();

    dSimplex<3, double> eq;
    eq.vertices = {{2553, 3528, 4184, 6413}}; // equal simplex
    eq.fingerprint();

    dSimplex<3, double> pe;
    pe.vertices = {{2553, 4184, 3528, 6413}}; // permuted simplex
    pe.fingerprint();

    dSimplex<3, double> nb;
    nb.vertices = {{3528, 4184, 5555, 6413}}; // neighbor simplex
    nb.fingerprint();

    dSimplex<3, double> ne;
    ne.vertices = {{2543, 3268, 4345, 6093}}; // not equal simplex
    ne.fingerprint();

    EXPECT_EQ(s.vertexFingerprint, eq.vertexFingerprint);

    // fingerprint should be equal
    EXPECT_EQ(s.vertexFingerprint, pe.fingerprint());

    // find neighbor
    bool foundNeighbor = false;
    for (uint i = 0; i < 4; ++i) {
        for (uint j = 0; j < 4; ++j) {
            tHashType a = s.faceFingerprint(i);
            tHashType b = nb.faceFingerprint(j);
            if (a == b)
                foundNeighbor = true;
        }
    }
    EXPECT_TRUE(foundNeighbor);

    // the other simplex
    EXPECT_NE(s.vertexFingerprint, ne.vertexFingerprint);

}

TEST(Geometry3D, SimpleInSphere) {

    dPoints<3, float> points;
    points.resize(7);

    const uint A = 0;
    const uint B = 1;
    const uint C = 2;

    const uint Du = 3;
    const uint Dd = 4;

    const uint I = 5;
    const uint O = 6;

    // A (0,0,0)
    points[A].id = A;
    points[A].coords = {{0, 0, 0}};

    // B (2,0,0)
    points[B].id = B;
    points[B].coords = {{2, 0, 0}};

    // C (1,2,0)
    points[C].id = C;
    points[C].coords = {{1, 2, 0}};

    // D_u (1,1,1)
    points[Du].id = Du;
    points[Du].coords = {{1, 1, 1}};

    // D_d (1,1,-1)
    points[Dd].id = Dd;
    points[Dd].coords = {{1, 1, -1}};

    // I (1,1,0)
    points[I].id = I;
    points[I].coords = {{1, 1, 0}};

    // O (5,5,5)
    points[O].id = O;
    points[O].coords = {{5, 5, 5}};

    auto test = [&](const std::array<uint, 4> &v) {
        dSimplex<3, float> s;
        s.vertices = v;

        EXPECT_TRUE(s.inSphere(points[I], points));
        EXPECT_FALSE(s.inSphere(points[O], points));

    };

    // triangle ABC is counter-clockwise
    // D_d is below -> orientation > 0
    test({{A, B, C, Dd}});

    // triangle ABC is counter-clockwise
    // D_u is above -> orientation < 0
    test({{A, B, C, Du}});

    // triangle ACB is clockwise
    // D_d is below -> orientation < 0
    test({{A, C, B, Dd}});

    // triangle ACB is clockwise
    // D_u is above -> orientation > 0
    test({{A, C, B, Du}});
}
