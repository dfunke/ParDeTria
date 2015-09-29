#include "Geometry.h"

#define GTEST_HAS_TR1_TUPLE 0

#include <gtest/gtest.h>
#include <bitset>

#include "utils/Generator.h"
#include "Partitioner.h"

TEST(Geometry3D, HashingTest) {

    dSimplex<3, double> s;
    s.vertices = {{2553, 3528, 4184, 6413}};
    s.genFingerprint();

    dSimplex<3, double> eq;
    eq.vertices = {{2553, 3528, 4184, 6413}}; // equal simplex
    eq.genFingerprint();

    dSimplex<3, double> pe;
    pe.vertices = {{2553, 4184, 3528, 6413}}; // permuted simplex
    pe.genFingerprint();

    dSimplex<3, double> nb;
    nb.vertices = {{3528, 4184, 5555, 6413}}; // neighbor simplex
    nb.genFingerprint();

    dSimplex<3, double> ne;
    ne.vertices = {{2543, 3268, 4345, 6093}}; // not equal simplex
    ne.genFingerprint();

    EXPECT_EQ(s.fingerprint(), eq.fingerprint());

    // fingerprint should be equal
    EXPECT_EQ(s.fingerprint(), pe.genFingerprint());

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
    EXPECT_NE(s.fingerprint(), ne.fingerprint());

}

TEST(Geometry3D, ZeroHash) {

    //typedef std::bitset<32> bs;

    dSimplex<3, double> s;
    s.vertices = {{97, 892, 1371, 3932}};
    s.genFingerprint();

    EXPECT_NE(s.fingerprint(), 0);
    //std::cout << "Fingerprint: " << s.fingerprint() << " - " << bs(s.fingerprint()).to_string() << std::endl;

    for(uint d = 0; d < 3 + 1; ++d) {
        EXPECT_NE(s.faceFingerprint(d), 0);
        //std::cout << "Face " << d << ": " << s.faceFingerprint(d) << " - " << bs(s.faceFingerprint(d)).to_string() << std::endl;
    }

}

TEST(Geometry3D, SimpleInSphere) {

    dPoints<3, float> points;
    points.resize(7);

    const tIdType A = 0;
    const tIdType B = 1;
    const tIdType C = 2;

    const tIdType Du = 3;
    const tIdType Dd = 4;

    const tIdType I = 5;
    const tIdType O = 6;

    // A (0,0,0)
    //points[A].id = A;
    points[A].coords = {{0, 0, 0}};

    // B (2,0,0)
    //points[B].id = B;
    points[B].coords = {{2, 0, 0}};

    // C (1,2,0)
    //points[C].id = C;
    points[C].coords = {{1, 2, 0}};

    // D_u (1,1,1)
    //points[Du].id = Du;
    points[Du].coords = {{1, 1, 1}};

    // D_d (1,1,-1)
    //points[Dd].id = Dd;
    points[Dd].coords = {{1, 1, -1}};

    // I (1,1,0)
    //points[I].id = I;
    points[I].coords = {{1, 1, 0}};

    // O (5,5,5)
    //points[O].id = O;
    points[O].coords = {{5, 5, 5}};

    auto test = [&](const std::array<tIdType, 4> &v) {
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

TEST(Geometry3D, PointStats) {

    dBox<3, double> bounds(
            dVector<3, double>( {{ 0, 0, 0 }}),
            dVector<3, double>({{ 100,100,100}})
    );

    tGenerator gen(START_SEED);
    auto dice = RandomFactory<double>::make('u', gen);
    auto points = genPoints(1e5, bounds, dice);

    dPointStats<3, double> seq = getPointStatsSeq(std::size_t(0), points.size(), points);
    dPointStats<3, double> par = getPointStats(points, points);

    for(uint d = 0; d < 3; ++d){
        EXPECT_EQ(seq.min[d], par.min[d]);
        EXPECT_EQ(seq.mid[d], par.mid[d]);
        EXPECT_EQ(seq.max[d], par.max[d]);
    }
}