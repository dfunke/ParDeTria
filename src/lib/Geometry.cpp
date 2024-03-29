#include "Geometry.h"

#include <algorithm>
#include <sstream>
#include <atomic>

#include <tbb/parallel_for.h>
#include <tbb/spin_mutex.h>
#include <tbb/enumerable_thread_specific.h>

#include "Triangulator.h"
#include "utils/Predicates.h"

// static variables
template<uint D, typename Precision> constexpr tIdType dPoint<D, Precision>::cINF;

template<uint D, typename Precision>
constexpr tIdType dSimplex<D, Precision>::cINF;

template<uint D, typename Precision>
std::atomic<tIdType> dSimplices<D, Precision>::simplexID(1);
//

/* routines for orientation and inSphere/inCircle test based on
 * JR Shewchuk's implementation predicates.c */

/*****************************************************************************/
/*                                                                           */
/*  Routines for Arbitrary Precision Floating-point Arithmetic               */
/*  and Fast Robust Geometric Predicates                                     */
/*  (predicates.c)                                                           */
/*                                                                           */
/*  May 18, 1996                                                             */
/*                                                                           */
/*  Placed in the public domain by                                           */
/*  Jonathan Richard Shewchuk                                                */
/*  School of Computer Science                                               */
/*  Carnegie Mellon University                                               */
/*  5000 Forbes Avenue                                                       */
/*  Pittsburgh, Pennsylvania  15213-3891                                     */
/*  jrs@cs.cmu.edu                                                           */
/*****************************************************************************/

// common operations
template<std::size_t D, typename Precision>
Precision dot(const dVector<D, Precision> &a, const dVector<D, Precision> &b){
    Precision r = 0;
    for(uint d = 0; d < D; ++d)
        r += a[d] * b[d];
    return r;
}

template<std::size_t D, typename Precision>
dVector<D, Precision> operator+(const dVector<D, Precision> &a, const dVector<D, Precision> &b){
    dVector<D, Precision> r = a;
    for(uint d = 0; d < D; ++d)
        r[d] += b[d];
    return r;
}

template<std::size_t D, typename Precision>
dVector<D, Precision> operator-(const dVector<D, Precision> &a, const dVector<D, Precision> &b){
    dVector<D, Precision> r = b;
    for(uint d = 0; d < D; ++d)
        r[d] -= a[d];
    return r;
}

// specializations

template<typename Precision>
class GeometryCore<2, Precision> {
public:
    /*
       * return > 0: abc are counter-clockwise
       * return < 0: abc clockwise
       * return   0: co-linear
       *
       * Return a positive value if the points pa, pb, and pc occur
       * in counterclockwise order; a negative value if they occur
       * in clockwise order; and zero if they are collinear.  The
       * result is also a rough approximation of twice the signed
       * area of the triangle defined by the three points.
    */

    static Precision orientation(const dPoint<2, Precision> &s0, const dPoint<2, Precision> &s1, const dPoint<2, Precision> &s2) {
//        Precision acx, bcx, acy, bcy;
//
//        acx = s0.coords[0] -
//              s2.coords[0];
//        bcx = s1.coords[0] -
//              s2.coords[0];
//        acy = s0.coords[1] -
//              s2.coords[1];
//        bcy = s1.coords[1] -
//              s2.coords[1];
//        return acx * bcy - acy * bcx;

        return Predicates<2, Precision>::orient(s0.coords.data(), s1.coords.data(), s2.coords.data());
    }

    static bool inSimplex(const dPoint<2, Precision> &p,
                          const dPoint<2, Precision> &s0, const dPoint<2, Precision> &s1, const dPoint<2, Precision> &s2) {
        // Compute vectors
        auto v1 = s1.coords - s0.coords;
        auto v2 = s2.coords - s0.coords;
        auto vp = p.coords  - s0.coords;

        // Compute dot products
        auto v1v1 = dot(v1, v1);
        auto v1v2 = dot(v1, v2);
        auto v1vp = dot(v1, vp);
        auto v2v2 = dot(v2, v2);
        auto v2vp = dot(v2, vp);

        /*
         * Solve[{v1vp == u*(v1v1) + v*(v1v2),
         *        v2vp == u*(v1v2) + v*(v2v2)},
         *        {u, v}]
         *
         * {{u->-((v1vp v2v2-v1v2 v2vp)  / (v1v2^2-v1v1 v2v2)),
         *   v->-((-v1v2 v1vp+v1v1 v2vp) / (v1v2^2-v1v1 v2v2))}}
         *
         * */

        // Compute barycentric coordinates
        auto invDenom = Precision(1) / (v1v2 * v1v2 - v1v1 * v2v2);

        auto u = -(( v1vp * v2v2 - v1v2 * v2vp) * invDenom);
        auto v = -((-v1v2 * v1vp + v1v1 * v2vp) * invDenom);

        // Check if point is in triangle
        return (u >= 0) && (v >= 0) && (u + v <= 1);

    }

    static Precision sqDistanceToFace(const dPoint<2, Precision> &p,
                                      const dPoint<2, Precision> &s0, const dPoint<2, Precision> &s1) {

        auto v = s1.coords - s0.coords;
        auto r = p.coords  - s0.coords;

        auto n = dVector<2, Precision>({ v[1], -v[0] });

        return std::abs(dot(n,r) / dot(n, n));
    }

    static bool inSphere(const dPoint<2, Precision> &p,
                         const dPoint<2, Precision> &s0, const dPoint<2, Precision> &s1, const dPoint<2, Precision> &s2) {
//        Precision adx, ady, bdx, bdy, cdx, cdy;
//        Precision abdet, bcdet, cadet;
//        Precision alift, blift, clift;
//        Precision det;
//
//        adx = s0.coords[0] - p.coords[0];
//        ady = s0.coords[1] - p.coords[1];
//        bdx = s1.coords[0] - p.coords[0];
//        bdy = s1.coords[1] - p.coords[1];
//        cdx = s2.coords[0] - p.coords[0];
//        cdy = s2.coords[1] - p.coords[1];
//
//        abdet = adx * bdy - bdx * ady;
//        bcdet = bdx * cdy - cdx * bdy;
//        cadet = cdx * ady - adx * cdy;
//        alift = adx * adx + ady * ady;
//        blift = bdx * bdx + bdy * bdy;
//        clift = cdx * cdx + cdy * cdy;
//
//        /*
//         * Return a positive value if the point pd lies inside the
//         * circle passing through pa, pb, and pc; a negative value if
//         * it lies outside; and zero if the four points are cocircular.
//         * The points pa, pb, and pc must be in counterclockwise
//         * order, or the sign of the result will be reversed.
//         */
//
//        /*
//         * det > 0: d inside  abc (abc counter-clockwise)
//         *          d outside abc (abc clockwise)
//         * det < 0: d outside abc (abc counter-clockwise)
//         *          d inside  abc (abc clockwise)
//         */
//
//        det = alift * bcdet + blift * cadet + clift * abdet;

        return orientation(s0, s1, s2) * Predicates<2, Precision>::insphere(s0.coords.data(), s1.coords.data(), s2.coords.data(), p.coords.data()) >= 0;
    }

    static dSphere<2, Precision>
    circumsphere(const dPoint<2, Precision> &s0, const dPoint<2, Precision> &s1, const dPoint<2, Precision> &s2) {

        dSphere<2, Precision> sphere;

        Precision den = 2 * (s0.coords[0] *
                             (s1.coords[1] -
                              s2.coords[1]) +
                             s1.coords[0] *
                             (s2.coords[1] -
                              s0.coords[1]) +
                             s2.coords[0] *
                             (s0.coords[1] -
                              s1.coords[1]));

        sphere.center[0] = ((pow(s0.coords[0], 2) +
                             pow(s0.coords[1], 2)) *
                            (s1.coords[1] -
                             s2.coords[1]) +
                            (pow(s1.coords[0], 2) +
                             pow(s1.coords[1], 2)) *
                            (s2.coords[1] -
                             s0.coords[1]) +
                            (pow(s2.coords[0], 2) +
                             pow(s2.coords[1], 2)) *
                            (s0.coords[1] -
                             s1.coords[1])) /
                           den;

        sphere.center[1] = ((pow(s0.coords[0], 2) +
                             pow(s0.coords[1], 2)) *
                            (s2.coords[0] -
                             s1.coords[0]) +
                            (pow(s1.coords[0], 2) +
                             pow(s1.coords[1], 2)) *
                            (s0.coords[0] -
                             s2.coords[0]) +
                            (pow(s2.coords[0], 2) +
                             pow(s2.coords[1], 2)) *
                            (s1.coords[0] -
                             s0.coords[0])) /
                           den;

        sphere.radius =
                sqrt(pow(sphere.center[0] - s0.coords[0], 2) +
                     pow(sphere.center[1] - s0.coords[1], 2));

        return sphere;
    }
};

template<typename Precision>
class GeometryCore<3, Precision> {
public:
    /*
     * return > 0: d below abc - abc counter-clockwise
     *             d above abc - abc clockwise
     * return < 0: d above abc - abc counter-clockwise
     *             d below abc - abc clockwise
     * return   0: co-planar

     * Return a positive value if the point pd lies below the
     * plane passing through pa, pb, and pc; "below" is defined so
     * that pa, pb, and pc appear in counterclockwise order when
     * viewed from above the plane.  Returns a negative value if
     * pd lies above the plane.  Returns zero if the points are
     * coplanar.  The result is also a rough approximation of six
     * times the signed volume of the tetrahedron defined by the
     * four points.
     */

    static Precision orientation(const dPoint<3, Precision> &s0, const dPoint<3, Precision> &s1, const dPoint<3, Precision> &s2, const dPoint<3, Precision> &s3) {

//        Precision adx, bdx, cdx;
//        Precision ady, bdy, cdy;
//        Precision adz, bdz, cdz;
//
//        adx = s0.coords[0] -
//              s3.coords[0];
//        bdx = s1.coords[0] -
//              s3.coords[0];
//        cdx = s2.coords[0] -
//              s3.coords[0];
//        ady = s0.coords[1] -
//              s3.coords[1];
//        bdy = s1.coords[1] -
//              s3.coords[1];
//        cdy = s2.coords[1] -
//              s3.coords[1];
//        adz = s0.coords[2] -
//              s3.coords[2];
//        bdz = s1.coords[2] -
//              s3.coords[2];
//        cdz = s2.coords[2] -
//              s3.coords[2];
//
//        return adx * (bdy * cdz - bdz * cdy) + bdx * (cdy * adz - cdz * ady) +
//               cdx * (ady * bdz - adz * bdy);

        return Predicates<3, Precision>::orient(s0.coords.data(), s1.coords.data(), s2.coords.data(), s3.coords.data());
    }

    static bool inSimplex(const dPoint<3, Precision> &p,
                          const dPoint<3, Precision> &s0, const dPoint<3, Precision> &s1, const dPoint<3, Precision> &s2, const dPoint<3, Precision> &s3) {
        // Compute vectors
        auto v1 = s1.coords - s0.coords;
        auto v2 = s2.coords - s0.coords;
        auto v3 = s3.coords - s0.coords;
        auto vp = p.coords  - s0.coords;

        // Compute dot products
        auto v1v1 = dot(v1, v1);
        auto v1v2 = dot(v1, v2);
        auto v1v3 = dot(v1, v3);
        auto v1vp = dot(v1, vp);
        auto v2v2 = dot(v2, v2);
        auto v2v3 = dot(v2, v3);
        auto v2vp = dot(v2, vp);
        auto v3v3 = dot(v3, v3);
        auto v3vp = dot(v3, vp);

        /*
         * Solve[{v1vp == u*(v1v1) + v*(v1v2) + z*(v1v3),
                  v2vp == u*(v1v2) + v*(v2v2) + z*(v2v3),
                  v3vp == u*(v1v3) + v*(v2v3) + z*(v3v3)},
                  {u, v, z}]
         *
         * {{u -> -((-v1vp v2v3^2 + v1v3 v2v3 v2vp + v1vp v2v2 v3v3 - v1v2 v2vp v3v3 - v1v3 v2v2 v3vp + v1v2 v2v3 v3vp)
         *         /(v1v3^2 v2v2 - 2 v1v2 v1v3 v2v3 + v1v1 v2v3^2 + v1v2^2 v3v3 - v1v1 v2v2 v3v3)),
         *   v -> -((v1v3 v1vp v2v3 - v1v3^2 v2vp - v1v2 v1vp v3v3 + v1v1 v2vp v3v3 + v1v2 v1v3 v3vp - v1v1 v2v3 v3vp)
         *         /(v1v3^2 v2v2 - 2 v1v2 v1v3 v2v3 + v1v1 v2v3^2 + v1v2^2 v3v3 - v1v1 v2v2 v3v3)),
         *   z -> -((-v1v3 v1vp v2v2 + v1v2 v1vp v2v3 + v1v2 v1v3 v2vp - v1v1 v2v3 v2vp - v1v2^2 v3vp + v1v1 v2v2 v3vp)
         *         /(v1v3^2 v2v2 - 2 v1v2 v1v3 v2v3 + v1v1 v2v3^2 + v1v2^2 v3v3 - v1v1 v2v2 v3v3))}}
         *
         * */

        // Compute barycentric coordinates
        auto invDenom = Precision(1) / (v1v3 * v1v3 * v2v2 - 2 * v1v2 * v1v3 * v2v3 + v1v1 * v2v3 * v2v3 + v1v2 * v1v2 * v3v3 - v1v1 * v2v2 * v3v3);

        auto u = -((-v1vp * v2v3 * v2v3 + v1v3 * v2v3 * v2vp + v1vp * v2v2 * v3v3 - v1v2 * v2vp * v3v3 - v1v3 * v2v2 * v3vp + v1v2 * v2v3 * v3vp) * invDenom);
        auto v = -(( v1v3 * v1vp * v2v3 - v1v3 * v1v3 * v2vp - v1v2 * v1vp * v3v3 + v1v1 * v2vp * v3v3 + v1v2 * v1v3 * v3vp - v1v1 * v2v3 * v3vp) * invDenom);
        auto z = -((-v1v3 * v1vp * v2v2 + v1v2 * v1vp * v2v3 + v1v2 * v1v3 * v2vp - v1v1 * v2v3 * v2vp - v1v2 * v1v2 * v3vp + v1v1 * v2v2 * v3vp) * invDenom);

        // Check if point is in triangle
        return (u >= 0) && (v >= 0) && (z >= 0) && (u + v + z <= 1);

    }

    static Precision sqDistanceToFace(const dPoint<3, Precision> &p,
                                      const dPoint<3, Precision> &s0, const dPoint<3, Precision> &s1, const dPoint<3, Precision> &s2) {

        auto u = s1.coords - s0.coords;
        auto v = s2.coords - s0.coords;
        auto w = p.coords  - s0.coords;

        auto n = dVector<3, Precision>({ (u[1] * v[2] - u[2] * v[1]),
                                         (u[2] * v[0] - u[0] * v[2]),
                                         (u[0] * v[1] - u[1] * v[0]) });

        return std::abs(dot(w,n) / dot(n,n));
    }

    static bool inSphere(const dPoint<3, Precision> &p,
                         const dPoint<3, Precision> &s0, const dPoint<3, Precision> &s1, const dPoint<3, Precision> &s2, const dPoint<3, Precision> &s3) {

//        Precision aex, bex, cex, dex;
//        Precision aey, bey, cey, dey;
//        Precision aez, bez, cez, dez;
//        Precision alift, blift, clift, dlift;
//        Precision ab, bc, cd, da, ac, bd;
//        Precision abc, bcd, cda, dab;
//        Precision det;
//
//        aex = s0.coords[0] - p.coords[0];
//        bex = s1.coords[0] - p.coords[0];
//        cex = s2.coords[0] - p.coords[0];
//        dex = s3.coords[0] - p.coords[0];
//        aey = s0.coords[1] - p.coords[1];
//        bey = s1.coords[1] - p.coords[1];
//        cey = s2.coords[1] - p.coords[1];
//        dey = s3.coords[1] - p.coords[1];
//        aez = s0.coords[2] - p.coords[2];
//        bez = s1.coords[2] - p.coords[2];
//        cez = s2.coords[2] - p.coords[2];
//        dez = s3.coords[2] - p.coords[2];
//
//        ab = aex * bey - bex * aey;
//        bc = bex * cey - cex * bey;
//        cd = cex * dey - dex * cey;
//        da = dex * aey - aex * dey;
//
//        ac = aex * cey - cex * aey;
//        bd = bex * dey - dex * bey;
//
//        abc = aez * bc - bez * ac + cez * ab;
//        bcd = bez * cd - cez * bd + dez * bc;
//        cda = cez * da + dez * ac + aez * cd;
//        dab = dez * ab + aez * bd + bez * da;
//
//        alift = aex * aex + aey * aey + aez * aez;
//        blift = bex * bex + bey * bey + bez * bez;
//        clift = cex * cex + cey * cey + cez * cez;
//        dlift = dex * dex + dey * dey + dez * dez;
//
//        /*
//         * det > 0: e inside  abcd (abcd positive orientation)
//         *          e outside abcd (abcd negative orientation)
//         * det < 0: e outside abcd (abcd positive orientation)
//         *          e inside  abcd (abcd negative orientation)
//         *
//         * Return a positive value if the point pe lies inside the
//         * sphere passing through pa, pb, pc, and pd; a negative value
//         * if it lies outside; and zero if the five points are
//         * cospherical.  The points pa, pb, pc, and pd must be ordered
//         * so that they have a positive orientation (as defined by
//         * orient3d()), or the sign of the result will be reversed.
//         */
//
//        det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);

        return orientation(s0, s1, s2, s3) * Predicates<3, Precision>::insphere(s0.coords.data(), s1.coords.data(), s2.coords.data(), s3.coords.data(), p.coords.data()) >= 0;
    }

    static dSphere<3, Precision>
    circumsphere(const dPoint<3, Precision> &s0, const dPoint<3, Precision> &s1, const dPoint<3, Precision> &s2, const dPoint<3, Precision> &s3) {

        dSphere<3, Precision> sphere;

        Precision den = 2 * (s0.coords[0] *
                             (s2.coords[1] *
                              s3.coords[2] +
                              s1.coords[1] *
                              (s2.coords[2] -
                               s3.coords[2]) -
                              s3.coords[1] *
                              s2.coords[2] -
                              (s2.coords[1] -
                               s3.coords[1]) *
                              s1.coords[2]) -
                             s1.coords[0] *
                             (s2.coords[1] *
                              s3.coords[2] -
                              s3.coords[1] *
                              s2.coords[2]) -
                             s0.coords[1] *
                             (s2.coords[0] *
                              s3.coords[2] +
                              s1.coords[0] *
                              (s2.coords[2] -
                               s3.coords[2]) -
                              s3.coords[0] *
                              s2.coords[2] -
                              (s2.coords[0] -
                               s3.coords[0]) *
                              s1.coords[2]) +
                             s1.coords[1] *
                             (s2.coords[0] *
                              s3.coords[2] -
                              s3.coords[0] *
                              s2.coords[2]) +
                             s0.coords[2] *
                             (s2.coords[0] *
                              s3.coords[1] +
                              s1.coords[0] *
                              (s2.coords[1] -
                               s3.coords[1]) -
                              s3.coords[0] *
                              s2.coords[1] -
                              (s2.coords[0] -
                               s3.coords[0]) *
                              s1.coords[1]) -
                             s1.coords[2] *
                             (s2.coords[0] *
                              s3.coords[1] -
                              s3.coords[0] *
                              s2.coords[1]));

        sphere.center[0] =
                (-s0.coords[1] *
                 (-s2.coords[2] *
                  (pow(s3.coords[2], 2) +
                   pow(s3.coords[1], 2) +
                   pow(s3.coords[0], 2)) -
                  s1.coords[2] *
                  (-pow(s3.coords[2], 2) +
                   pow(s2.coords[2], 2) -
                   pow(s3.coords[1], 2) +
                   pow(s2.coords[1], 2) -
                   pow(s3.coords[0], 2) +
                   pow(s2.coords[0], 2)) +
                  (pow(s2.coords[2], 2) +
                   pow(s2.coords[1], 2) +
                   pow(s2.coords[0], 2)) *
                  s3.coords[2] +
                  (pow(s1.coords[2], 2) +
                   pow(s1.coords[1], 2) +
                   pow(s1.coords[0], 2)) *
                  (s2.coords[2] -
                   s3.coords[2])) +
                 s1.coords[1] *
                 ((pow(s2.coords[2], 2) +
                   pow(s2.coords[1], 2) +
                   pow(s2.coords[0], 2)) *
                  s3.coords[2] -
                  s2.coords[2] *
                  (pow(s3.coords[2], 2) +
                   pow(s3.coords[1], 2) +
                   pow(s3.coords[0], 2))) +
                 s0.coords[2] *
                 (-s2.coords[1] *
                  (pow(s3.coords[2], 2) +
                   pow(s3.coords[1], 2) +
                   pow(s3.coords[0], 2)) -
                  s1.coords[1] *
                  (-pow(s3.coords[2], 2) +
                   pow(s2.coords[2], 2) -
                   pow(s3.coords[1], 2) +
                   pow(s2.coords[1], 2) -
                   pow(s3.coords[0], 2) +
                   pow(s2.coords[0], 2)) +
                  s3.coords[1] *
                  (pow(s2.coords[2], 2) +
                   pow(s2.coords[1], 2) +
                   pow(s2.coords[0], 2)) +
                  (s2.coords[1] -
                   s3.coords[1]) *
                  (pow(s1.coords[2], 2) +
                   pow(s1.coords[1], 2) +
                   pow(s1.coords[0], 2))) -
                 s1.coords[2] *
                 (s3.coords[1] *
                  (pow(s2.coords[2], 2) +
                   pow(s2.coords[1], 2) +
                   pow(s2.coords[0], 2)) -
                  s2.coords[1] *
                  (pow(s3.coords[2], 2) +
                   pow(s3.coords[1], 2) +
                   pow(s3.coords[0], 2))) +
                 (pow(s0.coords[2], 2) +
                  pow(s0.coords[1], 2) +
                  pow(s0.coords[0], 2)) *
                 (s2.coords[1] *
                  s3.coords[2] +
                  s1.coords[1] *
                  (s2.coords[2] -
                   s3.coords[2]) -
                  s3.coords[1] *
                  s2.coords[2] -
                  (s2.coords[1] -
                   s3.coords[1]) *
                  s1.coords[2]) -
                 (pow(s1.coords[2], 2) +
                  pow(s1.coords[1], 2) +
                  pow(s1.coords[0], 2)) *
                 (s2.coords[1] *
                  s3.coords[2] -
                  s3.coords[1] *
                  s2.coords[2])) /
                den;

        sphere.center[1] =
                (s0.coords[0] *
                 (-s2.coords[2] *
                  (pow(s3.coords[2], 2) +
                   pow(s3.coords[1], 2) +
                   pow(s3.coords[0], 2)) -
                  s1.coords[2] *
                  (-pow(s3.coords[2], 2) +
                   pow(s2.coords[2], 2) -
                   pow(s3.coords[1], 2) +
                   pow(s2.coords[1], 2) -
                   pow(s3.coords[0], 2) +
                   pow(s2.coords[0], 2)) +
                  (pow(s2.coords[2], 2) +
                   pow(s2.coords[1], 2) +
                   pow(s2.coords[0], 2)) *
                  s3.coords[2] +
                  (pow(s1.coords[2], 2) +
                   pow(s1.coords[1], 2) +
                   pow(s1.coords[0], 2)) *
                  (s2.coords[2] -
                   s3.coords[2])) -
                 s1.coords[0] *
                 ((pow(s2.coords[2], 2) +
                   pow(s2.coords[1], 2) +
                   pow(s2.coords[0], 2)) *
                  s3.coords[2] -
                  s2.coords[2] *
                  (pow(s3.coords[2], 2) +
                   pow(s3.coords[1], 2) +
                   pow(s3.coords[0], 2))) -
                 s0.coords[2] *
                 (-s2.coords[0] *
                  (pow(s3.coords[2], 2) +
                   pow(s3.coords[1], 2) +
                   pow(s3.coords[0], 2)) -
                  s1.coords[0] *
                  (-pow(s3.coords[2], 2) +
                   pow(s2.coords[2], 2) -
                   pow(s3.coords[1], 2) +
                   pow(s2.coords[1], 2) -
                   pow(s3.coords[0], 2) +
                   pow(s2.coords[0], 2)) +
                  s3.coords[0] *
                  (pow(s2.coords[2], 2) +
                   pow(s2.coords[1], 2) +
                   pow(s2.coords[0], 2)) +
                  (s2.coords[0] -
                   s3.coords[0]) *
                  (pow(s1.coords[2], 2) +
                   pow(s1.coords[1], 2) +
                   pow(s1.coords[0], 2))) +
                 s1.coords[2] *
                 (s3.coords[0] *
                  (pow(s2.coords[2], 2) +
                   pow(s2.coords[1], 2) +
                   pow(s2.coords[0], 2)) -
                  s2.coords[0] *
                  (pow(s3.coords[2], 2) +
                   pow(s3.coords[1], 2) +
                   pow(s3.coords[0], 2))) -
                 (pow(s0.coords[2], 2) +
                  pow(s0.coords[1], 2) +
                  pow(s0.coords[0], 2)) *
                 (s2.coords[0] *
                  s3.coords[2] +
                  s1.coords[0] *
                  (s2.coords[2] -
                   s3.coords[2]) -
                  s3.coords[0] *
                  s2.coords[2] -
                  (s2.coords[0] -
                   s3.coords[0]) *
                  s1.coords[2]) +
                 (pow(s1.coords[2], 2) +
                  pow(s1.coords[1], 2) +
                  pow(s1.coords[0], 2)) *
                 (s2.coords[0] *
                  s3.coords[2] -
                  s3.coords[0] *
                  s2.coords[2])) /
                den;

        sphere.center[2] =
                (-s0.coords[0] *
                 (-s2.coords[1] *
                  (pow(s3.coords[2], 2) +
                   pow(s3.coords[1], 2) +
                   pow(s3.coords[0], 2)) -
                  s1.coords[1] *
                  (-pow(s3.coords[2], 2) +
                   pow(s2.coords[2], 2) -
                   pow(s3.coords[1], 2) +
                   pow(s2.coords[1], 2) -
                   pow(s3.coords[0], 2) +
                   pow(s2.coords[0], 2)) +
                  s3.coords[1] *
                  (pow(s2.coords[2], 2) +
                   pow(s2.coords[1], 2) +
                   pow(s2.coords[0], 2)) +
                  (s2.coords[1] -
                   s3.coords[1]) *
                  (pow(s1.coords[2], 2) +
                   pow(s1.coords[1], 2) +
                   pow(s1.coords[0], 2))) +
                 s1.coords[0] *
                 (s3.coords[1] *
                  (pow(s2.coords[2], 2) +
                   pow(s2.coords[1], 2) +
                   pow(s2.coords[0], 2)) -
                  s2.coords[1] *
                  (pow(s3.coords[2], 2) +
                   pow(s3.coords[1], 2) +
                   pow(s3.coords[0], 2))) +
                 s0.coords[1] *
                 (-s2.coords[0] *
                  (pow(s3.coords[2], 2) +
                   pow(s3.coords[1], 2) +
                   pow(s3.coords[0], 2)) -
                  s1.coords[0] *
                  (-pow(s3.coords[2], 2) +
                   pow(s2.coords[2], 2) -
                   pow(s3.coords[1], 2) +
                   pow(s2.coords[1], 2) -
                   pow(s3.coords[0], 2) +
                   pow(s2.coords[0], 2)) +
                  s3.coords[0] *
                  (pow(s2.coords[2], 2) +
                   pow(s2.coords[1], 2) +
                   pow(s2.coords[0], 2)) +
                  (s2.coords[0] -
                   s3.coords[0]) *
                  (pow(s1.coords[2], 2) +
                   pow(s1.coords[1], 2) +
                   pow(s1.coords[0], 2))) -
                 s1.coords[1] *
                 (s3.coords[0] *
                  (pow(s2.coords[2], 2) +
                   pow(s2.coords[1], 2) +
                   pow(s2.coords[0], 2)) -
                  s2.coords[0] *
                  (pow(s3.coords[2], 2) +
                   pow(s3.coords[1], 2) +
                   pow(s3.coords[0], 2))) -
                 (s2.coords[0] *
                  s3.coords[1] -
                  s3.coords[0] *
                  s2.coords[1]) *
                 (pow(s1.coords[2], 2) +
                  pow(s1.coords[1], 2) +
                  pow(s1.coords[0], 2)) +
                 (s2.coords[0] *
                  s3.coords[1] +
                  s1.coords[0] *
                  (s2.coords[1] -
                   s3.coords[1]) -
                  s3.coords[0] *
                  s2.coords[1] -
                  (s2.coords[0] -
                   s3.coords[0]) *
                  s1.coords[1]) *
                 (pow(s0.coords[2], 2) +
                  pow(s0.coords[1], 2) +
                  pow(s0.coords[0], 2))) /
                den;

        sphere.radius =
                sqrt(pow(sphere.center[0] - s0.coords[0], 2) +
                     pow(sphere.center[1] - s0.coords[1], 2) +
                     pow(sphere.center[2] - s0.coords[2], 2));

        return sphere;
    }
};

//template<uint D, typename Precision>
//tIdType dSimplices<D, Precision>::countDuplicates(const Simplex_Ids & simplices) const {
//
//    std::atomic<tIdType> duplicates(0);
//    tbb::spin_mutex mtx;
//
//    tbb::enumerable_thread_specific<dSimplicesConstHandle<D, Precision>,
//            tbb::cache_aligned_allocator<dSimplicesConstHandle<D, Precision>>,
//            tbb::ets_key_usage_type::ets_key_per_instance> tsThisHandle(std::ref(*this));
//
//    tbb::parallel_for(simplices.range(), [&](const auto &r) {
//
//        auto thisHandle = tsThisHandle.local();
//
//        for (const auto &a : r) {
//
//            auto s = thisHandle.at(a);
//            auto b = simplices.begin();
//            b.setIdx(a - simplices.lowerBound() + 1, true);
//            for(; b < simplices.end(); ++b){
//                if(s.equalVertices(thisHandle.at(*b))){
//                    // we found a duplicate
//                    ++duplicates;
//
//                    tbb::spin_mutex::scoped_lock lock(mtx);
//                    LOG("found duplicate " << s << " and " << thisHandle.at(*b) << std::endl);
//                }
//            }
//        }
//    });
//
//    return duplicates;
//}

template<uint D, typename Precision>
CrossCheckReport<D, Precision> dSimplices<D, Precision>::crossCheck(
        const dSimplices<D, Precision> &realSimplices) const {
    CrossCheckReport<D, Precision> result;
    result.valid = true;

    std::atomic<tIdType> mySize(0);
    std::atomic<tIdType> realSize(0);

    //build hashmap for simplex lookup
    tbb::concurrent_unordered_multimap<tHashType, tIdType> simplexLookup(this->exact_size() + realSimplices.exact_size());

    tbb::parallel_for(this->range(), [&](const auto &r) {

        for (const auto &a : r) {
            ++mySize;
            simplexLookup.insert(std::make_pair(a.fingerprint(), a.id));
        }
    });

    tbb::parallel_for(realSimplices.range(), [&](const auto & r) {

        for (const auto &a : r) {
            ++realSize;
            simplexLookup.insert(std::make_pair(a.fingerprint(), a.id));
        }
    });

    // check whether sizes are equal
    if (mySize != realSize) {
        LOG("Wrong number of simplices - my size: " << mySize
            << " real size: " << realSize << std::endl);
        result.valid = false;
    }

    tbb::spin_mutex mtx;

    tbb::enumerable_thread_specific<dSimplicesConstHandle<D, Precision>,
            tbb::cache_aligned_allocator<dSimplicesConstHandle<D, Precision>>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsThisHandle(std::ref(*this));

    tbb::enumerable_thread_specific<dSimplicesConstHandle<D, Precision>,
            tbb::cache_aligned_allocator<dSimplicesConstHandle<D, Precision>>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsRealHandle(std::ref(realSimplices));

    // check whether all simplices of real DT are present
    tbb::parallel_for(realSimplices.range(), [&](const auto &r) {

        auto &thisHandle = tsThisHandle.local();
        auto &realHandle = tsRealHandle.local();

        for(const auto & realSimplex : r) {
            // find corresponding mySimplex for realSimplex
            // we can limit the search by using the face where-used ds
            auto hash = realSimplex.fingerprint();
            auto range = simplexLookup.equal_range(hash);
            auto mySimplex = std::find_if(range.first,
                                          range.second,
                                          [&](const auto &i) {
                                              return dSimplex<D, Precision>::isFinite(i.second)
                                                     && thisHandle.contains(i.second)
                                                     && thisHandle.at(i.second).equalVertices(realSimplex);
                                          });

            if (mySimplex == range.second) {
                tbb::spin_mutex::scoped_lock lock(mtx);
                LOG("did not find simplex " << realSimplex << std::endl);
                result.valid = false;
                result.missing.push_back(realSimplex);
                return;
            }

            // check neighbors
            for (const auto &n : realSimplex.neighbors) {
                if (!dSimplex<D, Precision>::isFinite(n))
                    continue;

                bool found = false;
                for (const auto &nn : thisHandle.at(mySimplex->second).neighbors) {
                    if (dSimplex<D, Precision>::isFinite(nn) &&
                        thisHandle.at(nn).equalVertices(realHandle[n])) {
                        found = true;
                        break;
                    }
                }

                if (!found) {
                    tbb::spin_mutex::scoped_lock lock(mtx);
                    LOG("did not find neighbor " << realHandle[n] << " of simplex "
                                                                 << realSimplex << std::endl);
                    result.valid = false;
                }
            }
        }
    });

    // check for own simplices that are not in real DT
    tbb::parallel_for(this->range(), [&](const auto &r) {

        auto &realHandle = tsRealHandle.local();

        for (const auto &mySimplex : r) {

            // again we use the face where-used ds for the lookup
            auto hash = mySimplex.fingerprint();
            auto range = simplexLookup.equal_range(hash);
            auto realSimplex = std::find_if(range.first,
                                            range.second,
                                            [&](const auto &i) {
                                                return dSimplex<D, Precision>::isFinite(i.second)
                                                       && realHandle.contains(i.second)
                                                       && realHandle.at(i.second).equalVertices(mySimplex);
                                            });

            if (realSimplex == range.second /*&& mySimplex.isFinite()*/) {

                tbb::spin_mutex::scoped_lock lock(mtx);
                LOG("simplex " << mySimplex << " does not exist in real triangulation"
                    << std::endl);
                result.valid = false;
                result.invalid.push_back(mySimplex);
            }
        }
    });

    // check whether sizes are equal
    //TODO no size operator yet - do we need it?
//    if (pt.simplices.size() != realPT.simplices.size()) {
//        LOG("my size: " + std::to_string(pt.simplices.size()) +
//            " real size: "
//            << std::to_string(realPT.simplices.size()) + "\n");
//        result.valid = false;
//    }

    LOG("cross check " << (result.valid ? "" : "NOT ") << "successful"
        << std::endl);

    return result;
}

template<uint D, typename Precision>
VerificationReport<D, Precision>
dSimplices<D, Precision>::verify(const dPoints<D, Precision> &points, const Point_Ids *pointIds) const {
    INDENT
    VerificationReport<D, Precision> result;
    result.valid = true;

    tbb::spin_mutex mtx;

    tbb::enumerable_thread_specific<dSimplicesConstHandle<D, Precision>,
            tbb::cache_aligned_allocator<dSimplicesConstHandle<D, Precision>>,
            tbb::ets_key_usage_type::ets_key_per_instance> tsThisHandle(std::ref(*this));

    // verify that every input point is used
    LOG("Checking points" << std::endl);
    std::unordered_set<tIdType> usedPoints;
    for (const auto &s : *this) {
        usedPoints.insert(s.vertices.begin(), s.vertices.end());
    }
    if (points != usedPoints) {
        // not all points of input used
        std::stringstream sNotUsed;

        if (pointIds) {
            for (const auto i : *pointIds) {
                const auto &p = points[i];
                if (usedPoints.count(i) != 1 && dPoint<D, Precision>::isFinite(i)) {
                    sNotUsed << i << " - " << p << " ";
                    result.valid = false;
                }
            }
        } else {
            for (std::size_t i = 0; i < points.finite_size(); ++i) {
                const auto &p = points[i];
                if (usedPoints.count(i) != 1 && dPoint<D, Precision>::isFinite(i)) {
                    sNotUsed << i << " - " << p << " ";
                    result.valid = false;
                }
            }
        }
        if (!result.valid) {
            LOG("Points of input not used: " << sNotUsed.str() << std::endl);
        }

        std::stringstream sInvalidP;
        if (pointIds) {
            for (const auto &p : usedPoints) {
                if (dPoint<D, Precision>::isFinite(p) && !pointIds->contains(p)) {
                    sInvalidP << p << " ";
                    result.valid = false;
                }
            }
        } else {
            for (const auto &p : usedPoints) {
                if (dPoint<D, Precision>::isFinite(p) && !points.contains(p)) {
                    sInvalidP << p << " ";
                    result.valid = false;
                }
            }
        }
        if (!result.valid) {
            LOG("Used points not in input: " << sInvalidP.str() << std::endl);
        }
    }

    // verify convex hull
    LOG("Checking convex-hull" << std::endl);
    tbb::parallel_for(this->range(), [&](const auto &r) {

        for (const auto &s : r) {

            if (!s.isFinite() && std::find(convexHull.begin(), convexHull.end(), s.id) == convexHull.end()) {
                // s is infinite but not part of the convex hull
                tbb::spin_mutex::scoped_lock lock(mtx);
                LOG("Infinite simplex " << s << " NOT in convex hull" << std::endl);
                result.valid = false;
            }
        }
    });

    tbb::parallel_for(tbb::blocked_range<std::size_t>(std::size_t(0), convexHull.size(), 1000), [&](const auto & r) {

        auto &thisHandle = tsThisHandle.local();

        for (auto i = r.begin(); i != r.end(); ++i) {

            const dSimplex<D, Precision> &s = thisHandle.at(convexHull[i]);
            if (s.isFinite()) {
                // s is finite but part of convex hull
                tbb::spin_mutex::scoped_lock lock(mtx);
                LOG("Finite simplex " << s << " IS in convex hull" << std::endl);
                result.valid = false;
            }
        }
    });

    //check for duplicates
    //build hashmap for simplex lookup
    LOG("Building simplex look-up map" << std::endl);
    tbb::concurrent_unordered_multimap<tHashType, tIdType> simplexLookup(this->exact_size());

    tbb::parallel_for(this->range(), [&](const auto &r) {

        for (const auto &a : r) {
            auto fp = a.fingerprint();
            simplexLookup.insert(std::make_pair(fp, a.id));
        }
    });

    LOG("Check for duplicates" << std::endl);
    tbb::parallel_for(this->range(), [&](const auto &r) {

        auto &thisHandle = tsThisHandle.local();

        for (const auto &a : r) {

            //check for any duplicates in the triangulation
            auto simplexHash = a.fingerprint();
            auto range = simplexLookup.equal_range(simplexHash);

            ASSERT(range.first != range.second); // we have at least ourself in the map

            for (auto it = range.first; it != range.second; ++it) {

                const auto &b = thisHandle.at(it->second);
                // a and b must have identical vertices and same id or different id and different vertices
                if (!(a.equalVertices(b) == (a.id == b.id))) {

                    tbb::spin_mutex::scoped_lock lock(mtx);
                    LOG("Found duplicates " << a << " and " << b
                        << std::endl);
                    result.duplicates.emplace_back(a, b);
                    result.valid = false;
                }
            }
        }
    });


    // verify that all simplices with a shared D-1 simplex are neighbors

    //build hashmap for face lookup
    LOG("Building facette look-up map" << std::endl);
    tbb::concurrent_unordered_multimap<tHashType, tIdType> faceLookup((D+1) * this->exact_size());

    tbb::parallel_for(this->range(), [&](const auto &r) {

        for (const auto &a : r) {
            for (uint i = 0; i < D + 1; ++i) {
                auto facetteHash = a.faceFingerprint(i);

                faceLookup.insert(std::make_pair(facetteHash, a.id));
            }
        }
    });

    LOG("Checking neighbors" << std::endl);
    tbb::parallel_for(this->range(), [&](const auto &r) {

        auto &thisHandle = tsThisHandle.local();

        for (const auto &a : r) {

            // we already verified that the face where-used data structure is correct
            // we can use it to verify the neighbor relation
            for (uint i = 0; i < D + 1; ++i) {
                auto faceHash = a.faceFingerprint(i);
                auto range = faceLookup.equal_range(faceHash);

                ASSERT(range.first != range.second); // we have at least ourself in the map

                for (auto it = range.first; it != range.second; ++it) {

                    const auto &b = thisHandle.at(it->second);
                    // a and b are neighbors: the neighbor property is symmetric and the
                    // corresponding simplices must be present in the neighbors arrays
                    // accordingly
                    if ((a.isNeighbor(b) &&
                         (!b.isNeighbor(a) ||
                          std::find(a.neighbors.begin(), a.neighbors.end(), b.id) == a.neighbors.end() ||
                          std::find(b.neighbors.begin(), b.neighbors.end(), a.id) == b.neighbors.end()))
                        // a and b are NOT neighbors: must be symmetric and simplices NOT be
                        // present in neighbors arrays
                        ||
                        (!a.isNeighbor(b) &&
                         (b.isNeighbor(a) ||
                          std::find(a.neighbors.begin(), a.neighbors.end(), b.id) != a.neighbors.end() ||
                          std::find(b.neighbors.begin(), b.neighbors.end(), a.id) != b.neighbors.end()))) {

                        tbb::spin_mutex::scoped_lock lock(mtx);
                        LOG("Wrong neighbor relation between " << a << " and " << b
                            << std::endl);
                        result.wrongNeighbors.emplace_back(a, b);
                        result.valid = false;
                    }
                }
            }
        }
    });

    // check in circle criterion
    LOG("Checking empty circle criterion" << std::endl);
    tbb::parallel_for(this->range(), [&](const auto &r) {

        auto &thisHandle = tsThisHandle.local();

        for (const auto &s : r) {

            // we have established that the neighboorhood is correctly set
            // we check for all neighbors whether the NOT shared point is in the circle

            for (const auto &n : s.neighbors) {
                if (!dSimplex<D, Precision>::isFinite(n))
                    continue;

                const dSimplex<D, Precision> &nn = thisHandle.at(n);
                for (uint d = 0; d < D + 1; ++d) {
                    if (dPoint<D,Precision>::isFinite(nn.vertices[d]) && !s.contains(nn.vertices[d])) {
                        // we have found the point of nn that is NOT shared with s
                        const auto &p = points[nn.vertices[d]];
                        if (s.inSphere(p, points)) {
                            LOG("Point " << nn.vertices[d] << " " << p << " is in circle of " << s << std::endl);

                            tbb::spin_mutex::scoped_lock lock(mtx);
                            result.valid = false;

                            result.inCircle[s].insert(nn.vertices[d]);
                        }
                    }
                }
            }
        }
    });
    DEDENT

    LOG("Triangulation is " << (result.valid ? "" : "NOT ") << "valid"
        << std::endl);

    return result;
}

template <uint D, typename Precision>
typename dSimplices<D, Precision>::tHash dSimplices<D, Precision>::genFingerprint() const {

    tHash hash;
    std::fill(hash.begin(), hash.end(), 0); //init with zero

    for(auto  s : *this){
        auto h = s.fingerprint();
        uint i = h & 0xFF; // get last byte

        ASSERT(i <= 256);
        hash[i] ^= h; //xor into hash array
    }

    return hash;

}

// specialiations

// float

template
class dPoint<2, float>;

template
class dPoint<3, float>;

template
class dSimplex<2, float>;

template
class dSimplex<3, float>;

template
class dSimplices<2, float>;

template
class dSimplices<3, float>;

template
class GeometryCore<2, float>;

template
class GeometryCore<3, float>;

// double

template
class dPoint<2, double>;

template
class dPoint<3, double>;

template
class dSimplex<2, double>;

template
class dSimplex<3, double>;

template
class dSimplices<2, double>;

template
class dSimplices<3, double>;

template
class GeometryCore<2, double>;

template
class GeometryCore<3, double>;
