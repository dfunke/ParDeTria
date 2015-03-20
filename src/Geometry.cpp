#include "Geometry.h"

#include <algorithm>
#include <sstream>
#include <atomic>

#include <tbb/parallel_for.h>

// static variables
template <uint D, typename Precision> constexpr uint dPoint<D, Precision>::cINF;

template <uint D, typename Precision>
constexpr uint dSimplex<D, Precision>::cINF;
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

template <uint D, typename Precision> class GeometryHelper {
public:
  static Precision orientation(const dSimplex<D, Precision> &simplex,
                               const dPoints<D, Precision> &points);
  static bool inSphere(const dSimplex<D, Precision> &simplex,
                       const dPoint<D, Precision> &p,
                       const dPoints<D, Precision> &points);
  static dSphere<D, Precision>
  circumsphere(const dSimplex<D, Precision> &simplex,
               const dPoints<D, Precision> &points);
};

// specializations

template <typename Precision> class GeometryHelper<2, Precision> {
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

  static Precision orientation(const dSimplex<2, Precision> &simplex,
                               const dPoints<2, Precision> &points) {
    Precision acx, bcx, acy, bcy;

    acx = points[simplex.vertices[0]].coords[0] -
          points[simplex.vertices[2]].coords[0];
    bcx = points[simplex.vertices[1]].coords[0] -
          points[simplex.vertices[2]].coords[0];
    acy = points[simplex.vertices[0]].coords[1] -
          points[simplex.vertices[2]].coords[1];
    bcy = points[simplex.vertices[1]].coords[1] -
          points[simplex.vertices[2]].coords[1];
    return acx * bcy - acy * bcx;
  }

  static bool inSphere(const dSimplex<2, Precision> &simplex,
                       const dPoint<2, Precision> &p,
                       const dPoints<2, Precision> &points) {
    Precision adx, ady, bdx, bdy, cdx, cdy;
    Precision abdet, bcdet, cadet;
    Precision alift, blift, clift;
    Precision det;

    adx = points[simplex.vertices[0]].coords[0] - p.coords[0];
    ady = points[simplex.vertices[0]].coords[1] - p.coords[1];
    bdx = points[simplex.vertices[1]].coords[0] - p.coords[0];
    bdy = points[simplex.vertices[1]].coords[1] - p.coords[1];
    cdx = points[simplex.vertices[2]].coords[0] - p.coords[0];
    cdy = points[simplex.vertices[2]].coords[1] - p.coords[1];

    abdet = adx * bdy - bdx * ady;
    bcdet = bdx * cdy - cdx * bdy;
    cadet = cdx * ady - adx * cdy;
    alift = adx * adx + ady * ady;
    blift = bdx * bdx + bdy * bdy;
    clift = cdx * cdx + cdy * cdy;

    /*
     * Return a positive value if the point pd lies inside the
     * circle passing through pa, pb, and pc; a negative value if
     * it lies outside; and zero if the four points are cocircular.
     * The points pa, pb, and pc must be in counterclockwise
     * order, or the sign of the result will be reversed.
     */

    /*
     * det > 0: d inside  abc (abc counter-clockwise)
     *          d outside abc (abc clockwise)
     * det < 0: d outside abc (abc counter-clockwise)
     *          d inside  abc (abc clockwise)
     */

    det = alift * bcdet + blift * cadet + clift * abdet;

    return orientation(simplex, points) * det >= 0;
  }

  static dSphere<2, Precision>
  circumsphere(const dSimplex<2, Precision> &simplex,
               const dPoints<2, Precision> &points) {

    dSphere<2, Precision> sphere;

    Precision den = 2 * (points[simplex.vertices[0]].coords[0] *
                             (points[simplex.vertices[1]].coords[1] -
                              points[simplex.vertices[2]].coords[1]) +
                         points[simplex.vertices[1]].coords[0] *
                             (points[simplex.vertices[2]].coords[1] -
                              points[simplex.vertices[0]].coords[1]) +
                         points[simplex.vertices[2]].coords[0] *
                             (points[simplex.vertices[0]].coords[1] -
                              points[simplex.vertices[1]].coords[1]));

    sphere.center[0] = ((pow(points[simplex.vertices[0]].coords[0], 2) +
                         pow(points[simplex.vertices[0]].coords[1], 2)) *
                            (points[simplex.vertices[1]].coords[1] -
                             points[simplex.vertices[2]].coords[1]) +
                        (pow(points[simplex.vertices[1]].coords[0], 2) +
                         pow(points[simplex.vertices[1]].coords[1], 2)) *
                            (points[simplex.vertices[2]].coords[1] -
                             points[simplex.vertices[0]].coords[1]) +
                        (pow(points[simplex.vertices[2]].coords[0], 2) +
                         pow(points[simplex.vertices[2]].coords[1], 2)) *
                            (points[simplex.vertices[0]].coords[1] -
                             points[simplex.vertices[1]].coords[1])) /
                       den;

    sphere.center[1] = ((pow(points[simplex.vertices[0]].coords[0], 2) +
                         pow(points[simplex.vertices[0]].coords[1], 2)) *
                            (points[simplex.vertices[2]].coords[0] -
                             points[simplex.vertices[1]].coords[0]) +
                        (pow(points[simplex.vertices[1]].coords[0], 2) +
                         pow(points[simplex.vertices[1]].coords[1], 2)) *
                            (points[simplex.vertices[0]].coords[0] -
                             points[simplex.vertices[2]].coords[0]) +
                        (pow(points[simplex.vertices[2]].coords[0], 2) +
                         pow(points[simplex.vertices[2]].coords[1], 2)) *
                            (points[simplex.vertices[1]].coords[0] -
                             points[simplex.vertices[0]].coords[0])) /
                       den;

    sphere.radius =
        sqrt(pow(sphere.center[0] - points[simplex.vertices[0]].coords[0], 2) +
             pow(sphere.center[1] - points[simplex.vertices[0]].coords[1], 2));

    return sphere;
  }
};

template <typename Precision> class GeometryHelper<3, Precision> {
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

  static Precision orientation(const dSimplex<3, Precision> &simplex,
                               const dPoints<3, Precision> &points) {

    Precision adx, bdx, cdx;
    Precision ady, bdy, cdy;
    Precision adz, bdz, cdz;

    adx = points[simplex.vertices[0]].coords[0] -
          points[simplex.vertices[3]].coords[0];
    bdx = points[simplex.vertices[1]].coords[0] -
          points[simplex.vertices[3]].coords[0];
    cdx = points[simplex.vertices[2]].coords[0] -
          points[simplex.vertices[3]].coords[0];
    ady = points[simplex.vertices[0]].coords[1] -
          points[simplex.vertices[3]].coords[1];
    bdy = points[simplex.vertices[1]].coords[1] -
          points[simplex.vertices[3]].coords[1];
    cdy = points[simplex.vertices[2]].coords[1] -
          points[simplex.vertices[3]].coords[1];
    adz = points[simplex.vertices[0]].coords[2] -
          points[simplex.vertices[3]].coords[2];
    bdz = points[simplex.vertices[1]].coords[2] -
          points[simplex.vertices[3]].coords[2];
    cdz = points[simplex.vertices[2]].coords[2] -
          points[simplex.vertices[3]].coords[2];

    return adx * (bdy * cdz - bdz * cdy) + bdx * (cdy * adz - cdz * ady) +
           cdx * (ady * bdz - adz * bdy);
  }

  static bool inSphere(const dSimplex<3, Precision> &simplex,
                       const dPoint<3, Precision> &p,
                       const dPoints<3, Precision> &points) {

    Precision aex, bex, cex, dex;
    Precision aey, bey, cey, dey;
    Precision aez, bez, cez, dez;
    Precision alift, blift, clift, dlift;
    Precision ab, bc, cd, da, ac, bd;
    Precision abc, bcd, cda, dab;
    Precision det;

    aex = points[simplex.vertices[0]].coords[0] - p.coords[0];
    bex = points[simplex.vertices[1]].coords[0] - p.coords[0];
    cex = points[simplex.vertices[2]].coords[0] - p.coords[0];
    dex = points[simplex.vertices[3]].coords[0] - p.coords[0];
    aey = points[simplex.vertices[0]].coords[1] - p.coords[1];
    bey = points[simplex.vertices[1]].coords[1] - p.coords[1];
    cey = points[simplex.vertices[2]].coords[1] - p.coords[1];
    dey = points[simplex.vertices[3]].coords[1] - p.coords[1];
    aez = points[simplex.vertices[0]].coords[2] - p.coords[2];
    bez = points[simplex.vertices[1]].coords[2] - p.coords[2];
    cez = points[simplex.vertices[2]].coords[2] - p.coords[2];
    dez = points[simplex.vertices[3]].coords[2] - p.coords[2];

    ab = aex * bey - bex * aey;
    bc = bex * cey - cex * bey;
    cd = cex * dey - dex * cey;
    da = dex * aey - aex * dey;

    ac = aex * cey - cex * aey;
    bd = bex * dey - dex * bey;

    abc = aez * bc - bez * ac + cez * ab;
    bcd = bez * cd - cez * bd + dez * bc;
    cda = cez * da + dez * ac + aez * cd;
    dab = dez * ab + aez * bd + bez * da;

    alift = aex * aex + aey * aey + aez * aez;
    blift = bex * bex + bey * bey + bez * bez;
    clift = cex * cex + cey * cey + cez * cez;
    dlift = dex * dex + dey * dey + dez * dez;

    /*
     * det > 0: e inside  abcd (abcd positive orientation)
     *          e outside abcd (abcd negative orientation)
     * det < 0: e outside abcd (abcd positive orientation)
     *          e inside  abcd (abcd negative orientation)
     *
     * Return a positive value if the point pe lies inside the
     * sphere passing through pa, pb, pc, and pd; a negative value
     * if it lies outside; and zero if the five points are
     * cospherical.  The points pa, pb, pc, and pd must be ordered
     * so that they have a positive orientation (as defined by
     * orient3d()), or the sign of the result will be reversed.
     */

    det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);

    return orientation(simplex, points) * det >= 0;
  }

  static dSphere<3, Precision>
  circumsphere(const dSimplex<3, Precision> &simplex,
               const dPoints<3, Precision> &points) {

    dSphere<3, Precision> sphere;

    Precision den = 2 * (points[simplex.vertices[0]].coords[0] *
                             (points[simplex.vertices[2]].coords[1] *
                                  points[simplex.vertices[3]].coords[2] +
                              points[simplex.vertices[1]].coords[1] *
                                  (points[simplex.vertices[2]].coords[2] -
                                   points[simplex.vertices[3]].coords[2]) -
                              points[simplex.vertices[3]].coords[1] *
                                  points[simplex.vertices[2]].coords[2] -
                              (points[simplex.vertices[2]].coords[1] -
                               points[simplex.vertices[3]].coords[1]) *
                                  points[simplex.vertices[1]].coords[2]) -
                         points[simplex.vertices[1]].coords[0] *
                             (points[simplex.vertices[2]].coords[1] *
                                  points[simplex.vertices[3]].coords[2] -
                              points[simplex.vertices[3]].coords[1] *
                                  points[simplex.vertices[2]].coords[2]) -
                         points[simplex.vertices[0]].coords[1] *
                             (points[simplex.vertices[2]].coords[0] *
                                  points[simplex.vertices[3]].coords[2] +
                              points[simplex.vertices[1]].coords[0] *
                                  (points[simplex.vertices[2]].coords[2] -
                                   points[simplex.vertices[3]].coords[2]) -
                              points[simplex.vertices[3]].coords[0] *
                                  points[simplex.vertices[2]].coords[2] -
                              (points[simplex.vertices[2]].coords[0] -
                               points[simplex.vertices[3]].coords[0]) *
                                  points[simplex.vertices[1]].coords[2]) +
                         points[simplex.vertices[1]].coords[1] *
                             (points[simplex.vertices[2]].coords[0] *
                                  points[simplex.vertices[3]].coords[2] -
                              points[simplex.vertices[3]].coords[0] *
                                  points[simplex.vertices[2]].coords[2]) +
                         points[simplex.vertices[0]].coords[2] *
                             (points[simplex.vertices[2]].coords[0] *
                                  points[simplex.vertices[3]].coords[1] +
                              points[simplex.vertices[1]].coords[0] *
                                  (points[simplex.vertices[2]].coords[1] -
                                   points[simplex.vertices[3]].coords[1]) -
                              points[simplex.vertices[3]].coords[0] *
                                  points[simplex.vertices[2]].coords[1] -
                              (points[simplex.vertices[2]].coords[0] -
                               points[simplex.vertices[3]].coords[0]) *
                                  points[simplex.vertices[1]].coords[1]) -
                         points[simplex.vertices[1]].coords[2] *
                             (points[simplex.vertices[2]].coords[0] *
                                  points[simplex.vertices[3]].coords[1] -
                              points[simplex.vertices[3]].coords[0] *
                                  points[simplex.vertices[2]].coords[1]));

    sphere.center[0] =
        (-points[simplex.vertices[0]].coords[1] *
             (-points[simplex.vertices[2]].coords[2] *
                  (pow(points[simplex.vertices[3]].coords[2], 2) +
                   pow(points[simplex.vertices[3]].coords[1], 2) +
                   pow(points[simplex.vertices[3]].coords[0], 2)) -
              points[simplex.vertices[1]].coords[2] *
                  (-pow(points[simplex.vertices[3]].coords[2], 2) +
                   pow(points[simplex.vertices[2]].coords[2], 2) -
                   pow(points[simplex.vertices[3]].coords[1], 2) +
                   pow(points[simplex.vertices[2]].coords[1], 2) -
                   pow(points[simplex.vertices[3]].coords[0], 2) +
                   pow(points[simplex.vertices[2]].coords[0], 2)) +
              (pow(points[simplex.vertices[2]].coords[2], 2) +
               pow(points[simplex.vertices[2]].coords[1], 2) +
               pow(points[simplex.vertices[2]].coords[0], 2)) *
                  points[simplex.vertices[3]].coords[2] +
              (pow(points[simplex.vertices[1]].coords[2], 2) +
               pow(points[simplex.vertices[1]].coords[1], 2) +
               pow(points[simplex.vertices[1]].coords[0], 2)) *
                  (points[simplex.vertices[2]].coords[2] -
                   points[simplex.vertices[3]].coords[2])) +
         points[simplex.vertices[1]].coords[1] *
             ((pow(points[simplex.vertices[2]].coords[2], 2) +
               pow(points[simplex.vertices[2]].coords[1], 2) +
               pow(points[simplex.vertices[2]].coords[0], 2)) *
                  points[simplex.vertices[3]].coords[2] -
              points[simplex.vertices[2]].coords[2] *
                  (pow(points[simplex.vertices[3]].coords[2], 2) +
                   pow(points[simplex.vertices[3]].coords[1], 2) +
                   pow(points[simplex.vertices[3]].coords[0], 2))) +
         points[simplex.vertices[0]].coords[2] *
             (-points[simplex.vertices[2]].coords[1] *
                  (pow(points[simplex.vertices[3]].coords[2], 2) +
                   pow(points[simplex.vertices[3]].coords[1], 2) +
                   pow(points[simplex.vertices[3]].coords[0], 2)) -
              points[simplex.vertices[1]].coords[1] *
                  (-pow(points[simplex.vertices[3]].coords[2], 2) +
                   pow(points[simplex.vertices[2]].coords[2], 2) -
                   pow(points[simplex.vertices[3]].coords[1], 2) +
                   pow(points[simplex.vertices[2]].coords[1], 2) -
                   pow(points[simplex.vertices[3]].coords[0], 2) +
                   pow(points[simplex.vertices[2]].coords[0], 2)) +
              points[simplex.vertices[3]].coords[1] *
                  (pow(points[simplex.vertices[2]].coords[2], 2) +
                   pow(points[simplex.vertices[2]].coords[1], 2) +
                   pow(points[simplex.vertices[2]].coords[0], 2)) +
              (points[simplex.vertices[2]].coords[1] -
               points[simplex.vertices[3]].coords[1]) *
                  (pow(points[simplex.vertices[1]].coords[2], 2) +
                   pow(points[simplex.vertices[1]].coords[1], 2) +
                   pow(points[simplex.vertices[1]].coords[0], 2))) -
         points[simplex.vertices[1]].coords[2] *
             (points[simplex.vertices[3]].coords[1] *
                  (pow(points[simplex.vertices[2]].coords[2], 2) +
                   pow(points[simplex.vertices[2]].coords[1], 2) +
                   pow(points[simplex.vertices[2]].coords[0], 2)) -
              points[simplex.vertices[2]].coords[1] *
                  (pow(points[simplex.vertices[3]].coords[2], 2) +
                   pow(points[simplex.vertices[3]].coords[1], 2) +
                   pow(points[simplex.vertices[3]].coords[0], 2))) +
         (pow(points[simplex.vertices[0]].coords[2], 2) +
          pow(points[simplex.vertices[0]].coords[1], 2) +
          pow(points[simplex.vertices[0]].coords[0], 2)) *
             (points[simplex.vertices[2]].coords[1] *
                  points[simplex.vertices[3]].coords[2] +
              points[simplex.vertices[1]].coords[1] *
                  (points[simplex.vertices[2]].coords[2] -
                   points[simplex.vertices[3]].coords[2]) -
              points[simplex.vertices[3]].coords[1] *
                  points[simplex.vertices[2]].coords[2] -
              (points[simplex.vertices[2]].coords[1] -
               points[simplex.vertices[3]].coords[1]) *
                  points[simplex.vertices[1]].coords[2]) -
         (pow(points[simplex.vertices[1]].coords[2], 2) +
          pow(points[simplex.vertices[1]].coords[1], 2) +
          pow(points[simplex.vertices[1]].coords[0], 2)) *
             (points[simplex.vertices[2]].coords[1] *
                  points[simplex.vertices[3]].coords[2] -
              points[simplex.vertices[3]].coords[1] *
                  points[simplex.vertices[2]].coords[2])) /
        den;

    sphere.center[1] =
        (points[simplex.vertices[0]].coords[0] *
             (-points[simplex.vertices[2]].coords[2] *
                  (pow(points[simplex.vertices[3]].coords[2], 2) +
                   pow(points[simplex.vertices[3]].coords[1], 2) +
                   pow(points[simplex.vertices[3]].coords[0], 2)) -
              points[simplex.vertices[1]].coords[2] *
                  (-pow(points[simplex.vertices[3]].coords[2], 2) +
                   pow(points[simplex.vertices[2]].coords[2], 2) -
                   pow(points[simplex.vertices[3]].coords[1], 2) +
                   pow(points[simplex.vertices[2]].coords[1], 2) -
                   pow(points[simplex.vertices[3]].coords[0], 2) +
                   pow(points[simplex.vertices[2]].coords[0], 2)) +
              (pow(points[simplex.vertices[2]].coords[2], 2) +
               pow(points[simplex.vertices[2]].coords[1], 2) +
               pow(points[simplex.vertices[2]].coords[0], 2)) *
                  points[simplex.vertices[3]].coords[2] +
              (pow(points[simplex.vertices[1]].coords[2], 2) +
               pow(points[simplex.vertices[1]].coords[1], 2) +
               pow(points[simplex.vertices[1]].coords[0], 2)) *
                  (points[simplex.vertices[2]].coords[2] -
                   points[simplex.vertices[3]].coords[2])) -
         points[simplex.vertices[1]].coords[0] *
             ((pow(points[simplex.vertices[2]].coords[2], 2) +
               pow(points[simplex.vertices[2]].coords[1], 2) +
               pow(points[simplex.vertices[2]].coords[0], 2)) *
                  points[simplex.vertices[3]].coords[2] -
              points[simplex.vertices[2]].coords[2] *
                  (pow(points[simplex.vertices[3]].coords[2], 2) +
                   pow(points[simplex.vertices[3]].coords[1], 2) +
                   pow(points[simplex.vertices[3]].coords[0], 2))) -
         points[simplex.vertices[0]].coords[2] *
             (-points[simplex.vertices[2]].coords[0] *
                  (pow(points[simplex.vertices[3]].coords[2], 2) +
                   pow(points[simplex.vertices[3]].coords[1], 2) +
                   pow(points[simplex.vertices[3]].coords[0], 2)) -
              points[simplex.vertices[1]].coords[0] *
                  (-pow(points[simplex.vertices[3]].coords[2], 2) +
                   pow(points[simplex.vertices[2]].coords[2], 2) -
                   pow(points[simplex.vertices[3]].coords[1], 2) +
                   pow(points[simplex.vertices[2]].coords[1], 2) -
                   pow(points[simplex.vertices[3]].coords[0], 2) +
                   pow(points[simplex.vertices[2]].coords[0], 2)) +
              points[simplex.vertices[3]].coords[0] *
                  (pow(points[simplex.vertices[2]].coords[2], 2) +
                   pow(points[simplex.vertices[2]].coords[1], 2) +
                   pow(points[simplex.vertices[2]].coords[0], 2)) +
              (points[simplex.vertices[2]].coords[0] -
               points[simplex.vertices[3]].coords[0]) *
                  (pow(points[simplex.vertices[1]].coords[2], 2) +
                   pow(points[simplex.vertices[1]].coords[1], 2) +
                   pow(points[simplex.vertices[1]].coords[0], 2))) +
         points[simplex.vertices[1]].coords[2] *
             (points[simplex.vertices[3]].coords[0] *
                  (pow(points[simplex.vertices[2]].coords[2], 2) +
                   pow(points[simplex.vertices[2]].coords[1], 2) +
                   pow(points[simplex.vertices[2]].coords[0], 2)) -
              points[simplex.vertices[2]].coords[0] *
                  (pow(points[simplex.vertices[3]].coords[2], 2) +
                   pow(points[simplex.vertices[3]].coords[1], 2) +
                   pow(points[simplex.vertices[3]].coords[0], 2))) -
         (pow(points[simplex.vertices[0]].coords[2], 2) +
          pow(points[simplex.vertices[0]].coords[1], 2) +
          pow(points[simplex.vertices[0]].coords[0], 2)) *
             (points[simplex.vertices[2]].coords[0] *
                  points[simplex.vertices[3]].coords[2] +
              points[simplex.vertices[1]].coords[0] *
                  (points[simplex.vertices[2]].coords[2] -
                   points[simplex.vertices[3]].coords[2]) -
              points[simplex.vertices[3]].coords[0] *
                  points[simplex.vertices[2]].coords[2] -
              (points[simplex.vertices[2]].coords[0] -
               points[simplex.vertices[3]].coords[0]) *
                  points[simplex.vertices[1]].coords[2]) +
         (pow(points[simplex.vertices[1]].coords[2], 2) +
          pow(points[simplex.vertices[1]].coords[1], 2) +
          pow(points[simplex.vertices[1]].coords[0], 2)) *
             (points[simplex.vertices[2]].coords[0] *
                  points[simplex.vertices[3]].coords[2] -
              points[simplex.vertices[3]].coords[0] *
                  points[simplex.vertices[2]].coords[2])) /
        den;

    sphere.center[2] =
        (-points[simplex.vertices[0]].coords[0] *
             (-points[simplex.vertices[2]].coords[1] *
                  (pow(points[simplex.vertices[3]].coords[2], 2) +
                   pow(points[simplex.vertices[3]].coords[1], 2) +
                   pow(points[simplex.vertices[3]].coords[0], 2)) -
              points[simplex.vertices[1]].coords[1] *
                  (-pow(points[simplex.vertices[3]].coords[2], 2) +
                   pow(points[simplex.vertices[2]].coords[2], 2) -
                   pow(points[simplex.vertices[3]].coords[1], 2) +
                   pow(points[simplex.vertices[2]].coords[1], 2) -
                   pow(points[simplex.vertices[3]].coords[0], 2) +
                   pow(points[simplex.vertices[2]].coords[0], 2)) +
              points[simplex.vertices[3]].coords[1] *
                  (pow(points[simplex.vertices[2]].coords[2], 2) +
                   pow(points[simplex.vertices[2]].coords[1], 2) +
                   pow(points[simplex.vertices[2]].coords[0], 2)) +
              (points[simplex.vertices[2]].coords[1] -
               points[simplex.vertices[3]].coords[1]) *
                  (pow(points[simplex.vertices[1]].coords[2], 2) +
                   pow(points[simplex.vertices[1]].coords[1], 2) +
                   pow(points[simplex.vertices[1]].coords[0], 2))) +
         points[simplex.vertices[1]].coords[0] *
             (points[simplex.vertices[3]].coords[1] *
                  (pow(points[simplex.vertices[2]].coords[2], 2) +
                   pow(points[simplex.vertices[2]].coords[1], 2) +
                   pow(points[simplex.vertices[2]].coords[0], 2)) -
              points[simplex.vertices[2]].coords[1] *
                  (pow(points[simplex.vertices[3]].coords[2], 2) +
                   pow(points[simplex.vertices[3]].coords[1], 2) +
                   pow(points[simplex.vertices[3]].coords[0], 2))) +
         points[simplex.vertices[0]].coords[1] *
             (-points[simplex.vertices[2]].coords[0] *
                  (pow(points[simplex.vertices[3]].coords[2], 2) +
                   pow(points[simplex.vertices[3]].coords[1], 2) +
                   pow(points[simplex.vertices[3]].coords[0], 2)) -
              points[simplex.vertices[1]].coords[0] *
                  (-pow(points[simplex.vertices[3]].coords[2], 2) +
                   pow(points[simplex.vertices[2]].coords[2], 2) -
                   pow(points[simplex.vertices[3]].coords[1], 2) +
                   pow(points[simplex.vertices[2]].coords[1], 2) -
                   pow(points[simplex.vertices[3]].coords[0], 2) +
                   pow(points[simplex.vertices[2]].coords[0], 2)) +
              points[simplex.vertices[3]].coords[0] *
                  (pow(points[simplex.vertices[2]].coords[2], 2) +
                   pow(points[simplex.vertices[2]].coords[1], 2) +
                   pow(points[simplex.vertices[2]].coords[0], 2)) +
              (points[simplex.vertices[2]].coords[0] -
               points[simplex.vertices[3]].coords[0]) *
                  (pow(points[simplex.vertices[1]].coords[2], 2) +
                   pow(points[simplex.vertices[1]].coords[1], 2) +
                   pow(points[simplex.vertices[1]].coords[0], 2))) -
         points[simplex.vertices[1]].coords[1] *
             (points[simplex.vertices[3]].coords[0] *
                  (pow(points[simplex.vertices[2]].coords[2], 2) +
                   pow(points[simplex.vertices[2]].coords[1], 2) +
                   pow(points[simplex.vertices[2]].coords[0], 2)) -
              points[simplex.vertices[2]].coords[0] *
                  (pow(points[simplex.vertices[3]].coords[2], 2) +
                   pow(points[simplex.vertices[3]].coords[1], 2) +
                   pow(points[simplex.vertices[3]].coords[0], 2))) -
         (points[simplex.vertices[2]].coords[0] *
              points[simplex.vertices[3]].coords[1] -
          points[simplex.vertices[3]].coords[0] *
              points[simplex.vertices[2]].coords[1]) *
             (pow(points[simplex.vertices[1]].coords[2], 2) +
              pow(points[simplex.vertices[1]].coords[1], 2) +
              pow(points[simplex.vertices[1]].coords[0], 2)) +
         (points[simplex.vertices[2]].coords[0] *
              points[simplex.vertices[3]].coords[1] +
          points[simplex.vertices[1]].coords[0] *
              (points[simplex.vertices[2]].coords[1] -
               points[simplex.vertices[3]].coords[1]) -
          points[simplex.vertices[3]].coords[0] *
              points[simplex.vertices[2]].coords[1] -
          (points[simplex.vertices[2]].coords[0] -
           points[simplex.vertices[3]].coords[0]) *
              points[simplex.vertices[1]].coords[1]) *
             (pow(points[simplex.vertices[0]].coords[2], 2) +
              pow(points[simplex.vertices[0]].coords[1], 2) +
              pow(points[simplex.vertices[0]].coords[0], 2))) /
        den;

    sphere.radius =
        sqrt(pow(sphere.center[0] - points[simplex.vertices[0]].coords[0], 2) +
             pow(sphere.center[1] - points[simplex.vertices[0]].coords[1], 2) +
             pow(sphere.center[2] - points[simplex.vertices[0]].coords[2], 2));

    return sphere;
  }
};

template <uint D, typename Precision>
Precision
dSimplex<D, Precision>::orientation(const dPoints<D, Precision> &points) const {
  return GeometryHelper<D, Precision>::orientation(*this, points);
}

template <uint D, typename Precision>
bool dSimplex<D, Precision>::inSphere(
    const dPoint<D, Precision> &p, const dPoints<D, Precision> &points) const {
  return GeometryHelper<D, Precision>::inSphere(*this, p, points);
}

template <uint D, typename Precision>
dSphere<D, Precision> dSimplex<D, Precision>::circumsphere(
    const dPoints<D, Precision> &points) const {
  return GeometryHelper<D, Precision>::circumsphere(*this, points);
}

template <uint D, typename Precision>
uint dSimplices<D, Precision>::countDuplicates() const {

  std::atomic<uint> duplicates(0);

  tbb::parallel_for(std::size_t(0), this->bucket_count(), [&](const uint i) {

    for (auto it = this->begin(i); it != this->end(i); ++it) {

      const dSimplex<D, Precision> &s = *it;
      auto simplices = findSimplices(s.vertices, true);
      duplicates += simplices.size() - 1;
    }
  });

  PLOG("Found " << duplicates << " duplicates" << std::endl);

  return duplicates;
}

template <uint D, typename Precision>
CrossCheckReport<D, Precision> dSimplices<D, Precision>::crossCheck(
    const dSimplices<D, Precision> &realDT) const {
  CrossCheckReport<D, Precision> result;
  result.valid = true;

  // check whether all simplices of real DT are present
  for (const auto &realSimplex : realDT) {
    // find my simplex, compares vertices ids
    auto mySimplex = std::find_if(this->begin(), this->end(),
                                  [&](const dSimplex<D, Precision> &s) {
                                    return s.equalVertices(realSimplex);
                                  });

    if (mySimplex == this->end()) {
      LOG("did not find simplex " << realSimplex << std::endl);
      result.valid = false;
      result.missing.insert(realSimplex);
      continue;
    }

    // check neighbors
    for (const auto &n : realSimplex.neighbors) {
      if (!dSimplex<D, Precision>::isFinite(n))
        continue;

      bool found = false;
      for (const auto &nn : mySimplex->neighbors) {
        if (dSimplex<D, Precision>::isFinite(nn) &&
            this->operator[](nn).equalVertices(realDT[n])) {
          found = true;
          break;
        }
      }

      if (!found) {
        LOG("did not find neighbor " << realDT[n] << " of simplex "
                                     << realSimplex << std::endl);
        result.valid = false;
      }
    }
  }

  // check for own simplices that are not in real DT
  for (const auto &mySimplex : *this) {
    auto realSimplex = std::find_if(realDT.begin(), realDT.end(),
                                    [&](const dSimplex<D, Precision> &s) {
                                      return s.equalVertices(mySimplex);
                                    });

    if (realSimplex == realDT.end() /*&& mySimplex.isFinite()*/) {
      LOG("simplex " << mySimplex << " does not exist in real triangulation"
                     << std::endl);
      result.valid = false;
      result.invalid.insert(mySimplex);
    }
  }

  // check whether sizes are equal
  if (dSimplices<D, Precision>::size() != realDT.size()) {
    LOG("my size: " + std::to_string(dSimplices<D, Precision>::size()) +
            " other size: "
        << std::to_string(realDT.size()) + "\n");
    result.valid = false;
  }

  LOG("cross check " << (result.valid ? "" : "NOT ") << "successful"
                     << std::endl);

  return result;
}
template <uint D, typename Precision>
VerificationReport<D, Precision>
dSimplices<D, Precision>::verify(const dPoints<D, Precision> &points) const {
  INDENT
  VerificationReport<D, Precision> result;
  result.valid = true;

  tbb::spin_mutex mtx;
  // verify that every input point is used
  LOG("Checking points" << std::endl);
  Ids usedPoints;
  for (const auto &s : *this) {
    usedPoints.insert(s.vertices.begin(), s.vertices.end());
  }
  if (points != usedPoints) {
    // not all points of input used
    std::stringstream sNotUsed;
    for (const auto &p : points) {
      if (usedPoints.count(p.id) != 1 && p.isFinite()) {
        sNotUsed << p << " ";
        result.valid = false;
      }
    }
    if (!result.valid) {
      LOG("Points of input not used: " << sNotUsed << std::endl);
    }

    std::stringstream sInvalidP;
    for (const auto &p : usedPoints) {
      if (!points.contains(p) && dPoint<D, Precision>::isFinite(p)) {
        sInvalidP << p << " ";
        result.valid = false;
      }
    }
    if (!result.valid) {
      LOG("Used points not in input: " << sInvalidP << std::endl);
    }
  }

  // verify where-used data structure
  LOG("Checking where-used relation" << std::endl);
  for (const auto &s : *this) {
    for (const auto &p : s.vertices) {
      // point p of s not correctly flagged as used in s
      if (std::find(points[p].simplices.begin(), points[p].simplices.end(),
                    s.id) == points[p].simplices.end()) {
        LOG("Point " << p << " NOT flagged as used in " << s << std::endl);
        LOGGER.logContainer(points[p].simplices, Logger::Verbosity::NORMAL,
                            "p.simplices:");
        result.valid = false;
      }
    }
  }
  for (const auto &p : points) {
    for (const auto &id : p.simplices) {
      if (!dSimplex<D, Precision>::isFinite(id) || !this->contains(id))
        continue; // simplex of another triangulation

      // p is flagged as being used in s, but its not
      const auto &s = this->operator[](id);
      if (!std::binary_search(s.vertices.begin(), s.vertices.end(), p.id)) {
        LOG("Point " << p << " SHOULD be used in " << s << std::endl);
        LOGGER.logContainer(p.simplices, Logger::Verbosity::NORMAL,
                            "p.simplices:");
        result.valid = false;
      }
    }
  }

  // verify that all simplices with a shared D-1 simplex are neighbors
  LOG("Checking neighbors" << std::endl);
  tbb::parallel_for(std::size_t(0), this->bucket_count(), [&](const uint i) {

    for (auto it = this->begin(i); it != this->end(i); ++it) {

      const dSimplex<D, Precision> &a = *it;
      for (const auto &b : *this) {
        // a and b are neighbors: the neighbor property is symmetric and the
        // corresponding simplices must be present in the neighbors arrays
        // accordingly
        if ((a.isNeighbor(b) &&
             (!b.isNeighbor(a) || a.neighbors.count(b.id) != 1 ||
              b.neighbors.count(a.id) != 1))
            // a and b are NOT neighbors: must be symmetric and simplices NOT be
            // present in neighbors arrays
            ||
            (!a.isNeighbor(b) &&
             (b.isNeighbor(a) || a.neighbors.count(b.id) != 0 ||
              b.neighbors.count(a.id) != 0))) {
          LOG("Wrong neighbor relation between " << a << " and " << b
                                                 << std::endl);
          tbb::spin_mutex::scoped_lock lock(mtx);
          result.wrongNeighbors.emplace_back(a, b);
          result.valid = false;
        }
      }
    }
  });

  // check in circle criterion
  LOG("Checking empty circle criterion" << std::endl);
  tbb::parallel_for(std::size_t(0), this->bucket_count(), [&](const uint i) {

    for (auto it = this->begin(i); it != this->end(i); ++it) {

      const dSimplex<D, Precision> &s = *it;
      for (const auto &p : points) {
        if (!p.isFinite())
          continue; // skip infinite points

        bool contains = s.contains(p);
        bool inCircle = s.inSphere(p, points);
        if (contains != inCircle) {
          LOG("Point " << p << " is " << (inCircle ? "" : "NOT ")
                       << "in circle of " << s << " but should "
                       << (contains ? "" : "NOT ") << "be" << std::endl);

          tbb::spin_mutex::scoped_lock lock(mtx);
          result.valid = false;

          result.inCircle[s].insert(p.id);
        }
      }
    }
  });
  DEDENT

  LOG("Triangulation is " << (result.valid ? "" : "NOT ") << "valid"
                          << std::endl);

  return result;
}

// specialiations

// float

template class dPoint<2, float>;
template class dPoint<3, float>;

template class dSimplex<2, float>;
template class dSimplex<3, float>;

template class dSimplices<2, float>;
template class dSimplices<3, float>;

// double

template class dPoint<2, double>;
template class dPoint<3, double>;

template class dSimplex<2, double>;
template class dSimplex<3, double>;

template class dSimplices<2, double>;
template class dSimplices<3, double>;
