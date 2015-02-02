#include "Geometry.h"

#include <algorithm>
#include <sstream>

// static variables
template <uint D> constexpr uint dPoint<D>::cINF;

template <uint D> constexpr uint dSimplex<D>::cINF;
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

// specializations

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

template <>
tCoordinate dSimplex<2>::orientation(const dPoints<2> &points) const {
  tCoordinate acx, bcx, acy, bcy;

  acx = points[vertices[0]].coords[0] - points[vertices[2]].coords[0];
  bcx = points[vertices[1]].coords[0] - points[vertices[2]].coords[0];
  acy = points[vertices[0]].coords[1] - points[vertices[2]].coords[1];
  bcy = points[vertices[1]].coords[1] - points[vertices[2]].coords[1];
  return acx * bcy - acy * bcx;
}

template <>
bool dSimplex<2>::inSphere(const dPoint<2> &p, const dPoints<2> &points) const {
  tCoordinate adx, ady, bdx, bdy, cdx, cdy;
  tCoordinate abdet, bcdet, cadet;
  tCoordinate alift, blift, clift;
  tCoordinate det;

  adx = points[vertices[0]].coords[0] - p.coords[0];
  ady = points[vertices[0]].coords[1] - p.coords[1];
  bdx = points[vertices[1]].coords[0] - p.coords[0];
  bdy = points[vertices[1]].coords[1] - p.coords[1];
  cdx = points[vertices[2]].coords[0] - p.coords[0];
  cdy = points[vertices[2]].coords[1] - p.coords[1];

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

  return this->orientation(points) * det >= 0;
}

template <>
dSphere<2> dSimplex<2>::circumsphere(const dPoints<2> &points) const {

  dSphere<2> sphere;

  tCoordinate den =
      2 * (points[vertices[0]].coords[0] *
               (points[vertices[1]].coords[1] - points[vertices[2]].coords[1]) +
           points[vertices[1]].coords[0] *
               (points[vertices[2]].coords[1] - points[vertices[0]].coords[1]) +
           points[vertices[2]].coords[0] *
               (points[vertices[0]].coords[1] - points[vertices[1]].coords[1]));

  sphere.center[0] =
      ((pow(points[vertices[0]].coords[0], 2) +
        pow(points[vertices[0]].coords[1], 2)) *
           (points[vertices[1]].coords[1] - points[vertices[2]].coords[1]) +
       (pow(points[vertices[1]].coords[0], 2) +
        pow(points[vertices[1]].coords[1], 2)) *
           (points[vertices[2]].coords[1] - points[vertices[0]].coords[1]) +
       (pow(points[vertices[2]].coords[0], 2) +
        pow(points[vertices[2]].coords[1], 2)) *
           (points[vertices[0]].coords[1] - points[vertices[1]].coords[1])) /
      den;

  sphere.center[1] =
      ((pow(points[vertices[0]].coords[0], 2) +
        pow(points[vertices[0]].coords[1], 2)) *
           (points[vertices[2]].coords[0] - points[vertices[1]].coords[0]) +
       (pow(points[vertices[1]].coords[0], 2) +
        pow(points[vertices[1]].coords[1], 2)) *
           (points[vertices[0]].coords[0] - points[vertices[2]].coords[0]) +
       (pow(points[vertices[2]].coords[0], 2) +
        pow(points[vertices[2]].coords[1], 2)) *
           (points[vertices[1]].coords[0] - points[vertices[0]].coords[0])) /
      den;

  sphere.radius =
      sqrt(pow(sphere.center[0] - points[vertices[0]].coords[0], 2) +
           pow(sphere.center[1] - points[vertices[0]].coords[1], 2));

  return sphere;
}

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

template <>
tCoordinate dSimplex<3>::orientation(const dPoints<3> &points) const {

  tCoordinate adx, bdx, cdx;
  tCoordinate ady, bdy, cdy;
  tCoordinate adz, bdz, cdz;

  adx = points[vertices[0]].coords[0] - points[vertices[3]].coords[0];
  bdx = points[vertices[1]].coords[0] - points[vertices[3]].coords[0];
  cdx = points[vertices[2]].coords[0] - points[vertices[3]].coords[0];
  ady = points[vertices[0]].coords[1] - points[vertices[3]].coords[1];
  bdy = points[vertices[1]].coords[1] - points[vertices[3]].coords[1];
  cdy = points[vertices[2]].coords[1] - points[vertices[3]].coords[1];
  adz = points[vertices[0]].coords[2] - points[vertices[3]].coords[2];
  bdz = points[vertices[1]].coords[2] - points[vertices[3]].coords[2];
  cdz = points[vertices[2]].coords[2] - points[vertices[3]].coords[2];

  return adx * (bdy * cdz - bdz * cdy) + bdx * (cdy * adz - cdz * ady) +
         cdx * (ady * bdz - adz * bdy);
}

template <>
bool dSimplex<3>::inSphere(const dPoint<3> &p, const dPoints<3> &points) const {

  tCoordinate aex, bex, cex, dex;
  tCoordinate aey, bey, cey, dey;
  tCoordinate aez, bez, cez, dez;
  tCoordinate alift, blift, clift, dlift;
  tCoordinate ab, bc, cd, da, ac, bd;
  tCoordinate abc, bcd, cda, dab;
  tCoordinate det;

  aex = points[vertices[0]].coords[0] - p.coords[0];
  bex = points[vertices[1]].coords[0] - p.coords[0];
  cex = points[vertices[2]].coords[0] - p.coords[0];
  dex = points[vertices[3]].coords[0] - p.coords[0];
  aey = points[vertices[0]].coords[1] - p.coords[1];
  bey = points[vertices[1]].coords[1] - p.coords[1];
  cey = points[vertices[2]].coords[1] - p.coords[1];
  dey = points[vertices[3]].coords[1] - p.coords[1];
  aez = points[vertices[0]].coords[2] - p.coords[2];
  bez = points[vertices[1]].coords[2] - p.coords[2];
  cez = points[vertices[2]].coords[2] - p.coords[2];
  dez = points[vertices[3]].coords[2] - p.coords[2];

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

  return this->orientation(points) * det >= 0;
}

template <>
dSphere<3> dSimplex<3>::circumsphere(const dPoints<3> &points) const {

  dSphere<3> sphere;

  tCoordinate den =
      2 *
      (points[vertices[0]].coords[0] *
           (points[vertices[2]].coords[1] * points[vertices[3]].coords[2] +
            points[vertices[1]].coords[1] * (points[vertices[2]].coords[2] -
                                             points[vertices[3]].coords[2]) -
            points[vertices[3]].coords[1] * points[vertices[2]].coords[2] -
            (points[vertices[2]].coords[1] - points[vertices[3]].coords[1]) *
                points[vertices[1]].coords[2]) -
       points[vertices[1]].coords[0] *
           (points[vertices[2]].coords[1] * points[vertices[3]].coords[2] -
            points[vertices[3]].coords[1] * points[vertices[2]].coords[2]) -
       points[vertices[0]].coords[1] *
           (points[vertices[2]].coords[0] * points[vertices[3]].coords[2] +
            points[vertices[1]].coords[0] * (points[vertices[2]].coords[2] -
                                             points[vertices[3]].coords[2]) -
            points[vertices[3]].coords[0] * points[vertices[2]].coords[2] -
            (points[vertices[2]].coords[0] - points[vertices[3]].coords[0]) *
                points[vertices[1]].coords[2]) +
       points[vertices[1]].coords[1] *
           (points[vertices[2]].coords[0] * points[vertices[3]].coords[2] -
            points[vertices[3]].coords[0] * points[vertices[2]].coords[2]) +
       points[vertices[0]].coords[2] *
           (points[vertices[2]].coords[0] * points[vertices[3]].coords[1] +
            points[vertices[1]].coords[0] * (points[vertices[2]].coords[1] -
                                             points[vertices[3]].coords[1]) -
            points[vertices[3]].coords[0] * points[vertices[2]].coords[1] -
            (points[vertices[2]].coords[0] - points[vertices[3]].coords[0]) *
                points[vertices[1]].coords[1]) -
       points[vertices[1]].coords[2] *
           (points[vertices[2]].coords[0] * points[vertices[3]].coords[1] -
            points[vertices[3]].coords[0] * points[vertices[2]].coords[1]));

  sphere.center[0] =
      (-points[vertices[0]].coords[1] *
           (-points[vertices[2]].coords[2] *
                (pow(points[vertices[3]].coords[2], 2) +
                 pow(points[vertices[3]].coords[1], 2) +
                 pow(points[vertices[3]].coords[0], 2)) -
            points[vertices[1]].coords[2] *
                (-pow(points[vertices[3]].coords[2], 2) +
                 pow(points[vertices[2]].coords[2], 2) -
                 pow(points[vertices[3]].coords[1], 2) +
                 pow(points[vertices[2]].coords[1], 2) -
                 pow(points[vertices[3]].coords[0], 2) +
                 pow(points[vertices[2]].coords[0], 2)) +
            (pow(points[vertices[2]].coords[2], 2) +
             pow(points[vertices[2]].coords[1], 2) +
             pow(points[vertices[2]].coords[0], 2)) *
                points[vertices[3]].coords[2] +
            (pow(points[vertices[1]].coords[2], 2) +
             pow(points[vertices[1]].coords[1], 2) +
             pow(points[vertices[1]].coords[0], 2)) *
                (points[vertices[2]].coords[2] -
                 points[vertices[3]].coords[2])) +
       points[vertices[1]].coords[1] *
           ((pow(points[vertices[2]].coords[2], 2) +
             pow(points[vertices[2]].coords[1], 2) +
             pow(points[vertices[2]].coords[0], 2)) *
                points[vertices[3]].coords[2] -
            points[vertices[2]].coords[2] *
                (pow(points[vertices[3]].coords[2], 2) +
                 pow(points[vertices[3]].coords[1], 2) +
                 pow(points[vertices[3]].coords[0], 2))) +
       points[vertices[0]].coords[2] *
           (-points[vertices[2]].coords[1] *
                (pow(points[vertices[3]].coords[2], 2) +
                 pow(points[vertices[3]].coords[1], 2) +
                 pow(points[vertices[3]].coords[0], 2)) -
            points[vertices[1]].coords[1] *
                (-pow(points[vertices[3]].coords[2], 2) +
                 pow(points[vertices[2]].coords[2], 2) -
                 pow(points[vertices[3]].coords[1], 2) +
                 pow(points[vertices[2]].coords[1], 2) -
                 pow(points[vertices[3]].coords[0], 2) +
                 pow(points[vertices[2]].coords[0], 2)) +
            points[vertices[3]].coords[1] *
                (pow(points[vertices[2]].coords[2], 2) +
                 pow(points[vertices[2]].coords[1], 2) +
                 pow(points[vertices[2]].coords[0], 2)) +
            (points[vertices[2]].coords[1] - points[vertices[3]].coords[1]) *
                (pow(points[vertices[1]].coords[2], 2) +
                 pow(points[vertices[1]].coords[1], 2) +
                 pow(points[vertices[1]].coords[0], 2))) -
       points[vertices[1]].coords[2] *
           (points[vertices[3]].coords[1] *
                (pow(points[vertices[2]].coords[2], 2) +
                 pow(points[vertices[2]].coords[1], 2) +
                 pow(points[vertices[2]].coords[0], 2)) -
            points[vertices[2]].coords[1] *
                (pow(points[vertices[3]].coords[2], 2) +
                 pow(points[vertices[3]].coords[1], 2) +
                 pow(points[vertices[3]].coords[0], 2))) +
       (pow(points[vertices[0]].coords[2], 2) +
        pow(points[vertices[0]].coords[1], 2) +
        pow(points[vertices[0]].coords[0], 2)) *
           (points[vertices[2]].coords[1] * points[vertices[3]].coords[2] +
            points[vertices[1]].coords[1] * (points[vertices[2]].coords[2] -
                                             points[vertices[3]].coords[2]) -
            points[vertices[3]].coords[1] * points[vertices[2]].coords[2] -
            (points[vertices[2]].coords[1] - points[vertices[3]].coords[1]) *
                points[vertices[1]].coords[2]) -
       (pow(points[vertices[1]].coords[2], 2) +
        pow(points[vertices[1]].coords[1], 2) +
        pow(points[vertices[1]].coords[0], 2)) *
           (points[vertices[2]].coords[1] * points[vertices[3]].coords[2] -
            points[vertices[3]].coords[1] * points[vertices[2]].coords[2])) /
      den;

  sphere.center[1] =
      (points[vertices[0]].coords[0] *
           (-points[vertices[2]].coords[2] *
                (pow(points[vertices[3]].coords[2], 2) +
                 pow(points[vertices[3]].coords[1], 2) +
                 pow(points[vertices[3]].coords[0], 2)) -
            points[vertices[1]].coords[2] *
                (-pow(points[vertices[3]].coords[2], 2) +
                 pow(points[vertices[2]].coords[2], 2) -
                 pow(points[vertices[3]].coords[1], 2) +
                 pow(points[vertices[2]].coords[1], 2) -
                 pow(points[vertices[3]].coords[0], 2) +
                 pow(points[vertices[2]].coords[0], 2)) +
            (pow(points[vertices[2]].coords[2], 2) +
             pow(points[vertices[2]].coords[1], 2) +
             pow(points[vertices[2]].coords[0], 2)) *
                points[vertices[3]].coords[2] +
            (pow(points[vertices[1]].coords[2], 2) +
             pow(points[vertices[1]].coords[1], 2) +
             pow(points[vertices[1]].coords[0], 2)) *
                (points[vertices[2]].coords[2] -
                 points[vertices[3]].coords[2])) -
       points[vertices[1]].coords[0] *
           ((pow(points[vertices[2]].coords[2], 2) +
             pow(points[vertices[2]].coords[1], 2) +
             pow(points[vertices[2]].coords[0], 2)) *
                points[vertices[3]].coords[2] -
            points[vertices[2]].coords[2] *
                (pow(points[vertices[3]].coords[2], 2) +
                 pow(points[vertices[3]].coords[1], 2) +
                 pow(points[vertices[3]].coords[0], 2))) -
       points[vertices[0]].coords[2] *
           (-points[vertices[2]].coords[0] *
                (pow(points[vertices[3]].coords[2], 2) +
                 pow(points[vertices[3]].coords[1], 2) +
                 pow(points[vertices[3]].coords[0], 2)) -
            points[vertices[1]].coords[0] *
                (-pow(points[vertices[3]].coords[2], 2) +
                 pow(points[vertices[2]].coords[2], 2) -
                 pow(points[vertices[3]].coords[1], 2) +
                 pow(points[vertices[2]].coords[1], 2) -
                 pow(points[vertices[3]].coords[0], 2) +
                 pow(points[vertices[2]].coords[0], 2)) +
            points[vertices[3]].coords[0] *
                (pow(points[vertices[2]].coords[2], 2) +
                 pow(points[vertices[2]].coords[1], 2) +
                 pow(points[vertices[2]].coords[0], 2)) +
            (points[vertices[2]].coords[0] - points[vertices[3]].coords[0]) *
                (pow(points[vertices[1]].coords[2], 2) +
                 pow(points[vertices[1]].coords[1], 2) +
                 pow(points[vertices[1]].coords[0], 2))) +
       points[vertices[1]].coords[2] *
           (points[vertices[3]].coords[0] *
                (pow(points[vertices[2]].coords[2], 2) +
                 pow(points[vertices[2]].coords[1], 2) +
                 pow(points[vertices[2]].coords[0], 2)) -
            points[vertices[2]].coords[0] *
                (pow(points[vertices[3]].coords[2], 2) +
                 pow(points[vertices[3]].coords[1], 2) +
                 pow(points[vertices[3]].coords[0], 2))) -
       (pow(points[vertices[0]].coords[2], 2) +
        pow(points[vertices[0]].coords[1], 2) +
        pow(points[vertices[0]].coords[0], 2)) *
           (points[vertices[2]].coords[0] * points[vertices[3]].coords[2] +
            points[vertices[1]].coords[0] * (points[vertices[2]].coords[2] -
                                             points[vertices[3]].coords[2]) -
            points[vertices[3]].coords[0] * points[vertices[2]].coords[2] -
            (points[vertices[2]].coords[0] - points[vertices[3]].coords[0]) *
                points[vertices[1]].coords[2]) +
       (pow(points[vertices[1]].coords[2], 2) +
        pow(points[vertices[1]].coords[1], 2) +
        pow(points[vertices[1]].coords[0], 2)) *
           (points[vertices[2]].coords[0] * points[vertices[3]].coords[2] -
            points[vertices[3]].coords[0] * points[vertices[2]].coords[2])) /
      den;

  sphere.center[2] =
      (-points[vertices[0]].coords[0] *
           (-points[vertices[2]].coords[1] *
                (pow(points[vertices[3]].coords[2], 2) +
                 pow(points[vertices[3]].coords[1], 2) +
                 pow(points[vertices[3]].coords[0], 2)) -
            points[vertices[1]].coords[1] *
                (-pow(points[vertices[3]].coords[2], 2) +
                 pow(points[vertices[2]].coords[2], 2) -
                 pow(points[vertices[3]].coords[1], 2) +
                 pow(points[vertices[2]].coords[1], 2) -
                 pow(points[vertices[3]].coords[0], 2) +
                 pow(points[vertices[2]].coords[0], 2)) +
            points[vertices[3]].coords[1] *
                (pow(points[vertices[2]].coords[2], 2) +
                 pow(points[vertices[2]].coords[1], 2) +
                 pow(points[vertices[2]].coords[0], 2)) +
            (points[vertices[2]].coords[1] - points[vertices[3]].coords[1]) *
                (pow(points[vertices[1]].coords[2], 2) +
                 pow(points[vertices[1]].coords[1], 2) +
                 pow(points[vertices[1]].coords[0], 2))) +
       points[vertices[1]].coords[0] *
           (points[vertices[3]].coords[1] *
                (pow(points[vertices[2]].coords[2], 2) +
                 pow(points[vertices[2]].coords[1], 2) +
                 pow(points[vertices[2]].coords[0], 2)) -
            points[vertices[2]].coords[1] *
                (pow(points[vertices[3]].coords[2], 2) +
                 pow(points[vertices[3]].coords[1], 2) +
                 pow(points[vertices[3]].coords[0], 2))) +
       points[vertices[0]].coords[1] *
           (-points[vertices[2]].coords[0] *
                (pow(points[vertices[3]].coords[2], 2) +
                 pow(points[vertices[3]].coords[1], 2) +
                 pow(points[vertices[3]].coords[0], 2)) -
            points[vertices[1]].coords[0] *
                (-pow(points[vertices[3]].coords[2], 2) +
                 pow(points[vertices[2]].coords[2], 2) -
                 pow(points[vertices[3]].coords[1], 2) +
                 pow(points[vertices[2]].coords[1], 2) -
                 pow(points[vertices[3]].coords[0], 2) +
                 pow(points[vertices[2]].coords[0], 2)) +
            points[vertices[3]].coords[0] *
                (pow(points[vertices[2]].coords[2], 2) +
                 pow(points[vertices[2]].coords[1], 2) +
                 pow(points[vertices[2]].coords[0], 2)) +
            (points[vertices[2]].coords[0] - points[vertices[3]].coords[0]) *
                (pow(points[vertices[1]].coords[2], 2) +
                 pow(points[vertices[1]].coords[1], 2) +
                 pow(points[vertices[1]].coords[0], 2))) -
       points[vertices[1]].coords[1] *
           (points[vertices[3]].coords[0] *
                (pow(points[vertices[2]].coords[2], 2) +
                 pow(points[vertices[2]].coords[1], 2) +
                 pow(points[vertices[2]].coords[0], 2)) -
            points[vertices[2]].coords[0] *
                (pow(points[vertices[3]].coords[2], 2) +
                 pow(points[vertices[3]].coords[1], 2) +
                 pow(points[vertices[3]].coords[0], 2))) -
       (points[vertices[2]].coords[0] * points[vertices[3]].coords[1] -
        points[vertices[3]].coords[0] * points[vertices[2]].coords[1]) *
           (pow(points[vertices[1]].coords[2], 2) +
            pow(points[vertices[1]].coords[1], 2) +
            pow(points[vertices[1]].coords[0], 2)) +
       (points[vertices[2]].coords[0] * points[vertices[3]].coords[1] +
        points[vertices[1]].coords[0] *
            (points[vertices[2]].coords[1] - points[vertices[3]].coords[1]) -
        points[vertices[3]].coords[0] * points[vertices[2]].coords[1] -
        (points[vertices[2]].coords[0] - points[vertices[3]].coords[0]) *
            points[vertices[1]].coords[1]) *
           (pow(points[vertices[0]].coords[2], 2) +
            pow(points[vertices[0]].coords[1], 2) +
            pow(points[vertices[0]].coords[0], 2))) /
      den;

  sphere.radius =
      sqrt(pow(sphere.center[0] - points[vertices[0]].coords[0], 2) +
           pow(sphere.center[1] - points[vertices[0]].coords[1], 2) +
           pow(sphere.center[2] - points[vertices[0]].coords[2], 2));

  return sphere;
}

template <uint D>
CrossCheckReport<D>
dSimplices<D>::crossCheck(const dSimplices<D> &realDT) const {
  CrossCheckReport<D> result;
  result.valid = true;

  // check whether all simplices of real DT are present
  for (const auto &realSimplex : realDT) {
    // find my simplex, compares simplex id or vertices ids
    auto mySimplex = std::find_if(
        this->begin(), this->end(),
        [&](const dSimplex<D> &s) { return s.equalVertices(realSimplex); });

    if (mySimplex == this->end()) {
      LOG << "did not find simplex " << realSimplex << std::endl;
      result.valid = false;
      result.missing.insert(realSimplex);
      continue;
    }

    // check neighbors
    for (const auto &n : realSimplex.neighbors) {
      if (!dSimplex<D>::isFinite(n))
        continue;

      bool found = false;
      for (const auto &nn : mySimplex->neighbors) {
        if (dSimplex<D>::isFinite(nn) &&
            this->operator[](nn).equalVertices(realDT[n])) {
          found = true;
          break;
        }
      }

      if (!found) {
        LOG << "did not find neighbor " << realDT[n] << " of simplex "
            << realSimplex << std::endl;
        result.valid = false;
      }
    }
  }

  // check for own simplices that are not in real DT
  uint infSimplices = 0;
  for (const auto &mySimplex : *this) {
    auto realSimplex = std::find_if(
        realDT.begin(), realDT.end(),
        [&](const dSimplex<D> &s) { return s.equalVertices(mySimplex); });

    infSimplices += !mySimplex.isFinite();

    if (realSimplex == realDT.end() && mySimplex.isFinite()) {
      LOG << "simplex " << mySimplex << " does not exist in real triangulation"
          << std::endl;
      result.valid = false;
      result.invalid.insert(mySimplex);
    }
  }

  // check whether sizes are equal, account for infinite simplices
  if (dSimplices<D>::size() - infSimplices != realDT.size()) {
    LOG << "my size: " << dSimplices<D>::size() << " - " << infSimplices
        << " other size: " << realDT.size() << std::endl;
    result.valid = false;
  }

  LOG << "cross check " << (result.valid ? "" : "NOT ") << "successful"
      << std::endl;

  return result;
}
template <uint D>
VerificationReport<D> dSimplices<D>::verify(const dPoints<D> &points) const {
  INDENT
  VerificationReport<D> result;
  result.valid = true;

  // verify that every input point is used
  LOG << "Checking points" << std::endl;
  std::set<uint> usedPoints;
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
      LOG << "Points of input not used: " << sNotUsed << std::endl;
    }

    std::stringstream sInvalidP;
    for (const auto &p : usedPoints) {
      if (!points.contains(p) && dPoint<D>::isFinite(p)) {
        sInvalidP << p << " ";
        result.valid = false;
      }
    }
    if (!result.valid) {
      LOG << "Used points not in input: " << sInvalidP << std::endl;
    }
  }

  // verify where-used data structure
  LOG << "Checking where-used relation" << std::endl;
  for (const auto &s : *this) {
    for (const auto &p : s.vertices) {
      // point p of s not correctly flagged as used in s
      if (points[p].simplices.count(s.id) != 1) {
        LOG << "Point " << p << " NOT flagged as used in " << s << std::endl;
        ;
        Logger::getInstance().logContainer(
            points[p].simplices, Logger::Verbosity::NORMAL, "p.simplices");
        result.valid = false;
      }
    }
  }
  for (const auto &p : points) {
    for (const auto &id : p.simplices) {
      if (!this->contains(id))
        continue; // simplex of another triangulation

      // p is flagged as being used in s, but its not
      const auto &s = this->operator[](id);
      if (std::find(s.vertices.begin(), s.vertices.end(), p.id) ==
          s.vertices.end()) {
        LOG << "Point " << p << " SHOULD be used in " << s << std::endl;
        Logger::getInstance().logContainer(
            p.simplices, Logger::Verbosity::NORMAL, "p.simplices");
        result.valid = false;
      }
    }
  }

  // verify that all simplices with a shared D-1 simplex are neighbors
  LOG << "Checking neighbors" << std::endl;
  for (const auto &a : *this) {
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
        LOG << "Wrong neighbor relation between " << a << " and " << b
            << std::endl;
        result.valid = false;
      }
    }
  }

  // check in circle criterion
  LOG << "Checking empty circle criterion" << std::endl;
  for (const auto &s : *this) {
    for (const auto &p : points) {
      if (!p.isFinite())
        continue; // skip infinite points

      bool contains = s.contains(p);
      bool inCircle = s.inSphere(p, points);
      if (contains != inCircle) {
        LOG << "Point " << p << " is " << (inCircle ? "" : "NOT ")
            << "in circle of " << s << " but should "
            << (contains ? "" : "NOT ") << "be" << std::endl;
        result.valid = false;

        result.inCircle[s].insert(p.id);
      }
    }
  }
  DEDENT

  LOG << "Triangulation is " << (result.valid ? "" : "NOT ") << "valid"
      << std::endl;

  return result;
}

// specialiations

template class dPoint<2>;
template class dPoint<3>;

template class dSimplex<2>;
template class dSimplex<3>;

template class dSimplices<2>;
template class dSimplices<3>;
