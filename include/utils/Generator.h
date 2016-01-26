#pragma once

#include "Geometry.h"
#include "Random.h"
#include "SpatialSort.h"

template<uint D, typename Precision>
class PointGenerator {

public:
    dPoints<D, Precision> generate(const tIdType n, const dBox<D, Precision> &bounds,
                                     tGenerator &gen) const {

        auto points = _init(n);
        _gen(points, n, bounds, gen);
        _sort(points);

        return points;

    }

protected:

    dPoints<D, Precision> _init(const tIdType n) const {

        dPoints<D, Precision> points;
        points.reserveUpToIdx(n + 1);

//        for (uint d = 0; d < D; ++d) {
//            points[0].coords[d] = 0;
//        }

        return points;
    }

    void _sort(dPoints<D, Precision> &points) const {
        CGALSpatialSorter<D, Precision> sort;
        sort.sort(points);
    }

    virtual void _gen(dPoints<D, Precision> &points,
                      const tIdType n,
                      const dBox<D, Precision> &bounds,
                      tGenerator &gen) const = 0;

};

template<uint D, typename Precision>
class UniformPointGenerator : public PointGenerator<D, Precision> {

protected:

    void _gen(dPoints<D, Precision> &points,
              const tIdType n,
              const dBox<D, Precision> &bounds,
              tGenerator &gen) const {

        std::uniform_real_distribution<Precision> distribution(0, 1);
        auto dice = std::bind(distribution, gen);

        dVector<D, Precision> dim = bounds.high;
        for (uint d = 0; d < D; ++d) {
            dim[d] -= bounds.low[d];
        }

        for (tIdType i = 1; i <= n; ++i) {
            dPoint<D, Precision> p;

            for (uint d = 0; d < D; ++d) {
                p.coords[d] = bounds.low[d] + dim[d] * dice();
            }


            points[i] = std::move(p);
        }

        // TODO checks for colliding points, straight lines, etc.
    }
};

template<uint D, typename Precision>
class NormalPointGenerator : public PointGenerator<D, Precision> {

protected:

    void _gen(dPoints<D, Precision> &points,
              const tIdType n,
              const dBox<D, Precision> &bounds,
              tGenerator &gen) const {

        std::normal_distribution<Precision> distribution(0, 1);
        auto dice = std::bind(distribution, gen);

        dVector<D, Precision> midpoint = bounds.low;
        dVector<D, Precision> span = bounds.high;

        for (uint d = 0; d < D; ++d) {
            midpoint[d] += (bounds.high[d] - bounds.low[d]) / 2;
            span[d] -= bounds.low[d];
        }

        for (tIdType i = 1; i <= n; ++i) {
            dPoint<D, Precision> p;
            p.coords = midpoint;

            for (uint d = 0; d < D; ++d) {
                p.coords[d] += span[d] * dice();
            }


            points[i] = std::move(p);
        }

        // TODO checks for colliding points, straight lines, etc.
    }
};

template<uint D, typename Precision>
class BubblePointGenerator : public PointGenerator<D, Precision> {

protected:

    void _gen(dPoints<D, Precision> &points,
              const tIdType n,
              const dBox<D, Precision> &bounds,
              tGenerator &gen) const {

        std::normal_distribution<Precision> distribution(0, 1);
        auto dice = std::bind(distribution, gen);

        dVector<D, Precision> gMidpoint = bounds.low;

        for (uint d = 0; d < D; ++d) {
            gMidpoint[d] += (bounds.high[d] - bounds.low[d]) / 2;
        }

        std::vector<dBox<D, Precision>> bubbles;
        bubbles.resize(pow(2, D));

        for (uint i = 0; i < pow(2, D); ++i) {
            for (uint d = 0; d < D; ++d) {
                bubbles[i].low[d] =
                        i & (1 << d) ? gMidpoint[d] : bounds.low[d];
                bubbles[i].high[d] =
                        i & (1 << d) ? bounds.high[d] : gMidpoint[d];
            }
        }

        tIdType ppB = n / bubbles.size(); // points per bubble
        for(uint b = 0; b < bubbles.size(); ++b) {
            const auto & bubble = bubbles[b];

            dVector<D, Precision> midpoint = bubble.low;
            dVector<D, Precision> span = bubble.high;

            for (uint d = 0; d < D; ++d) {
                midpoint[d] += (bubble.high[d] - bubble.low[d]) / 2;
                span[d] -= bubble.low[d];
            }

            for (tIdType i = b*ppB + 1; i <= (b+1)*ppB; ++i) {
                dPoint<D, Precision> p;
                p.coords = midpoint;

                for (uint d = 0; d < D; ++d) {
                    p.coords[d] += span[d] * dice();
                }


                points[i] = std::move(p);
            }
        }

        // TODO checks for colliding points, straight lines, etc.
    }
};

template<uint D, typename Precision>
class EllipsoidPointGenerator : public PointGenerator<D, Precision> {

protected:

    void _gen(dPoints<D, Precision> &points,
              const tIdType n,
              const dBox<D, Precision> &bounds,
              tGenerator &gen) const {

        std::normal_distribution<Precision> distribution(0, 1);
        auto dice = std::bind(distribution, gen);

        dVector<D, Precision> radius;
        dVector<D, Precision> center = bounds.low;

        for (uint d = 0; d < D; ++d) {
            center[d] += (bounds.high[d] - bounds.low[d]) / 2;
            radius[d] = .45 * (bounds.high[d] - bounds.low[d]); //give it some space at the edge
        }

        // we generate normal distributed points and map them to points on a unit circle/sphere
        // these points are uniform distributed on the circle/sphere
        // we then map the circle/sphere points to the ellipse/ellipsoid
        // the resulting points are not uniform distributed, but close enough

        for (tIdType i = 1; i <= n; ++i) {
            dPoint<D, Precision> p;

            //normal distributed point around origin
            Precision r = 0; //radius
            for (uint d = 0; d < D; ++d) {
                p.coords[d] = dice();
                r += p.coords[d] * p.coords[d];
            }

            r = std::sqrt(r);

            // now map them to the ellipse/ellipsoid
            for (uint d = 0; d < D; ++d) {
                //            origin      distortion  point on circle/sphere
                p.coords[d] = center[d] + radius[d] * (p.coords[d] / r);
            }


            points[i] = std::move(p);
        }

        // TODO checks for colliding points, straight lines, etc.
    }
};

template<uint D, typename Precision>
class SkewLinePointGenerator : public PointGenerator<D, Precision> {

protected:

    void _gen(dPoints<D, Precision> &points,
              const tIdType n,
              const dBox<D, Precision> &bounds,
              tGenerator &gen) const {
        _gen(points, n, bounds, gen, .1);
    }

    void _gen(dPoints<D, Precision> &points,
              const tIdType n,
              const dBox<D, Precision> &bounds,
              tGenerator &gen,
              Precision skew) const {

        std::uniform_real_distribution<Precision> distribution(0, 1);
        auto dice = std::bind(distribution, gen);

        dVector<D, Precision> a = bounds.low;
        a[1] += skew * (bounds.high[1] - bounds.low[1]);
        dVector<D, Precision> x = bounds.high;
        for (uint d = 0; d < D; ++d) {
            x[d] -= a[d];
        }

        dVector<D, Precision> b = bounds.low;
        b[0] = bounds.high[0]; //opposite corner
        dVector<D, Precision> y = bounds.high;
        y[0] = bounds.low[0]; //opposite corner
        for (uint d = 0; d < D; ++d) {
            y[d] -= b[d];
        }


        // we have two lines a + rx and b + ry
        // a is skewed
        // b crosses the opposite corner

        for (tIdType i = 1; i < n; i += 2) {
            dPoint<D, Precision> p1, p2;

            Precision r = dice();

            for (uint d = 0; d < D; ++d) {
                p1.coords[d] = a[d] + r * x[d];
                p2.coords[d] = b[d] + r * y[d];
            }


            points[i] = std::move(p1);
            points[i + 1] = std::move(p2);
        }

        // TODO checks for colliding points, straight lines, etc.
    }
};

template<uint D, typename Precision>
class GeneratorFactory {

public:
    static std::unique_ptr<PointGenerator<D, Precision>> make(const unsigned char type) {

        switch (type) {
            case 'e':
                return std::make_unique<EllipsoidPointGenerator<D, Precision>>();
            case 'l':
                return std::make_unique<SkewLinePointGenerator<D, Precision>>();
            case 'n':
                return std::make_unique<NormalPointGenerator<D, Precision>>();
            case 'b':
                return std::make_unique<BubblePointGenerator<D, Precision>>();
            case 'u':
            default:
                return std::make_unique<UniformPointGenerator<D, Precision>>();
        }
    }

};