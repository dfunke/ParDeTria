#pragma once

#include <fstream>

// boost
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>

template<class Object>
void storeObject(Object &o, const std::string &file) {

    std::ofstream f(file, std::ios::trunc | std::ios::out | std::ios::binary);
    boost::archive::binary_oarchive oa(f);
    oa << o;
}

template<class Object>
Object loadObject(const std::string &file) {

    std::ifstream f(file, std::ios::in | std::ios::binary);
    boost::archive::binary_iarchive ia(f);

    Object o;
    ia >> o;

    return o;
}

namespace boost {
    namespace serialization {

        template<class Archive, uint D, typename Precision>
        void serialize(Archive &ar, dVector<D, Precision> &v,
                       __attribute((unused)) const unsigned int version) {
            ar & v;
        }

        template<class Archive, uint D, typename Precision>
        void serialize(Archive &ar, dBox<D, Precision> &b,
                       __attribute((unused)) const unsigned int version) {
            ar & b.low;
            ar & b.high;
        }

        template<class Archive, uint D, typename Precision>
        void serialize(Archive &ar, dSphere<D, Precision> &s,
                       __attribute((unused)) const unsigned int version) {
            ar & s.center;
            ar & s.radius;
        }

//        template<class Archive, typename V>
//        void serialize(Archive &ar, IndexedVector<V> &v,
//                       __attribute((unused)) const unsigned int version) {
//            ar & boost::serialization::base_object<std::vector<V>>(v);
//        }

// points

        template<class Archive, uint D, typename Precision>
        void serialize(Archive &ar, dPoint<D, Precision> &p,
                       __attribute((unused)) const unsigned int version) {
            //ar & p.id;
            ar & p.coords;
            // ar &p.simplices;
        }

        template<class Archive, uint D, typename Precision>
        void serialize(Archive &ar, dPointStats<D, Precision> &s,
                       __attribute((unused)) const unsigned int version) {
            ar & s.min;
            ar & s.mid;
            ar & s.max;
        }

        template<class Archive, uint D, typename Precision>
        void serialize(Archive &ar, dPoints<D, Precision> &ps,
                       __attribute((unused)) const unsigned int version) {
            ar & boost::serialization::base_object<VectorAdapter2<dPoint<D, Precision>>>(
                    ps);
        }

// simplices

        template<class Archive, uint D, typename Precision>
        void serialize(Archive &ar, dSimplex<D, Precision> &s,
                       __attribute((unused)) const unsigned int version) {
            ar & s.id;
            ar & s.vertices;
            ar & s.neighbors;
            ar & s.vFingerprint;
            ar & s.mark;
        }

//        template<class Archive, uint D, typename Precision>
//        void serialize(Archive &ar, dSimplices<D, Precision> &ss,
//                       __attribute((unused)) const unsigned int version) {
//            ar & boost::serialization::base_object<Concurrent_BlockedArray<dSimplex<D, Precision>>>(
//                    ss);
//        }

// partitions

        template<class Archive, uint D, typename Precision>
        void serialize(Archive &ar, Partition<D, Precision> &p,
                       __attribute((unused)) const unsigned int version) {
            ar & p.id;
            ar & p.bounds;
            ar & p.points;
        }

//        template<class Archive, uint D, typename Precision>
//        void serialize(Archive &ar, Partitioning<D, Precision> &ps,
//                       __attribute((unused)) const unsigned int version) {
//            ar & boost::serialization::base_object<IndexedVector<Partition<D, Precision>>>(
//                    ps);
//        }
    }
}
