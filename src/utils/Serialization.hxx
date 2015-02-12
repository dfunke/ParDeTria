#pragma once

#include <fstream>

// boost
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/base_object.hpp>

template <class Object> void storeObject(Object &o, const std::string &file) {

  std::ofstream f(file, std::ios::trunc | std::ios::out | std::ios::binary);
  boost::archive::binary_oarchive oa(f);
  oa << o;
}

template <class Object> Object loadObject(const std::string &file) {

  std::ifstream f(file, std::ios::in | std::ios::binary);
  boost::archive::binary_iarchive ia(f);

  Object o;
  ia >> o;

  return o;
}

namespace boost {
namespace serialization {

template <class Archive, uint D, typename Precision>
void serialize(Archive &ar, dVector<D, Precision> &v,
               __attribute((unused)) const unsigned int version) {
  ar &v;
}

template <class Archive, uint D, typename Precision>
void serialize(Archive &ar, dBox<D, Precision> &b,
               __attribute((unused)) const unsigned int version) {
  ar &b.coords;
  ar &b.dim;
}

template <class Archive, uint D, typename Precision>
void serialize(Archive &ar, dSphere<D, Precision> &s,
               __attribute((unused)) const unsigned int version) {
  ar &s.center;
  ar &s.radius;
}

template <class Archive, typename V, typename K>
void serialize(Archive &ar, IndexedVector<V, K> &v,
               __attribute((unused)) const unsigned int version) {
  ar &boost::serialization::base_object<std::map<K, V>>(v);
}

// points

template <class Archive, uint D, typename Precision>
void serialize(Archive &ar, dPoint<D, Precision> &p,
               __attribute((unused)) const unsigned int version) {
  ar &p.id;
  ar &p.coords;
  ar &p.simplices;
}

template <class Archive, uint D, typename Precision>
void serialize(Archive &ar, dPoints<D, Precision> &ps,
               __attribute((unused)) const unsigned int version) {
  ar &boost::serialization::base_object<IndexedVector<dPoint<D, Precision>>>(
      ps);
}

// simplices

template <class Archive, uint D, typename Precision>
void serialize(Archive &ar, dSimplex<D, Precision> &s,
               __attribute((unused)) const unsigned int version) {
  ar &s.id;
  ar &s.vertices;
  ar &s.neighbors;
}

template <class Archive, uint D, typename Precision>
void serialize(Archive &ar, dSimplices<D, Precision> &ss,
               __attribute((unused)) const unsigned int version) {
  ar &boost::serialization::base_object<IndexedVector<dSimplex<D, Precision>>>(
      ss);
}

// partitions

template <class Archive, uint D, typename Precision>
void serialize(Archive &ar, Partition<D, Precision> &p,
               __attribute((unused)) const unsigned int version) {
  ar &p.id;
  ar &p.bounds;
  ar &p.points;
}

template <class Archive, uint D, typename Precision>
void serialize(Archive &ar, Partitioning<D, Precision> &ps,
               __attribute((unused)) const unsigned int version) {
  ar &boost::serialization::base_object<IndexedVector<Partition<D, Precision>>>(
      ps);
}
}
}
