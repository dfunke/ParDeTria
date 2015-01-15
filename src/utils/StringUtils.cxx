#include "Geometry.h"
#include "Partitioner.h"

#include <sstream>

template <uint D> std::string to_string(const dPoint<D> &p) {
  std::stringstream o;
  o << p.id << "-[" << p.coords[0];
  for (uint i = 1; i < D; ++i)
    o << ", " << p.coords[i];
  o << "]";
  return o.str();
}

template <uint D> std::string to_string(const dSimplex<D> &p) {
  std::stringstream o;
  o << p.id << ": V = [" << p.vertices[0];
  for (uint i = 1; i < D + 1; ++i)
    o << ", " << p.vertices[i];
  o << "] N = [";
  for (auto it = p.neighbors.begin(); it != p.neighbors.end(); ++it) {
    if (it != p.neighbors.begin())
      o << ", ";
    o << *it;
  }
  o << "]";
  return o.str();
}

template <uint D> std::string to_string(const dBox<D> &b) {
  std::stringstream o;
  o << "[" << b.coords[0] << " - " << b.coords[0] + b.dim[0];
  for (uint i = 1; i < D; ++i)
    o << ", " << b.coords[i] << " - " << b.coords[i] + b.dim[i];
  o << "]";
  return o.str();
}

template <uint D> std::string to_string(const Partition<D> &p) {
  std::stringstream o;
  o << p.id << " " << p.bounds << " [";
  for (auto it = p.points.begin(); it != p.points.end(); ++it) {
    if (it != p.points.begin())
      o << ", ";
    o << *it;
  }
  o << "]";
  return o.str();
}

template <uint D>
std::ostream &operator<<(std::ostream &o, const dPoint<D> &p) {
  return o << to_string(p);
}

template <uint D>
std::ostream &operator<<(std::ostream &o, const dSimplex<D> &p) {
  return o << to_string(p);
}

template <uint D> std::ostream &operator<<(std::ostream &o, const dBox<D> &b) {
  return o << to_string(b);
}

template <uint D>
std::ostream &operator<<(std::ostream &o, const Partition<D> &p) {
  return o << to_string(p);
}

// specializations

// 2D
template std::string to_string(const dPoint<2> &p);
template std::string to_string(const dSimplex<2> &p);
template std::string to_string(const dBox<2> &p);
template std::string to_string(const Partition<2> &p);

template std::ostream &operator<<(std::ostream &o, const dPoint<2> &p);
template std::ostream &operator<<(std::ostream &o, const dSimplex<2> &p);
template std::ostream &operator<<(std::ostream &o, const dBox<2> &p);
template std::ostream &operator<<(std::ostream &o, const Partition<2> &p);

// 3D
template std::string to_string(const dPoint<3> &p);
template std::string to_string(const dSimplex<3> &p);
template std::string to_string(const dBox<3> &p);
template std::string to_string(const Partition<3> &p);

template std::ostream &operator<<(std::ostream &o, const dPoint<3> &p);
template std::ostream &operator<<(std::ostream &o, const dSimplex<3> &p);
template std::ostream &operator<<(std::ostream &o, const dBox<3> &p);
template std::ostream &operator<<(std::ostream &o, const Partition<3> &p);
