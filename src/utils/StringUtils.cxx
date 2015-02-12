#include "Geometry.h"
#include "Partitioner.h"

#include <sstream>

template <uint D, typename Precision>
std::string to_string(const dPoint<D, Precision> &p) {
  std::stringstream o;
  o << p.id << "-[" << p.coords[0];
  for (uint i = 1; i < D; ++i)
    o << ", " << p.coords[i];
  o << "]";
  return o.str();
}

template <uint D, typename Precision>
std::string to_string(const dSimplex<D, Precision> &p) {
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

template <uint D, typename Precision>
std::string to_string(const dBox<D, Precision> &b) {
  std::stringstream o;
  o << "[" << b.low[0] << " - " << b.high[0];
  for (uint i = 1; i < D; ++i)
    o << ", " << b.low[i] << " - " << b.high[i];
  o << "]";
  return o.str();
}

template <uint D, typename Precision>
std::string to_string(const dSphere<D, Precision> &b) {
  std::stringstream o;
  o << "O = [" << b.center[0];
  for (uint i = 1; i < D; ++i)
    o << ", " << b.center[i];
  o << "] - r = " << b.radius;
  return o.str();
}

template <uint D, typename Precision>
std::string to_string(const Partition<D, Precision> &p) {
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

template <uint D, typename Precision>
std::ostream &operator<<(std::ostream &o, const dPoint<D, Precision> &p) {
  return o << to_string(p);
}

template <uint D, typename Precision>
std::ostream &operator<<(std::ostream &o, const dSimplex<D, Precision> &p) {
  return o << to_string(p);
}

template <uint D, typename Precision>
std::ostream &operator<<(std::ostream &o, const dBox<D, Precision> &b) {
  return o << to_string(b);
}

template <uint D, typename Precision>
std::ostream &operator<<(std::ostream &o, const dSphere<D, Precision> &b) {
  return o << to_string(b);
}

template <uint D, typename Precision>
std::ostream &operator<<(std::ostream &o, const Partition<D, Precision> &p) {
  return o << to_string(p);
}

// specializations

// 2D
template std::string to_string(const dPoint<2, float> &p);
template std::string to_string(const dSimplex<2, float> &p);
template std::string to_string(const dBox<2, float> &p);
template std::string to_string(const dSphere<2, float> &p);
template std::string to_string(const Partition<2, float> &p);

template std::ostream &operator<<(std::ostream &o, const dPoint<2, float> &p);
template std::ostream &operator<<(std::ostream &o, const dSimplex<2, float> &p);
template std::ostream &operator<<(std::ostream &o, const dSphere<2, float> &p);
template std::ostream &operator<<(std::ostream &o, const dBox<2, float> &p);
template std::ostream &operator<<(std::ostream &o,
                                  const Partition<2, float> &p);

// 3D
template std::string to_string(const dPoint<3, float> &p);
template std::string to_string(const dSimplex<3, float> &p);
template std::string to_string(const dBox<3, float> &p);
template std::string to_string(const dSphere<3, float> &p);
template std::string to_string(const Partition<3, float> &p);

template std::ostream &operator<<(std::ostream &o, const dPoint<3, float> &p);
template std::ostream &operator<<(std::ostream &o, const dSimplex<3, float> &p);
template std::ostream &operator<<(std::ostream &o, const dBox<3, float> &p);
template std::ostream &operator<<(std::ostream &o, const dSphere<3, float> &p);
template std::ostream &operator<<(std::ostream &o,
                                  const Partition<3, float> &p);
