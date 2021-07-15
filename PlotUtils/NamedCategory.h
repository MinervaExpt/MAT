//==============================================================================
// Basically just a pair<T,string> class, with the possibility of expanding.
// Used as an input to the HistFolio.
//==============================================================================
#ifndef NAMEDCATEGORY_H
#define NAMEDCATEGORY_H

namespace MAT {
template <class value_t>
struct NamedCategory {
  constexpr NamedCategory(const value_t vals, const char* name)
      : m_value(vals), m_name(name) {}

  const value_t m_value;

  const std::string m_name;

  // Need a comparison operator to insert into maps.
  // This restricts the types that can be used as a template param to ones
  // that have '<' defined.
  bool operator<(const NamedCategory& rhs) const {
    return (m_value < rhs.m_value);
  }
};
}  // namespace MAT

#endif  // NAMEDCATEGORY_H
