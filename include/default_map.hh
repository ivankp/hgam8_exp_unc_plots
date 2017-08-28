#ifndef IVANP_DEFAULT_MAP_HH
#define IVANP_DEFAULT_MAP_HH

template <typename Map, typename F=void>
class default_map: Map {
  F f;
public:
  using key_type = typename Map::key_type;
  using mapped_type = typename Map::mapped_type;
  using value_type = typename Map::value_type;

  mapped_type operator[](const key_type& key) const {
    try { return Map::at(key); } catch (...) { return f(key); }
  }

  using Map::emplace;

  default_map() = default;
  template <typename U>
  default_map(std::initializer_list<value_type> init, U&& f)
  : Map(init), f(std::forward<U>(f)) { }
};

template <typename Map>
class default_map<Map,void>: Map {
public:
  using key_type = typename Map::key_type;
  using mapped_type = typename Map::mapped_type;
  using value_type = typename Map::value_type;

  const mapped_type& operator[](const key_type& key) const {
    try { return Map::at(key); } catch (...) { return key; }
  }

  using Map::emplace;

  default_map() = default;
  default_map(std::initializer_list<value_type> init): Map(init) { }
};

template <typename Key, typename T, typename F>
inline default_map<std::unordered_map<Key,T>,F> make_default_map(
  std::initializer_list<std::pair<const Key, T>> init, F&& f
) {
  return { init, std::forward<F>(f) };
}
template <typename Key, typename T = Key>
inline default_map<std::unordered_map<Key,T>,void> make_default_map(
  std::initializer_list<std::pair<const Key, T>> init
) {
  return { init };
}

#endif
