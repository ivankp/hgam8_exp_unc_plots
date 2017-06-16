#ifndef IVANP_LISTS_HH
#define IVANP_LISTS_HH

#include "detect.hh"

template <class T, class... Args>
using has_emplace_back_t = decltype(
  std::declval<T>().emplace_back(std::declval<Args>()...));

template <class T, class... Args>
using has_emplace_t = decltype(
  std::declval<T>().emplace(std::declval<Args>()...));

template <typename List, typename T, std::enable_if_t<
          is_detected<has_emplace_back_t,List,T>::value >* = nullptr>
inline decltype(auto) operator<<(List&& list, T&& x) {
  list.emplace_back(std::forward<T>(x));
  return std::forward<List>(list);
}

template <typename List, typename T, std::enable_if_t<
        ! is_detected<has_emplace_back_t,List,T>::value
       && is_detected<has_emplace_t,List,T>::value >* = nullptr>
inline decltype(auto) operator<<(List&& list, T&& x) {
  list.emplace(std::forward<T>(x));
  return std::forward<List>(list);
}

#endif
