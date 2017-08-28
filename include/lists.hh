#ifndef IVANP_LISTS_HH
#define IVANP_LISTS_HH

// #include "detect.hh"
#include "meta.hh"

template <class T, class... Args>
using has_emplace_back_t = decltype(
  std::declval<T>().emplace_back(std::declval<Args>()...));

template <class T, class... Args>
using has_emplace_front_t = decltype(
  std::declval<T>().emplace_front(std::declval<Args>()...));

template <class T, class... Args>
using has_emplace_t = decltype(
  std::declval<T>().emplace(std::declval<Args>()...));

// emplace back =====================================================

template <typename List, typename T, std::enable_if_t<
          ivanp::is_detected<has_emplace_back_t,List,T>::value >* = nullptr>
inline decltype(auto) operator<<(List&& list, T&& x) {
  list.emplace_back(std::forward<T>(x));
  return std::forward<List>(list);
}

template <typename List, typename T, std::enable_if_t<
        ! ivanp::is_detected<has_emplace_back_t,List,T>::value
       && ivanp::is_detected<has_emplace_t,List,T>::value >* = nullptr>
inline decltype(auto) operator<<(List&& list, T&& x) {
  list.emplace(std::forward<T>(x));
  return std::forward<List>(list);
}

// emplace front ====================================================

template <typename List, typename T, std::enable_if_t<
          ivanp::is_detected<has_emplace_front_t,List,T>::value >* = nullptr>
inline decltype(auto) operator>>(T&& x, List&& list) {
  list.emplace_front(std::forward<T>(x));
  return std::forward<List>(list);
}

template <typename List, typename T, std::enable_if_t<
        ! ivanp::is_detected<has_emplace_front_t,List,T>::value
       && ivanp::is_detected<has_emplace_t,List,
                      typename List::const_iterator,T>::value >* = nullptr>
inline decltype(auto) operator>>(T&& x, List&& list) {
  list.emplace(list.begin(),std::forward<T>(x));
  return std::forward<List>(list);
}

#endif
