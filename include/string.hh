#ifndef IVANP_STRING_HH
#define IVANP_STRING_HH

#include <string>
#include <sstream>
#include <utility>

namespace ivanp {

namespace detail {

template <typename S, typename T>
inline void cat_impl(S& s, const T& t) { s << t; }

template <typename S, typename T, typename... TT>
inline void cat_impl(S& s, const T& t, const TT&... tt) {
  s << t;
  cat_impl(s,tt...);
}

}

template <typename... TT>
inline std::string cat(const TT&... tt) {
  std::stringstream ss;
  detail::cat_impl(ss,tt...);
  return ss.str();
}

template <typename Str, unsigned N>                                             
bool starts_with(const Str& str, const char(&prefix)[N]) {                      
  for (unsigned i=0; i<N-1; ++i)                                                
    if (str[i]=='\0' || str[i]!=prefix[i]) return false;                        
  return true;                                                                  
}

}

#endif

