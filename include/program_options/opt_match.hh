#ifndef IVANP_OPT_MATCH_HH
#define IVANP_OPT_MATCH_HH

#ifdef PROGRAM_OPTIONS_STD_REGEX
#include <regex>
#elif defined(PROGRAM_OPTIONS_BOOST_REGEX)
#include <boost/regex.hpp>
#endif

#include <sstream>

namespace ivanp { namespace po {
namespace detail {

// Matchers ---------------------------------------------------------
// These represent rules for matching program arguments with argument
// definitions

struct opt_match_base {
  virtual bool operator()(const char* arg) const noexcept = 0;
  virtual ~opt_match_base() { }
  virtual std::string str() const noexcept = 0;
};

template <typename T> using oss_det = decltype(
  std::declval<std::ostringstream&>() << std::declval<const T&>() );
template <typename U> constexpr bool can_print =
  is_detected_v<oss_det,U> &&
  // this is here because lambdas without closure can convert to bool
  !(std::is_convertible<U,bool>::value && !std::is_pointer<U>::value);

template <typename T>
class opt_match final : public opt_match_base {
  T m; // matching rule
  template <typename U = T>
  inline std::enable_if_t<!can_print<U>,std::string>
  str_impl() const noexcept { return "?"; }
  // { return cat('[',type_str<U>(),']'); }
  template <typename U = T>
  inline std::enable_if_t<can_print<U>,std::string>
  str_impl() const noexcept {
    std::ostringstream ss;
    ss << m;
    return ss.str();
  };
public:
  template <typename... Args>
  opt_match(Args&&... args): m(std::forward<Args>(args)...) { }
  inline bool operator()(const char* arg) const noexcept { return m(arg); }
  inline std::string str() const noexcept { return str_impl(); }
};

template <>
inline bool opt_match<char>::operator()(const char* arg) const noexcept {
  return arg[1]==m;
}
template <>
inline std::string opt_match<char>::str() const noexcept {
  return {'-',m,'\0'};
}

template <> // defined in .cc
bool opt_match<const char*>::operator()(const char* arg) const noexcept;
template <>
inline std::string opt_match<const char*>::str() const noexcept { return m; }

template <>
inline bool opt_match<std::string>::operator()(const char* arg) const noexcept {
  return m == arg;
}
template <>
inline std::string opt_match<std::string>::str() const noexcept { return m; }

#ifdef _GLIBCXX_REGEX
template <>
inline bool opt_match<std::regex>::operator()(const char* arg) const noexcept {
  return std::regex_match(arg,m);
}
#endif

#ifdef BOOST_RE_REGEX_HPP
template <>
inline bool opt_match<boost::regex>::operator()(const char* arg) const noexcept {
  return boost::regex_match(arg,m);
}
#endif

// Argument type ----------------------------------------------------

enum opt_type { long_opt, short_opt, context_opt };

opt_type get_opt_type(const char* arg) noexcept;
inline opt_type get_opt_type(const std::string& arg) noexcept {
  return get_opt_type(arg.c_str());
}

// Matcher factories ------------------------------------------------

using opt_match_type = std::pair<const opt_match_base*,opt_type>;
template <typename T> struct opt_match_tag { using type = T; };

template <typename T>
inline opt_match_type make_opt_match(T&& x) {
  return make_opt_match_impl(std::forward<T>(x),
    opt_match_tag<std::decay_t<T>>{});
}

template <typename T, typename Tag>
opt_match_type make_opt_match_impl(T&& x, Tag) {
  using type = typename Tag::type;
  return { new opt_match<type>( std::forward<T>(x) ), context_opt };
}
template <typename T>
opt_match_type make_opt_match_impl(T&& x, opt_match_tag<char>) noexcept {
  return { new opt_match<char>( x ), short_opt };
}
template <typename T, typename TagT>
std::enable_if_t<std::is_convertible<TagT,std::string>::value,opt_match_type>
make_opt_match_impl(T&& x, opt_match_tag<TagT>) {
  const opt_type t = get_opt_type(x);
#if defined(PROGRAM_OPTIONS_STD_REGEX) || defined(PROGRAM_OPTIONS_BOOST_REGEX)
  using regex_t =
# ifdef PROGRAM_OPTIONS_STD_REGEX
    std::regex;
# else
    boost::regex;
# endif
  if (t==context_opt)
    return { new opt_match<regex_t>( std::forward<T>(x) ), t };
  else
#endif
  if (t==short_opt) {
    if (x[2]!='\0') throw po::error(
      "short arg "+std::string(x)+" defined with more than one char");
    return { new opt_match<char>( x[1] ), t };
  } else {
    return { new opt_match<TagT>( std::forward<T>(x) ), t };
  }
}

}
}}

#endif
