// Written by Ivan Pogrebnyak

#ifndef IVANP_MATH_HH
#define IVANP_MATH_HH

template <typename T> [[ gnu::const ]]
constexpr auto sq(T x) noexcept { return x*x; }
template <typename T, typename... TT> [[ gnu::const ]]
constexpr auto sq(T x, TT... xx) noexcept { return sq(x)+sq(xx...); }

template <typename... TT> [[ gnu::const ]]
constexpr auto qadd(TT... xx) noexcept { return std::sqrt(sq(xx...)); }

constexpr auto prod() noexcept { return 1; }
template <typename T> [[ gnu::const ]]
constexpr T prod(T x) noexcept { return x; }
template <typename T, typename... TT> [[ gnu::const ]]
constexpr auto prod(T x, TT... xx) noexcept { return x*prod(xx...); }

constexpr auto sum() noexcept { return 0; }
template <typename T> [[ gnu::const ]]
constexpr T sum(T x) noexcept { return x; }
template <typename T, typename... TT> [[ gnu::const ]]
constexpr auto sum(T x, TT... xx) noexcept { return x+sum(xx...); }

constexpr bool eq() noexcept { return true; }
template <typename T> [[ gnu::const ]]
constexpr bool eq(const T&) noexcept { return true; }
template <typename T1, typename T> [[ gnu::const ]]
constexpr bool eq(const T& x1, const T& x) noexcept { return x1==x; }
template <typename T1, typename T, typename... TT> [[ gnu::const ]]
constexpr bool eq(const T& x1, const T& x, const TT&... xx) noexcept {
  return (x1==x) && eq(x1,xx...);
}


#endif
