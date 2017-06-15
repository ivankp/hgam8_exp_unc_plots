#ifndef IVANP_ALGEBRA_HH
#define IVANP_ALGEBRA_HH

#include <utility>
#include <tuple>
#include <array>
#include <vector>
#include <iterator>

#include "math.hh"

namespace ivanp { namespace math {

// ==================================================================

template <typename B>
constexpr bool increment_ranges(const B&) { return true; }
template <typename R1, typename R2, typename B>
bool increment_ranges(const B& begins, std::pair<R1,R2>& r) {
  ++r.first;
  if (r.first == r.second) return true;
  return false;
}
template <typename R1, typename R2, typename... RR1, typename... RR2, typename B>
bool increment_ranges(const B& begins, std::pair<R1,R2>& r, std::pair<RR1,RR2>&... rr) {
  ++r.first;
  if (r.first == r.second) {
    r.first = std::get<std::tuple_size<B>::value-sizeof...(rr)-1>(begins);
    return increment_ranges(begins,rr...);
  }
  return false;
}

template <typename F, typename... R1, typename... R2,
          typename Ret = decltype(std::declval<F>()(std::declval<R1>()...))>
auto cartesian_product(F&& f, std::pair<R1,R2>... ranges)
-> std::enable_if_t<std::is_void<Ret>::value,void>
{
  const auto begins = std::make_tuple(ranges.first...);
  for (;;) {
    f(*ranges.first...); // apply
    if (increment_ranges(begins, ranges...)) break;
  }
}

template <typename F, typename... R1, typename... R2,
          typename Ret = decltype(std::declval<F>()(std::declval<R1>()...))>
auto cartesian_product(F&& f, std::pair<R1,R2>... ranges)
-> std::enable_if_t<!std::is_void<Ret>::value,std::vector<Ret>>
{
  std::vector<Ret> ret;
  ret.reserve(sum(std::distance(ranges.first,ranges.second)...));
  const auto begins = std::make_tuple(ranges.first...);
  for (;;) {
    ret.emplace_back(f(*ranges.first...)); // apply
    if (increment_ranges(begins, ranges...)) break;
  }
  return ret;
}

// ==================================================================

template <typename F, typename InputIt1, typename InputIt2, typename... InputIts,
          typename Ret = decltype(std::declval<F>()(
            *std::declval<InputIt1>(),*std::declval<InputIts>()...))>
auto direct_product(F&& f, InputIt1 first, InputIt2 last, InputIts... firsts)
-> std::enable_if_t<std::is_void<Ret>::value,void>
{
  for (; first!=last; ++first) f(*first,*(firsts++)...);
}

template <typename F, typename InputIt1, typename InputIt2, typename... InputIts,
          typename Ret = decltype(std::declval<F>()(
            *std::declval<InputIt1>(),*std::declval<InputIts>()...))>
auto direct_product(F&& f, InputIt1 first, InputIt2 last, InputIts... firsts)
-> std::enable_if_t<!std::is_void<Ret>::value,std::vector<Ret>>
{
  std::vector<Ret> ret;
  ret.reserve(std::distance(first,last));
  for (; first!=last; ++first)
    ret.emplace_back(f(*first,*(firsts++)...));
  return ret;
}

namespace detail {

template <typename... Args, typename Pred, size_t... I>
inline auto apply_direct_product(const std::tuple<Args...>& args, Pred&& f,
  std::index_sequence<I...>
) {
  return direct_product(std::forward<Pred>(f),
    std::get<0>(args).begin(), std::get<0>(args).end(),
    std::get<I+1>(args).begin()...
  );
}

}

// ==================================================================

template <typename T, size_t N, typename Pred, size_t... I>
inline auto map(const std::array<T,N>& in, Pred f,
                            std::index_sequence<I...>)
-> std::array<decltype(f(std::declval<T>())),N> {
  return { f(std::get<I>(in))... };
}
template <typename T, size_t N, typename Pred>
inline auto map(const std::array<T,N>& in, Pred f) {
  return map(in,f,std::make_index_sequence<N>{});
}

template <typename... T, typename Pred, size_t... I>
inline auto map(const std::tuple<T...>& in, Pred f,
                            std::index_sequence<I...>
) {
  return std::make_tuple(f(std::get<I>(in))...);
}
template <typename... T, typename Pred>
inline auto map(const std::tuple<T...>& in, Pred f) {
  return map(in,f,std::index_sequence_for<T...>{});
}

template <typename Cont, typename Pred>
auto map(const Cont& in, Pred f) {
  std::vector<std::decay_t<decltype(f(*std::begin(in)))>> out;
  out.reserve(in.size());
  for (const auto& x : in) out.emplace_back(f(x));
  return std::move(out);
}

// ==================================================================

}} // end namespace

template <typename T, typename Pred>
auto operator|(const std::vector<T>& in, Pred f) {
  std::vector<std::decay_t<decltype(f(*std::begin(in)))>> out;
  out.reserve(in.size());
  for (const auto& x : in) out.emplace_back(f(x));
  return std::move(out);
}

template <typename... Args, typename Pred>
inline auto operator*(const std::tuple<Args...>& args, Pred&& f) {
  return ivanp::math::detail::apply_direct_product(
    args, std::forward<Pred>(f),
    std::make_index_sequence<sizeof...(Args)-1>{});
}

// template <typename Cont, typename Pred>
// inline auto operator|(const Cont& in, Pred f) { ivanp::math::map(in,f); }

#endif
