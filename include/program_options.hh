#ifndef IVANP_PROGRAM_OPTIONS_HH
#define IVANP_PROGRAM_OPTIONS_HH

#include <string>
#include <vector>
#include <queue>
#include <array>
#include <memory>
#include <type_traits>
#include <stdexcept>

#include "type.hh"
#include "type_traits.hh"
#include "string.hh"
#include "tuple_alg.hh"
#include "seq_alg.hh"

namespace ivanp { namespace po {
struct error : std::runtime_error {
  using std::runtime_error::runtime_error;
};
}}

#include "program_options/opt_match.hh"
#include "program_options/opt_def.hh"

namespace ivanp { namespace po {

class program_options {
  std::vector<std::unique_ptr<detail::opt_def>> opt_defs;
  std::array<std::vector<std::pair<
    std::unique_ptr<const detail::opt_match_base>,
    detail::opt_def*
  >>,3> matchers;
  std::queue<detail::opt_def*> pos;
  std::vector<detail::opt_def*> req, default_init;

  template <typename T, typename... Props>
  inline auto* add_opt(T* x, std::string&& descr, Props&&... p) {
    static_assert( !std::is_const<T>::value,
      "\033[33mpointer to const in program option definition\033[0m");

    using props_types = std::tuple<std::decay_t<Props>...>;
    auto props = std::forward_as_tuple(std::forward<Props>(p)...);

#define UNIQUE_PROP_ASSERT(NAME) \
    using NAME##_i = indices_of<_::is_##NAME, props_types>; \
    static_assert( NAME##_i::size() <= 1, \
      "\033[33mrepeated \"" #NAME "\" in program option definition\033[0m");

    UNIQUE_PROP_ASSERT(switch_init)
    UNIQUE_PROP_ASSERT(default_init)
    UNIQUE_PROP_ASSERT(named)
    UNIQUE_PROP_ASSERT(pos)
    UNIQUE_PROP_ASSERT(npos)
    UNIQUE_PROP_ASSERT(req)
    UNIQUE_PROP_ASSERT(multi)

#undef UNIQUE_PROP_ASSERT

    using parser_i = indices_of<_::is_parser<T>::template type, props_types>;
    static_assert( parser_i::size() <= 1,
      "\033[33mrepeated parser in program argument definition\033[0m");

    static_assert( !(pos_i::size() && npos_i::size()),
      "\033[33monly one positional property can be specified\033[0m");

    using prop_seq = seq::join_t<
      parser_i,
      switch_init_i,
      default_init_i,
      named_i,
      pos_i,
      npos_i,
      req_i,
      multi_i
    >;

    static_assert( prop_seq::size() == sizeof...(Props),
      "\033[33munrecognized argument in program option definition;"
      " check parser call signature\033[0m");

    auto *opt = detail::make_opt_def(
      x, std::move(descr), std::move(props), prop_seq{});
    opt_defs.emplace_back(opt);

    opt->set_name(std::move(props),seq::head_t<named_i>{});

    if (pos_i::size() || npos_i::size()) {
      if (pos.size() && pos_i::size() && pos.back()->is_pos_end())
        throw error("only one indefinite positional option can be specified");
      pos.push(opt);
    }

    if (req_i::size()) req.push_back(opt);
    if (default_init_i::size()) default_init.push_back(opt);

    return opt;
  }

  template <typename Matcher>
  inline void add_match(
    Matcher&& matcher, detail::opt_def* opt
  ) {
    auto&& m = detail::make_opt_match(std::forward<Matcher>(matcher));
    matchers[m.second].emplace_back(std::move(m.first),opt);
    if (!opt->is_named()) {
      std::string& name = opt->name;
      if (name.size()) name += ',';
      name += m.first->str();
    }
  }
  template <typename... M, size_t... I>
  inline void add_matches(
    const std::tuple<M...>& matchers, detail::opt_def* opt,
    std::index_sequence<I...>
  ) {
#ifdef __cpp_fold_expressions
    (add_match(std::get<I>(matchers),opt),...);
#else
    using discard = const char[];
    (void)discard{(add_match(std::get<I>(matchers),opt),'\0')...};
#endif
  }

public:
  void parse(int argc, char const * const * argv);

  template <typename T, typename... Props>
  program_options& operator()(T* x,
    std::initializer_list<const char*> matchers,
    std::string descr={}, Props&&... p
  ) {
    auto *opt = add_opt(x,std::move(descr),std::forward<Props>(p)...);
    for (const char* m : matchers) add_match(m,opt);
    return *this;
  }

  template <typename T, typename Matcher, typename... Props>
  std::enable_if_t<!is_tuple<std::decay_t<Matcher>>::value,program_options&>
  operator()(T* x,
    Matcher&& matcher,
    std::string descr={}, Props&&... p
  ) {
    auto *opt = add_opt(x,std::move(descr),std::forward<Props>(p)...);
    add_match(matcher,opt);
    return *this;
  }

  template <typename T, typename... Matchers, typename... Props>
  program_options& operator()(T* x,
    const std::tuple<Matchers...>& matchers,
    std::string descr={}, Props&&... p
  ) {
    static_assert( sizeof...(Matchers) > 0,
      "\033[33mempty tuple in program argument definition\033[0m");
    auto *opt = add_opt(x,std::move(descr),std::forward<Props>(p)...);
    add_matches(matchers,opt,std::index_sequence_for<Matchers...>{});
    return *this;
  }
};

}} // end namespace ivanp

#endif
