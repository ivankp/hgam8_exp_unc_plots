#include <iostream>
#include <cstring>
#include <stdexcept>

#include "program_options.hh"

using std::cout;
using std::cerr;
using std::endl;
using namespace std::string_literals;

bool is_number(const char* str) noexcept {
// https://stackoverflow.com/q/4654636/2640636
#ifdef PROGRAM_OPTIONS_BOOST_LEXICAL_CAST
  static double d;
  return boost::conversion::try_lexical_convert(str,d);
#else
  char* p;
  std::strtod(str,&p);
  return !*p;
#endif
}

namespace ivanp { namespace po {

namespace detail {

template <>
bool opt_match<const char*>::operator()(const char* arg) const noexcept {
  int i = 0;
  for (; m[i]!='\0' && arg[i]!='\0'; ++i)
    if ( arg[i]!=m[i] ) return false;
  return m[i]=='\0' && arg[i]=='\0';
}

opt_type get_opt_type(const char* arg) noexcept {
  unsigned char n = 0;
  for (char c=arg[n]; c=='-'; c=arg[++n]) ;
  switch (n) {
    case  1: return is_number(arg) ? context_opt : short_opt;
    case  2: return long_opt;
    default: return context_opt;
  }
}

}

void check_count(detail::opt_def* opt) {
  if (!opt->is_multi() && opt->count)
    throw error("too many options " + opt->name);
}

void program_options::parse(int argc, char const * const * argv) {
  using namespace ::ivanp::po::detail;

  opt_def *opt = nullptr;
  const char* val = nullptr;
  std::string tmp; // use view here
  bool last_was_val = false;

  for (int i=1; i<argc; ++i) {
    const char* arg = argv[i];
    last_was_val = false;

    const auto opt_type = get_opt_type(arg);

    // ==============================================================

    if (opt_type!=context_opt) {
      if (opt) {
        if (!opt->count) opt->as_switch();
        opt = nullptr;
      }
      if (opt_type==long_opt) { // long: split by '='
        if ((val = strchr(arg,'='))) arg = tmp.assign(arg,val).c_str(), ++val;
      } else { // short: allow spaceless
        if (arg[2]!='\0') val = arg+2;
      }
    }

    // ==============================================================

    if (!opt || (opt->is_multi() && opt->count)) {
      for (auto& m : matchers[opt_type]) {
        if ((*m.first)(arg)) { // match
          opt = m.second;
          check_count(opt);
          if (opt_type==context_opt) val = arg;
          if (opt->is_switch()) {
            if (val) throw po::error(
              "switch " + opt->name + " does not take arguments");
            opt->as_switch(), opt = nullptr;
          } else if (val) {
            opt->parse(val), val = nullptr;
            last_was_val = true;
            if (!opt->is_multi()) opt = nullptr;
          }
          goto next_arg;
        }
      }
    }

    if (opt) {
      opt->parse(arg);
      last_was_val = true;
      if (!opt->is_multi()) opt = nullptr;
      continue;
    }

    // handle positional options
    if (opt_type==context_opt && !opt && pos.size()) {
      auto *pos_opt = pos.front();
      check_count(pos_opt);
      pos_opt->parse(arg);
      last_was_val = true;
      if (!pos_opt->is_pos_end()) pos.pop();
      continue;
    }

    throw po::error("unexpected option "s + arg);
    next_arg: ;
  } // end arg loop
  if (opt) {
    if (!opt->count) opt->as_switch();
    else if (!last_was_val) throw po::error("dangling option " + opt->name);
    opt = nullptr;
  }

  for (opt_def *opt : req) // check required passed
    if (!opt->count) throw error("missing required option "+opt->name);

  for (opt_def *opt : default_init) // init with default values
    if (!opt->count) opt->default_init();
}

inline bool streq_ignorecase(const char* str, const char* s1) {
  const char *a=str, *b=s1;
  for (; *a!='\0' && *b!='\0'; ++a, ++b)
    if (toupper(*a) != toupper(*b)) return false;
  return *a=='\0' && *b=='\0';
}
template <typename S1>
inline bool streq_any_ignorecase(const char* str, S1 s1) {
  return streq_ignorecase(str,s1);
}
template <typename S1, typename... Ss>
inline bool streq_any_ignorecase(const char* str, S1 s1, Ss... ss) {
  return streq_ignorecase(str,s1) || streq_any_ignorecase(str,ss...);
}

namespace detail {

template <>
void arg_parser_impl<bool>(const char* arg, bool& var) {
  if (streq_any_ignorecase(arg,"1","TRUE","YES","ON","Y")) var = true;
  else if (streq_any_ignorecase(arg,"0","FALSE","NO","OFF")) var = false;
  else throw po::error(cat('\"',arg,"\" cannot be interpreted as bool"));
}

}

}} // end namespace ivanp
