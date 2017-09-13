#ifndef IVANP_OPT_PARSER_FWD_HH
#define IVANP_OPT_PARSER_FWD_HH

#if __has_include(<boost/lexical_cast/try_lexical_convert.hpp>)
#define PROGRAM_OPTIONS_BOOST_LEXICAL_CAST
#include <boost/lexical_cast/try_lexical_convert.hpp>
#else
#include <sstream>
#endif

namespace ivanp { namespace po {
namespace detail {

void arg_parser_impl_bool(const char* arg, bool& var);

}
}}

#endif
