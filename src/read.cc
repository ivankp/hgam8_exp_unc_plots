#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <stdexcept>

#include "program_options.hh"

#define TEST(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;
using ivanp::cat;

template <size_t N> // find N delimeters
inline auto findn(const std::string& s, char c, size_t pos=0) {
  std::array<size_t,N> ps{};
  for (size_t i=0; i<N; ++i) {
    pos = s.find(c,pos);
    if (pos==std::string::npos) break;
    ps[i] = pos;
    ++pos;
  }
  return ps;
}

template <typename T>
struct ref {
  const T* p;
  inline bool operator==(const ref& r) const { return *p == *r.p; }
  inline bool operator<(const ref& r) const { return *p < *r.p; }
  inline const T& operator*() const { return *p; }
  inline const T* operator->() const { return p; }
};
template <typename T>
inline std::ostream& operator<<(std::ostream& os, const ref<T>& rf) {
  return os << *rf;
}

template <typename C> // print values in container to string
inline std::string cont_str(const C& cont) {
  std::ostringstream ss;
  for (const auto& x : cont) ss << ' ' << x;
  return ss.str();
}

int main(int argc, char* argv[]) {
  const char* data_file_name;
  bool no_warnings = false,
       prt_bins = false, prt_modes = false, prt_vals = false;
  std::vector<const char*> vals;

  try {
    using namespace ivanp::po;
    if (program_options()
      (data_file_name,'f',"",req(),pos(1))
      (vals,{"-v","--vals"},"",pos(),multi())
      (prt_bins,"--prt-bins")
      (prt_modes,"--prt-modes")
      (prt_vals,"--prt-vals")
      (no_warnings,"--no-warnings")
      .parse(argc,argv,true)) return 0;
  } catch (const std::exception& e) {
    cerr <<"\033[31m"<< e.what() <<"\033[0m"<< endl;
    return 1;
  }

  std::unordered_map<
    std::string,
    std::unordered_map<
      std::string,
      std::unordered_map<
        std::string,
        std::vector<double>
      >
    >
  > data;

  // ================================================================
  { std::ifstream data_file(data_file_name);
  for (std::string line; std::getline(data_file,line); ) {
    static size_t line_i = 0;
    ++line_i; // count lines

    // skip blank lines and comments
    if (line.size()==0) continue;
    unsigned not_space = 0;
    while (std::isspace(line[not_space])) ++not_space;
    if (line.size()==not_space || line[not_space]=='#') continue;

    // find delimiters
    const auto d1s = findn<2>(line,'.',not_space);
    const auto d2s = findn<1>(line,':',d1s.back()+1);
    try {
      for (auto d : d1s) if (d==0) throw 0;
      if (d2s[0]==0) throw 0;
    } catch(...) {
      if (!no_warnings)
        cerr << "\033[33mLine " << line_i
             << ": unexpected formatting:\033[0m\n"
             << line << endl;
      continue;
    }

    // organize values in maps
    auto& var = data[line.substr(d1s[0]+1,d1s[1]-d1s[0]-1)];
    auto& mode = var[line.substr(0,d1s[0])];
    auto& val = mode[line.substr(d1s[1]+1,d2s[0]-d1s[1]-1)];

    if (val.size()) {
      if (!no_warnings)
        cerr << "\033[33mLine " << line_i
             << ": duplicate entry for:\033[0m\n"
             << line.substr(0,d2s[0]) << endl;
      continue;
    }

    std::istringstream ss(line.substr(d2s[0]+1));
    for (double x; ss >> x; ) val.push_back(x);

  }} // end lines loop
  // ================================================================

  // check binning for consistency ----------------------------------
  std::unordered_map<
    std::string,
    const std::vector<double>*
  > bins;

  for (const auto& var : data) {
    static const std::vector<double>* v0 = nullptr;
    for (const auto& mode : var.second) {
      const auto* v = &mode.second.at("bins");
      if (v0) {
        // compare vectors
        if (*v != *v0) {
          cerr << "\033[31mInconsistent binning at:\033[0m "
               << mode.first << '.' << var.first << ".bins" << endl;
          return 1;
        }
      } else v0 = v;
    }
    bins[var.first] = v0;
    v0 = nullptr;
  }

  if (prt_bins) { // option to print bins
    for (const auto& var : bins) {
      cout << var.first << ':';
      for (const auto& x : *var.second) {
        cout << ' ' << x;
      }
      cout << '\n';
    }
    cout << endl;
  }

  // check modes for consistency ------------------------------------
  std::set<ref<std::string>> modes;

  for (const auto& var : data) {
    static const auto* var1 = &var.first;
    decltype(modes) m;
    for (const auto& mode : var.second) {
      m.insert({&mode.first});
    }
    if (modes.size()) {
      if (m != modes) {
        cerr << "\033[31mInconsistent modes:\033[0m\n"
             << *var1 << ": " << cont_str(modes) << '\n'
             << var.first << ": " << cont_str(m) << endl;
        return 1;
      }
    } else modes = std::move(m);
  }

  if (prt_modes) {
    for (const auto& m : modes) cout << m << '\n';
    cout << endl;
  }

  // print values ---------------------------------------------------
  if (prt_vals) {
    std::set<ref<std::string>> vals_set;

    for (const auto& var : data) {
      for (const auto& mode : var.second)
        for (const auto& val : mode.second)
          vals_set.insert({&val.first});
    }

    std::vector<ref<std::string>> vals;
    vals.reserve(vals_set.size());
    for (auto& v : vals_set) vals.push_back(v);

    std::sort(vals.begin(),vals.end());
    for (auto& v : vals) cout << v << '\n';
    cout << endl;
  }

  // sum over production modes --------------------------------------
  if (vals.size()) {

    std::unordered_map<
      std::string,
      std::unordered_map<
        std::string,
        std::vector<double>
    >> sums;
    for (auto v : vals) sums[v];

    for (const auto& var : data) {
      for (auto& sum : sums) {
        auto& xs = sum.second[var.first];
        for (const auto& mode : var.second) {
          try {
            const auto& v = mode.second.at(sum.first);
            const auto n = xs.size();
            if (n) {
              if (v.size()!=n) {
                cerr << "\033[31mUnequal number of sumues for:\033[0m "
                     << mode.first << '.' << var.first << '.' << sum.first
                     << endl;
                return 1;
              }
              for (unsigned i=0; i<n; ++i) xs[i] += v[i];
            } else xs = v;
          } catch(const std::out_of_range& e) {
            cerr << "\033[31m" << mode.first << '.' << var.first
                 << " has no value " << sum.first << "\033[0m" << endl;
            return 1;
          }
        }
      }
    }

    for (const auto& val : sums) {
      cout << "\033[0;1m" << val.first << "\033[0m\n";
      for (const auto& var : val.second) {
        cout << "  " << var.first;
        for (const auto& x : var.second) cout << ' ' << x;
        cout << '\n';
      }
    }
    cout << endl;

  }

}
