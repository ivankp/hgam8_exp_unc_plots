#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <cmath>

#include "string.hh"
#include "math.hh"

using std::cout;
using std::cerr;
using std::endl;
using namespace ivanp;

int main() {
  std::ifstream hepdata("HGamEFTScanner/ATLAS_Run2_v2.HepData");

  struct bin {
    double min, max, xsec, stat;
    std::map<std::string,double> unc;
  };
  std::map<std::string,std::vector<bin>> vars;

  auto it = vars.end();
  std::string line, tok;
  unsigned line_n = 0;
  while (std::getline(hepdata,line)) {
    ++line_n;
    if (it==vars.end()) {
      if (starts_with(line,"*dataset:")) {
        const auto emp = vars.emplace(std::piecewise_construct,
          std::forward_as_tuple(line.substr(line.rfind('/')+1)),
          std::tie()
        );
        if (!emp.second) {
          cerr << "repeated variable: " << emp.first->first << endl;
          continue;
        }
        it = emp.first;
      }
    } else {
      const bool star = starts_with(line,"*");
      if (it->second.size()==0 && star) continue;
      if (line.size() && !star) {
        it->second.emplace_back();
        bin& b = it->second.back();

        const auto d1 = line.find(';');
        tok = line.substr(0,d1);
        if (starts_with(tok,">=")) tok.erase(0,2);
        std::stringstream ss(tok);
        ss >> b.min >> tok >> b.max;
        if (tok!="TO") b.max = b.min;
        
        const auto d2 = line.find('(',d1+1);
        ss.clear();
        ss.str(line.substr(d1+1,d2-d1-1));
        ss >> b.xsec >> tok >> b.stat;
        if (tok!="+-") throw std::runtime_error("missing +- in bin line");

        const char* cstr = line.c_str();
        for (size_t i=d2+1; line[i]!=';'; ) {
          // cout << cstr+i << endl;
          if (!starts_with(cstr+i,"DSYS")) throw std::runtime_error(
            "missing DSYS in bin line");
          const auto eq  = line.find('=',i);
          const auto col = line.find(':',eq+1);
          auto end = line.find(',',col+1);
          if (end==std::string::npos) end = line.find(')',col+1);
          const size_t sep = std::find(cstr+eq+1,cstr+col,',')-cstr;

          const auto unc = line.substr(col+1,end-col-1);
          if (!b.unc.emplace(
            unc,
            sep==col // one value
            ? std::stod(line.substr(eq+1,col-eq-1))
            : std::max( std::abs(std::stod(line.substr(eq +1,sep-eq -1))),
                        std::abs(std::stod(line.substr(sep+1,col-sep-1))) )
          ).second) throw std::runtime_error(cat(
            "duplicate uncert source \'",unc,"\' on line ",line_n));

          i = end+1;
        }

      } else it = vars.end();
    }
  }

  std::ofstream bands("bands.dat");

  for (auto& var : vars) {
    cout << var.first << endl;
    bands << "var " << var.first << endl;
    for (auto& b : var.second) {
      if (b.min==-100) continue;
      cout << b.min << ':' << b.max << ' ' << b.xsec << "±" << b.stat << endl;
      const double lumi = b.unc.at("lumi");
      b.unc.erase("lumi");
      cout << " lumi=" << lumi;
      const double fit = b.unc.at("fit"); // signal extraction
      b.unc.erase("fit");
      cout << " fit=" << fit;
      double corr_fac = 0.;
      for (const auto& unc : b.unc) {
        cout << ' ' << unc.first << '=' << unc.second;
        corr_fac += sq(unc.second);
      }
      corr_fac = std::sqrt(corr_fac);
      cout << endl;

      cout << lumi <<' '<< fit <<' '<< corr_fac <<' '<< b.stat << endl;

      bands << b.min
            << ' ' << qadd(lumi,corr_fac,fit,b.stat)/b.xsec
            << ' ' << qadd(lumi,corr_fac,fit)/b.xsec
            << ' ' << qadd(lumi,corr_fac)/b.xsec
            << ' ' << lumi/b.xsec
            << endl;
    }
    bands << var.second.back().max << " 0 0 0 0" << endl;
    bands << endl;
    cout << endl;
  }

  return 0;
}
