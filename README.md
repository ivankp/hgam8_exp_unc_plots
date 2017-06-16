Compilation: `make`

Executable: `bin/plot`

Program options:
1. The first option is always the input `HepData` file name
2. `corr` -- produce plots of 'correction factor' uncertainties.
   By default, plots showing lumi, correction factor, fit, and stat
   uncertainties are produced.
3. `burst` -- instead of a single file, output plots in individual files for
   each variable.

Output:
* Without `burst`, `uncert.pdf` file is produced.
* With `burst`, files named `VAR.pdf` are produced, where `VAR` is the name of
  the variable.
* With `corr`, suffix `_corr` is added to the file(s) name(s).

Usage examples:
```
./bin/plot HGamEFTScanner/ATLAS_Run2_v2.HepData burst
./bin/plot HGamEFTScanner/ATLAS_Run2_v2.HepData burst corr
```
