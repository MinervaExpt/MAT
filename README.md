# MINERvA Analysis Toolkit

An error analysis toolkit for particle physics described in [Arxiv 2103.08677](https://arxiv.org/abs/2103.08677).  Provides boilerplate for organizing systematic uncertainties so they can all be treated on the same footing during the event loop stage.  Systematics-aware histograms like `MnvH1D` propagate uncertainties through addition and division so that calculations like cross sections can be easily treated with the **many universe** systematic uncertainty method.

## Installation
If you are a MINERvA collaborator, just install [MAT-MINERvA](https://github.com/MinervaExpt/MAT-MINERvA) instead.  It installs this package too and keeps their builds synchronized.

### CMake build (for non-MINERvA experiments)
Requires at least CMake 3.

```
git clone https://github.com/MinervaExpt/MAT.git
mkdir -p opt/build && cd opt/build
cmake ../../MAT -DCMAKE_INSTALL_PREFIX=`pwd`/.. -DCMAKE_BUILD_TYPE=Release
make install #-j 4
```

## How to use MAT

### Learn to use MAT systematics

MINERvA's data preservation project provides a complete example of how to do an analysis using the MINERvA Analysis Toolkit:
- Install this package
- [MAT-MINERvA]() extends the MAT with MINERvA-specific systematic uncertainties
- The [2021 MINERvA 101 Tutorial](https://github.com/MinervaExpt/MINERvA-101-Cross-Section) walks you through a simple inclusive cross section analysis in detail

### Use it in your own analysis
- CMake: Works automatically if you install in the same `opt` prefix.  Just create a separate `buildYourPackage` area under `opt`.
- Makefile: Your package can link against libMAT-MINERvA in `opt/lib`.  Headers to include are in `opt/include`.
- PyROOT: MAT-MINERvA automatically generates python bindings thanks to ROOT.  Use the `.rootlogon.C` in the interpreter instructions below, and `from ROOT import PlotUtils`.
- ROOT's CLING Interpreter: Install this `.rootlogon.C` in either your home area or your working area
```
{
  if( gSystem->Getenv("PLOTUTILSROOT") )
  {
    string newpath = string(gROOT->GetMacroPath()) + ":" + string("${PLOTUTILSROOT}/../bin" );
    gROOT->SetMacroPath( newpath.c_str() );
    gInterpreter->AddIncludePath( "${PLOTUTILSROOT}/../include" );
    gInterpreter->AddIncludePath( "${PLOTUTILSROOT}/../include/PlotUtils" );
    std::vector<std::string> packages = { "MAT", "MAT-MINERvA" };
    for(const std::string& package: packages)
    {
      gSystem->Load( gSystem->ExpandPathName(("$PLOTUTILSROOT/lib" + package + ".so").c_str()) );
    }
  }
}
```



## Contributing

### New Adopters:
Feel free to fork this package.  If you want to contribute to the MAT in its shared form, please contact the MINERvA experiment to get write permissions.  We may also consider pull requests through github, but we'd rather keep contributing as simple for everyone as possible.

### MINERvA Analyzers:
- Commit directly to the main branch for now, effectively the way CVS used to work.  If we get more breaking changes than we can handle, the respository maintainers will move to a pull-request-powered contribution workflow.
- Use `git merge` to fold in changes to your working area and other branches.  We're keeping the commit procedure as simple as possible at the cost of a little commit history bloat.
- For big changes that you want to collaborate on, create a new branch with "feature/" at the beginning of its name.
