//File: Table.cpp
//Brief: A table of values converted to strings that can be printed to a
//       file or the terminal in markdown.  I like to turn its output into
//       a PDF with pandoc.
//Andrew Olivier aolivier@ur.rochester.edu

#ifndef UTIL_TABLE_CXX
#define UTIL_TABLE_CXX

#include "PlotUtils/Table.h"

namespace util
{
  template <int nCols>
  Table<nCols>::Table(const std::array<std::string, nCols>& colNames, const int nDigitsMax): fColMaxSizes{0}, fRows{colNames}, fNDigitsMax(nDigitsMax)
  {
    for(size_t whichCol = 0; whichCol < nCols; ++whichCol) fColMaxSizes[whichCol] = colNames[whichCol].length();
  }

  template <int nCols>
  std::ostream& Table<nCols>::print(std::ostream& destination) const
  {
    //Header
    destination << "|";
    for(size_t whichCol = 0; whichCol < nCols; ++whichCol)
    {
      destination << " " << std::setw(fColMaxSizes[whichCol]) << fRows[0][whichCol] << " |";
    }

    //Separator between header and entries 
    destination << "\n|" << std::setfill('-');
    for(size_t whichCol = 0; whichCol < nCols; ++whichCol)
    {
      destination << std::setw(fColMaxSizes[whichCol] + 3) << "|";
    }
    destination << std::setfill(' ');

    //Print each row
    for(auto row = fRows.cbegin() + 1; row != fRows.cend(); ++row)
    {
      destination << "\n|";
      for(size_t whichCol = 0; whichCol < nCols; ++whichCol)
      {
        destination << " " << std::setw(fColMaxSizes[whichCol]) << (*row)[whichCol] << " |";
      }
    }

    return destination;
  }
}

#endif //UTIL_TABLE_CXX
