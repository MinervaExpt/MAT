//File: Table.h
//Brief: A square table to be printed in markdown.  Use pandoc to
//       turn it into a PDF like this: pandoc -s -o table.pdf table.md
//       where table.md is the output from Table::print(std::cout).
//       Works with arbitrary ostream objects, so it can be print()ed
//       to a file too with ofstream.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef UTIL_TABLE_H
#define UTIL_TABLE_H
#ifndef __GCCXML__
//c++ includes
#include <sstream>
#include <array>
#include <ostream>
#include <vector>
#include <iomanip>

namespace util
{
  //A "2-pass" table that holds its data in memory as strings
  template <int nCols>
  class Table
  {
    public:
      Table(const std::array<std::string, nCols>& colNames, const int nDigitsMax = 5);

      //Add a row to this Table
      template <class ...ARGS>
      int appendRow(ARGS... args)
      {
        const std::array<std::string, nCols> newRow{stringify(args)...};
        for(size_t whichCol = 0; whichCol < nCols; ++whichCol)
        {
          fColMaxSizes[whichCol] = std::max(fColMaxSizes[whichCol], newRow[whichCol].length());
        }
        fRows.push_back(newRow);
        return fRows.size();
      }

      //Print this table to the terminal or a file
      std::ostream& print(std::ostream& destination) const;

    private:
      std::array<size_t, nCols> fColMaxSizes;
      std::vector<std::array<std::string, nCols>> fRows;
      int fNDigitsMax; //Maximum number of digits that will be displayed for an entry.  Think of std::setprecision.

      template <class T>
      std::string stringify(const T& value)
      {
        std::stringstream ss;
        ss << std::setprecision(fNDigitsMax) << value;
        return ss.str();
      }
  };
}
//The "Ben Maneuver": make a class template look more familiar by putting its implementation in a different file
#include "Table.cxx"
#endif //__GCCXML__

#endif //UTIL_TABLE_H
