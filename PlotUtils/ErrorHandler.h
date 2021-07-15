//File: ErrorHandler.h
//Brief: Override how ROOT handles errors so that the user has
//       the option to catch them as exceptions.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef ROOT_ERRORHANDLER_H
#define ROOT_ERRORHANDLER_H

//ROOT includes
#include "TError.h" //For SetErrorHandler()

//c++ includes
#include <stdexcept>
#include <string>

//Exception specifiers, like throw() applied to a function, are deprecated in c++11,
//but Reflex, which only knows about c++03, uses a version of std::exception that
//has an exception speficier on its destructor.  This behavior is required of all
//destructors in c++11 anyway.  So, I'm using a preprocessor macro to only show this
//exception specification to the gccxml tool used by Reflex to build our ROOT dictionary.
//I'll hide it from "real" compilers like the one ROOT 6 might be using.
#ifdef __GCC_XML__
  #define DESTRUCTOR_NOEXCEPT throw()
#else
  #define DESTRUCTOR_NOEXCEPT
#endif

namespace MAT
{
  //Pyroot interface to register an error handler that throws python exceptions on ROOT warning messages
  void HandleErrorsInPython();
}

namespace ROOT
{
  //ErrorHandlerFunc_t to convert TObject::Warning() and friends into
  //exceptions.
  //N.B.: TError.h's default error handler is declared extern, but my
  //      Googling has revealed that this is the default behavior for functions
  //      in c++ anyway.  If it were extern "C", that would be another story
  //      entirely, but (the vast majority of) ROOT is linked as c++.
  void errorHandler(int level, Bool_t abort, const char* location, const char* msg);

  //Throw a python warning or error instead of a c++ exception on any message from ROOT.
  //void pythonErrorHandler(int /*level*/, Bool_t abort, const char* location, const char* msg);

  //Exception classes to let the user separate "warnings" and "errors" from ROOT.
  //They derive from a common base class to both group functionality and let the
  //user catch them together.
  class exception: public std::runtime_error
  {
    public:
      exception(const int level, const std::string& where, const char* what);
      virtual ~exception() DESTRUCTOR_NOEXCEPT {};

      const std::string where; //Nominally what class/function caused this
                               //exception.  We're at the mercy of ROOT though.
      const int level; //How severe the error message was?  It's the level
                       //parameter to ErrorHandlerFunc_t
  };

  //An error is a ROOT::exception on which ROOT tried to abort the program.
  class error: public exception
  {
    public:
      error(const int level, const std::string& where, const char* what);
      virtual ~error() DESTRUCTOR_NOEXCEPT {};
  };

  //A warning is a ROOT::exception that ROOT did not think was severe enough
  //to abort the program.
  class warning: public exception
  {
    public:
      warning(const int level, const std::string& where, const char* what);
      virtual ~warning() DESTRUCTOR_NOEXCEPT {};
  };

  //Machinery to force the error handler to be overridden by just including
  //this file.
  //I got this idea from https://stackoverflow.com/questions/10897552/call-a-function-before-main
  namespace detail //i.e. don't you dare write code that uses this!
  {
    struct beforeMain
    {
      beforeMain();
    };
  }
}

#endif //ROOT_ERRORHANDLER_H
