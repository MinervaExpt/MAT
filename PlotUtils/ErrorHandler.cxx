//File: ErrorHandler.cxx
//Brief: Override how ROOT handles errors so that the user has
//       the option to catch them as exceptions.
//Usage: #include "ErrorHandler.h"
//       //and link against the library you build this into
//Author: Andrew Olivier aolivier@ur.rochester.edu

//Access the python API to communicate with an already-setup interpreter.
//This has to be included before anything else to avoid warning messages
//about scary-looking preprocessor macros being redefined.
//#include "Python.h"
//Include the header
#include "ErrorHandler.h"

namespace MAT
{
  void HandleErrorsInPython()
  {
    //SetErrorHandler(ROOT::pythonErrorHandler); //This is hard to compile on OS X as of February 2021.  I think it has to do with needing python development headers.
    SetErrorHandler(ROOT::errorHandler);
  }
}

namespace ROOT
{
  //ErrorHandlerFunc_t to convert TObject::Warning() and friends into
  //exceptions.
  void errorHandler(int level, Bool_t abort, const char* location, const char* msg)
  {
    if(abort) throw error(level, location, msg);
    if(level >= gErrorIgnoreLevel) throw warning(level, location, msg);
  }

  //Throw a python warning or error instead of a c++ exception on any message from ROOT.
  /*void pythonErrorHandler(int level, Bool_t abort, const char* location, const char* msg)
  {
    if(abort)
    {
      PyErr_SetString(PyExc_SystemExit, (std::string(location) + ": " + msg).c_str());
    }
    else if(level >= gErrorIgnoreLevel)
    {
      PyErr_SetString(PyExc_RuntimeWarning, (std::string(location) + ": " + msg).c_str());
    }
  }*/
 
  exception::exception(const int howBad, const std::string& location, const char* what): std::runtime_error(location + ": " + what), where(location), level(howBad)
  {
  }

  error::error(const int level, const std::string& where, const char* what): exception(level, where, what)
  {
  }

  warning::warning(const int level, const std::string& where, const char* what): exception(level, where, what)
  {
  }

  namespace detail //i.e. don't you dare write code that uses this!
  {
    //Machinery to force the error handler to be overridden by just including
    //this file.
    //I got this idea from https://stackoverflow.com/questions/10897552/call-a-function-before-main
    /*struct beforeMain
    {
      beforeMain()
      {
        SetErrorHandler(errorHandler);
      }
    };*/

    beforeMain::beforeMain()
    {
      SetErrorHandler(errorHandler);
    }
  }
}

//ROOT::detail::beforeMain forceErrorHandling;
