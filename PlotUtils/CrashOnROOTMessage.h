//File: CrashOnROOTMessage.h
//Brief: Include this header to cause your jobs to crash whenever ROOT
//       would have normally printed a message about a warning or an
//       error.  This actually throws a ROOT::exception that you can
//       catch.  For an expert-level interface for reacting to different
//       kinds of messages, see ErrorHandler.h
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef CRASHONROOTMESSAGE_H
#define CRASHONROOTMESSAGE_H

//Plotutils includes
#include "PlotUtils/ErrorHandler.h"

//Force override of ROOT's default error handler by just including this header
ROOT::detail::beforeMain forceErrorHandlerOverride;

#endif //CRASHONROOTMESSAGE_H
