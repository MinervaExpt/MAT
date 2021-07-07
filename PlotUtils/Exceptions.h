/*
 * Exceptions.h
 *
 *  Created on: Jul 30, 2014
 *      Author: J. Wolcott <jwolcott@fnal.gov>
 */

#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_

#include <stdexcept>

namespace PlotUtils
{
  /// Exception type: failed trying to adjust universe weights
  class BadUnivWgtError: public std::runtime_error
  {
    public:
      BadUnivWgtError(const std::string & what) :
        std::runtime_error(what) {};
  };

  /// Exception type: failed while trying to load constraint from file
  class ConstraintLoadError: public std::runtime_error
  {
    public:
      ConstraintLoadError(const std::string & what) :
        std::runtime_error(what) {};
  };

  /// Exception type: couldn't access constraint
  class ConstraintAccessError: public std::runtime_error
  {
    public:
      ConstraintAccessError(const std::string & what) :
        std::runtime_error(what) {};
  };

  /// Exception type: flux universe needed for constraint is missing
  class MissingFluxUnivError: public std::runtime_error
  {
    public:
      MissingFluxUnivError(const std::string & what) :
        std::runtime_error(what) {};
  };

  /// Exception type: strategy for dealing with errors for spectator error band was not specified
  class MissingSpectatorStrategyError: public std::runtime_error
  {
    public:
      MissingSpectatorStrategyError(const std::string & what) :
        std::runtime_error(what) {};
  };

  /// Exception type: no flux universes for constraint
  class NoFluxUnivError: public std::runtime_error
  {
    public:
      NoFluxUnivError(const std::string & what) :
        std::runtime_error(what) {};
  };


  /// Exception type: can't do spread errors with weighted universes
  class NoWgtdSpreadError: public std::runtime_error
  {
    public:
      NoWgtdSpreadError(const std::string & what) :
        std::runtime_error(what) {};
  };

  /// Exception type: an error band cannot both be a constraint and a spectator
  class SpectatorConstraintCollisionError: public std::runtime_error
  {
    public:
      SpectatorConstraintCollisionError(const std::string & what) :
        std::runtime_error(what) {};
  };

  /// Exception type: a vert/lat error band's CV histo is empty.
  class ErrBandEmptyCVError: public std::runtime_error
  {
    public:
      ErrBandEmptyCVError(const std::string & what) :
        std::runtime_error(what) {};
  };



}

#endif /* EXCEPTIONS_H_ */
