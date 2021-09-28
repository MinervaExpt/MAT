#ifndef MNV_MnvH2DLOG_H
#define MNV_MnvH2DLOG_H 1

#include "TObject.h"
#include "TString.h"
#include "TH2D.h"
#include "TVectorD.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvVertErrorBand2D.h"
#include "PlotUtils/MnvLatErrorBand2D.h"
#include <string>
#include <vector>
#include <map>
#include <cmath>

 using namespace PlotUtils;
 
  /* To use this, you need to instantiate as Clone of an existing MnvH2D
  
  MnvH2DLog * me = (MnvH2DLog) mnvh2d->Clone();
  
  // then you must (and only once) applyt the logify function
  
  me->logify();
   
  Currently not protected against doing this twice as I had trouble figure out how to initialize a variable in a clone.

  */
  
  
  class MnvH2DLog: public MnvH2D
  {
    public:
      //! Default constructor
  
     //! Deep copy constructor (the real copy)
      MnvH2DLog( const MnvH2DLog& h );
 
  
    private:
      //! A helper function which set variables for the deep copy and assignment
      void DeepCopy( const MnvH2DLog& h );
  
  void logifyOne(TH2* h,const double lambda=0.0);  // lambda is the Box-Cox parameter

  inline double BoxCoxVal(const double y, const double lambda){
    if ( lambda == 0.0){
      return log(y);
    }
    else{
      return (std::pow(y,lambda)-1.)/lambda;
    }
  }
  
  inline double BoxCoxDeriv(const double y, const double lambda){
    if ( lambda == 0.0){
      if (y!=0)      return 1/y;
      return 0.0;
    }
    else{
      return std::pow(y,lambda-1.);
    }
  }
  public:
  void logify(const double lambda = 0.0);

 


  }; //end of MnvH2DLog


//}//end PlotUtils


#endif
