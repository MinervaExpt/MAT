#ifndef HISTOGRAMUTILS_H
#define HISTOGRAMUTILS_H

//---------------------------------------------------------
// Very general functions for manipulating root histograms 
//----------------------------------------------------------

#include "MnvH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TMatrix.h"

#include<iomanip>
#include<iostream>
#include<stdlib.h>
#include<fstream>
#include<math.h>
#include<stdio.h>
#include<string>
#include<sstream>

using std::string;
using std::cout;
using std::flush;
using std::setw;
using std::setprecision;
using std::endl;

namespace MAT {



  class MnvHist {
    public:

    //! Do not Fill anything when the shift has this value.
    static const double NotPhysicalShiftNumber;
    //! Return true if the shift agrees with the NotPhysicalShiftNumber to 6 decimals
    static bool IsNotPhysicalShift( double shift );

    //! If histogram has this as max or min, then let the histogram decide what the min/max should be
    static const double AutoAxisLimit;  
    //! Return true if the number agrees with ROOT's use auto axis limit to 6 decimals
    static bool IsAutoAxisLimit( double x );


    //! Get the interquartile range of a vector
    static double GetInterquartileRange( const std::vector<double>& vec );
    //! Convert Interquartile range to sigma
    static const double InterquartileRangeToSigma;

    static TH3D*   divide3D( const TH3D *num, const TH3D *den );
    static TH2D*   divide2D( const TH2D *num, const TH2D *den );
    static MnvH1D* divide1D( const MnvH1D *num, const MnvH1D *den );

    static TH2D*   average3D_to_2D( const TH3D *object, const char *axis1, const char *axis2 );
    static MnvH1D* average2D_to_1D( const TH2D *object, const char *axis1 );

    static MnvH1D* integrate2D_to_1D( const TH2D *object, const char *axis1, const double mult = 1, const char *option = "none" );
    static double  integrate2D( const TH2D* object, const double mult = 1, const char *option = "none" );

    static void printHisto( TH2D *hist, string name = "2D histo" );
    static void printMatrix( TMatrix matrix, string name = "matrix" );

    //! add the contents of b to a in quadrature
    static int AddInQuadrature( TH1* a, const TH1* b);

    //! Put histogram bin errors into a diagonal matrix
    static TMatrixD GetErrorsAsMatrix( const TH1D *h );

  }; //end of MnvHist

}//end of PlotUtils

#endif


