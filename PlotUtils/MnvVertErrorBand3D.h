#ifndef MNV_MnvVertErrorBand3D_H
#define MNV_MnvVertErrorBand3D_H 1

#include "TObject.h"
#include "TString.h"
#include "TError.h"
#include "TNamed.h"
#include "TH3D.h"
#include "TMatrixD.h"
#include "TMatrixDBase.h"

#include <assert.h>
#include <vector>
#include <string>
#include <iostream>
#include <math.h>

namespace PlotUtils
{

class MnvVertErrorBand3D : public TH3D
{
  public:
    //! Default constructor 
    MnvVertErrorBand3D( ) : TH3D() {};

    MnvVertErrorBand3D( const std::string& name, const TH3D* base, const unsigned int nHists = 1000 );

    //! Add a new Constructor for already Filled vector of Histogram Error Bands
    MnvVertErrorBand3D( const std::string& name, const TH3D* base, const std::vector<TH3D*>& hists );

    //! Deep copy constructor
    MnvVertErrorBand3D( const MnvVertErrorBand3D& h );

    //! Deep assignment operator
    MnvVertErrorBand3D& operator=( const MnvVertErrorBand3D& h );

  private:
    //! A helper function which sets variables for the deep copy and assignment
    void DeepCopy( const MnvVertErrorBand3D& h );

  public:
    //!Destructor (note: root cleans up histograms for us)
    virtual ~MnvVertErrorBand3D() {};

    //! Fill the CV histo and all the universes' histos
    virtual bool Fill( const double xval, const double yval, const double zval, const double *weights, const double cvweight = 1, double cvWeightFromMe = 1.);

    virtual bool Fill( const double xval, const double yval, const double zval, const double weightDown, const double weightUp, const double cvweight = 1.0, double cvWeightFromMe = 1. );

    //! Get the error band histogram
    virtual TH3D GetErrorBand(bool asFrac = false , bool cov_area_normalize = false) const;
    
    //! Will the error band come from max spread?
    bool GetUseSpreadError() const { return fUseSpreadError; };

    //! Set the error band to come from max spread or standard dev
    void SetUseSpreadError( bool use ) { fUseSpreadError = use; };

    //! How many histograms (universes) does this have?
    unsigned int GetNHists() const { return fNHists; };

    //! Get the universes' histograms (const)
    const std::vector<TH3D*>& GetHists() const { return fHists; };

    //! Get a specific universe's histogram (const)
    const TH3D* GetHist(const unsigned int i) const;

    //! Get a specific universe's histogram (nonconst)
    TH3D* GetHist(const unsigned int i);

    //! Get the universes' histograms (nonconst)
    std::vector<TH3D*> GetHists() { return fHists; };

    //! Calculate Covariance Matrix
    TMatrixD CalcCovMx(bool area_normalize = false, bool asFrac = false) const;

    //! Add h1*c1 to this error band
    Bool_t Add( const MnvVertErrorBand3D* h1, const Double_t c1 = 1. );

    //! Replace this MnvVertErrorBand3D's contents with the result of a multiplication of two other MnvVertErrorBands3D
    Bool_t Multiply( const MnvVertErrorBand3D* h1, const MnvVertErrorBand3D* h2, Double_t c1 = 1, Double_t c2 = 1 );

    /*! Replace this MnvVertErrorBand3D's contents with the result of a division of one MnvVertErrorBand3D by a TH2
     */
    Bool_t DivideSingle( const MnvVertErrorBand3D* h1, const TH3* h2, Double_t c1 = 1, Double_t c2 = 1, Option_t* option="" );

    /*! Replace this MnvVertErrorBand3D's contents with the result of a division of two other MnvVertErrorBand3Ds
     */
    Bool_t Divide( const MnvVertErrorBand3D* h1, const MnvVertErrorBand3D* h2, Double_t c1 = 1, Double_t c2 = 1, Option_t* option="" );
    
    //! Scale all universes by some constant using TH3D::Scale
    void Scale( Double_t c1 = 1., Option_t *option = "", Bool_t allUniv = true );


  protected:
    unsigned int fNHists;         ///< Number of histograms (universes)
    bool fUseSpreadError;         ///< Are we using spread in histos to get the error
    std::vector<TH3D*> fHists;    ///< Vector of histograms for the universes
    std::vector<int> fGoodColors; ///< Current list of good colors to use for universe histos

  private:
    //!define a class named MnvVertErrorBand3D, at version 1 
    ClassDef( MnvVertErrorBand3D, 1 );
}; //end of MnvVertErrorBand3D

} //end of PlotUtils

#endif
