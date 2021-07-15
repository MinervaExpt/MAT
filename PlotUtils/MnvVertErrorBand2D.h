#ifndef MNV_MnvVertErrorBand2D_H
#define MNV_MnvVertErrorBand2D_H 1

#include "TObject.h"
#include "TString.h"
#include "TError.h"
#include "TNamed.h"
#include "TH2D.h"
#include "TMatrixD.h"
#include "TMatrixDBase.h"

#include <assert.h>
#include <vector>
#include <string>
#include <iostream>
#include <math.h>

namespace MAT
{

class MnvVertErrorBand2D : public TH2D
{
  public:
    //! Default constructor 
    MnvVertErrorBand2D( ) : TH2D() {};

    MnvVertErrorBand2D( const std::string& name, const TH2D* base, const unsigned int nHists = 1000 );

    //! Add a new Constructor for already Filled vector of Histogram Error Bands
    MnvVertErrorBand2D( const std::string& name, const TH2D* base, const std::vector<TH2D*>& hists );

    //! Deep copy constructor
    MnvVertErrorBand2D( const MnvVertErrorBand2D& h );

    //! Deep assignment operator
    MnvVertErrorBand2D& operator=( const MnvVertErrorBand2D& h );

  private:
    //! A helper function which sets variables for the deep copy and assignment
    void DeepCopy( const MnvVertErrorBand2D& h );

  public:
    //!Destructor
    virtual ~MnvVertErrorBand2D();

    //! Fill the CV histo and all the universes' histos
    virtual bool Fill( const double xval, const double yval, const double *weights, const double cvweight = 1, double cvWeightFromMe = 1. );

    virtual bool Fill( const double xval, const double yval, const double weightDown, const double weightUp, const double cvweight = 1.0, double cvWeightFromMe = 1. );

    //! Get the error band histogram
    virtual TH2D GetErrorBand(bool asFrac = false , bool cov_area_normalize = false) const;
    
    //! Will the error band come from max spread?
    bool GetUseSpreadError() const { return fUseSpreadError; };

    //! Set the error band to come from max spread or standard dev
    void SetUseSpreadError( bool use ) { fUseSpreadError = use; };

    //! How many histograms (universes) does this have?
    unsigned int GetNHists() const { return fNHists; };

    //! Get the universes' histograms (const)
    const std::vector<TH2D*>& GetHists() const { return fHists; };

    //! Get a specific universe's histogram (const)
    const TH2D* GetHist(const unsigned int i) const;

    //! Get a specific universe's histogram (nonconst)
    TH2D* GetHist(const unsigned int i);

    //! Retrieve the weight for a single universe.  If no universes have weights present, will return -1.
    double GetUnivWgt(unsigned int uid) const;

    //! Retrieve all the weights together.  If no universes have weights present, will return NULL.
    const std::vector<double> * GetUnivWgts() const;

    //! Set the weight for a single universe.  (If other universes' weights are as-yet unspecified, they are set to 1.)
    void SetUnivWgt(unsigned int uid, double wgt);

    //! Set all universe weights simultaneously.  Must set exactly fNHists universes (or an exception will be thrown).
    void SetUnivWgts(const std::vector<double>& wgts);

    //! Get the universes' histograms (nonconst)
    std::vector<TH2D*> GetHists() { return fHists; };

    //! Calculate Covariance Matrix
    TMatrixD CalcCovMx(bool area_normalize = false, bool asFrac = false) const;

    //! Add h1*c1 to this error band
    Bool_t Add( const MnvVertErrorBand2D* h1, const Double_t c1 = 1. );

    //! Replace this MnvVertErrorBand2D's contents with the result of a multiplication of two other MnvVertErrorBands2D
    Bool_t Multiply( const MnvVertErrorBand2D* h1, const MnvVertErrorBand2D* h2, Double_t c1 = 1, Double_t c2 = 1 );

    //! Replace this MnvVertErrorBand's contents with the result of a multiplication of another MnvLatErrorband and a TH1
    Bool_t MultiplySingle( const MnvVertErrorBand2D* h1, const TH2* h2, Double_t c1 = 1, Double_t c2 = 1 );

    /*! Replace this MnvVertErrorBand2D's contents with the result of a division of one MnvVertErrorBand2D by a TH2
     */
    Bool_t DivideSingle( const MnvVertErrorBand2D* h1, const TH2* h2, Double_t c1 = 1, Double_t c2 = 1, Option_t* option="" );

    /*! Replace this MnvVertErrorBand2D's contents with the result of a division of two other MnvVertErrorBand2Ds
     */
    Bool_t Divide( const MnvVertErrorBand2D* h1, const MnvVertErrorBand2D* h2, Double_t c1 = 1, Double_t c2 = 1, Option_t* option="" );

    //! Reset all histograms known to this error band
    void Reset( Option_t *option = "" );
    
    //! Scale the CV, and optionally all universes, by some constant using TH2D::Scale
    void Scale( Double_t c1 = 1., Option_t *option = "", Bool_t allUniv = true );


  protected:
    unsigned int fNHists;          ///< Number of histograms (universes)
    bool fUseSpreadError;          ///< Are we using spread in histos to get the error
    std::vector<double> fUnivWgts; ///< Weights for the universes.  1 is used for each universe if this goes unfilled.
    std::vector<TH2D*> fHists;     ///< Vector of histograms for the universes
    std::vector<int> fGoodColors;  ///< Current list of good colors to use for universe histos

  private:
    // adds the necessary ROOT garbledygook to make this class persistable
    ClassDef( MnvVertErrorBand2D, 2 );
}; //end of MnvVertErrorBand2D

} //end of PlotUtils

#endif
