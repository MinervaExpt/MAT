#ifndef MNV_MnvH2D_H
#define MNV_MnvH2D_H 1


#include "TObject.h"
#include "TString.h"
#include "TH2D.h"
#include "TVectorD.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvVertErrorBand2D.h"
#include "PlotUtils/MnvLatErrorBand2D.h"
#include <string>
#include <vector>
#include <map>

 // let this loose to debug problems with operations

namespace PlotUtils
{
  class MnvH2D: public TH2D
  {
    public:
      //! Default constructor
      MnvH2D();

      //! Default Constructor with specified normalized bin width
      explicit MnvH2D( Double_t normBinWidthX, Double_t normBinWidthY);

      //==== Copy Constructors ====//
      //! Copy constructor and Sumw2
      MnvH2D( const TH2D& h2d );

      //==== Copy Constructors with specified normalized bin width ====//
      //! Copy constructor and Sumw2
      MnvH2D( const TH2D& h2d, Double_t fnormBinWidthX, Double_t fnormBinWidthY );

      //==== Constructors with default normalized bin width ====//
      //! Construct with variable bin sizes and Sumw2
      MnvH2D( const char* name, const char* title, Int_t nbinsx, const Float_t* xbins, Int_t nbinxy, const Float_t* ybins);

      //! Construct with variable bin sizes and Sumw2
      MnvH2D( const char* name, const char* title, Int_t nbinsx, const Double_t* xbins, Int_t nbinsy, const Double_t* ybins);

      //! Construct with constant bin sizes and Sumw2
      MnvH2D( const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup);

      //! Deep copy constructor (the real copy)
      MnvH2D( const MnvH2D& h );

      //! Deep assignment
      MnvH2D& operator=( const MnvH2D& h );
  
    MnvH2D* Clone(const char* name = "" )const;

    protected:
      //! A helper function which set variables for the deep copy and assignment
      void DeepCopy( const MnvH2D& h );

    public:
      //=== Constructors which specify a default normalization bin width ===//
      //! Construct with variable bin sizes and use this binWidth to normalize bins by default
      MnvH2D( const char* name, const char* title, Int_t nbinsx, const Float_t* xbins, Int_t nbinsy, const Float_t* ybins, Double_t normBinWidthX, Double_t normBinWidthY );

      //! Construct with variable bin sizes and use this binWidth to normalize bins by default
      MnvH2D( const char* name, const char* title, Int_t nbinsx, const Double_t* xbins, Int_t nbinsy, const Double_t* ybins, Double_t normBinWidthX, Double_t normBinWidthY);

      //! Construct with variable bin sizes and use this binWidth to normalize bins by default
      MnvH2D( const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, Double_t normBinWidthX, Double_t normBinWidthY);
      //=== end of constructors ===/

      /*! Get an MnvH2D which has its bin content and errors normalized to bin width so it looks smooth
      @param[in] normBinWidthX,normBinWidthY bin width for X and Y axis to normalize to.  normalization = normBinWidthX * normBinWidthY / thisBinArea.
      where thisBinArea = thisbinWidthX*thisbinWidthY (nonpositive means use MnvH2D's default if set)
      @return A copy of this MnvH2D which has its bin content/error normalized to bin width
      */
      MnvH2D GetBinNormalizedCopy( Double_t normBinWidthX = -1., Double_t normBinWidthY = -1.) const;
      
      MnvH2D GetAreaNormalizedCopy() const;

      MnvH1D *Projection(const char* name, bool projectionX, Int_t firstbin, Int_t lastbin, Option_t* option) const;

      MnvH1D *ProjectionX(const char* name = "_px", Int_t firstybin = 0, Int_t lastybin = -1, Option_t* option = "") const;
          
      MnvH1D *ProjectionY(const char* name = "_py", Int_t firstxbin = 0, Int_t lastxbin = -1, Option_t* option = "") const;

      bool HasVertErrorBand( const std::string& name ) const;
      bool HasLatErrorBand( const std::string& name ) const;

      //! Check for the existence of an error band (of any type)
      bool HasErrorBand( const std::string& name ) const;
      
      bool HasErrorMatrix( const std::string& name ) const;
  
      bool AddVertErrorBand( const std::string& name, const int nhists = -1 );
      //! Add a customed MnvVertErrorBand
      bool AddVertErrorBand( const std::string& name, const std::vector<TH2D*>& base );
      //! Add a new MnvVertErrorBand and fill its universes with the CV 
      bool AddVertErrorBandAndFillWithCV( const std::string& name, const unsigned int nhists );

      bool AddLatErrorBand( const std::string& name, const int nhists = -1 );
      //! Add a customed MnvLatErrorBand
      bool AddLatErrorBand( const std::string& name, const std::vector<TH2D*>& base );
      //! Add a new MnvLattErrorBand and fill its universes with the CV 
      bool AddLatErrorBandAndFillWithCV( const std::string& name, const unsigned int nhists );

      //! Add missing MnvVertErrorBands or MnvLatErrorBands and fill their universes with the CV (templated)
      // templated because then it can handle MnvH1Ds too.  implementation is below, in THIS file.
      template <typename MnvHType>
      bool AddMissingErrorBandsAndFillWithCV( const MnvHType& ref );

      //! MnvH1D version
      bool AddMissingErrorBandsAndFillWithCV( const MnvH1D& ref )  { return AddMissingErrorBandsAndFillWithCV<MnvH1D>(ref); };
      
      //! MnvH2D version
      bool AddMissingErrorBandsAndFillWithCV( const MnvH2D& ref )  { return AddMissingErrorBandsAndFillWithCV<MnvH2D>(ref); };

      //! Rename histos and error bands obviously     
      void RenameHistosAndErrorBands( const std::string& newname);
 
      //! Delete all Error Bands 
      void ClearAllErrorBands();

      //! Transfer all lat and vert error bands to hist
      bool TransferErrorBands( MnvH2D *hist, bool removeFromMe = true );
      //! Transfer a lat error band to hist
      bool TransferVertErrorBand( MnvH2D *hist, const std::string& errName, bool removeFromMe = true );
      //! Transfer a lat error band to hist
      bool TransferLatErrorBand( MnvH2D *hist, const std::string& errName, bool removeFromMe = true );
  
  //! Transfer a lat error band to hist
      bool TransferSysErrorMatrix( MnvH2D *hist, const std::string& errName, bool removeFromMe = true );

      //! Remove the error band from the MnvH1D and return the pointer to the object (caller owns the memory)
      MnvVertErrorBand2D* PopVertErrorBand( const std::string& name );
      //! Remove the error band from the MnvH1D and return the pointer to the object (caller owns the memory)
      MnvLatErrorBand2D* PopLatErrorBand( const std::string& name );
      //! Remove the cov matrix band from the MnvH1D and return the pointer to the object (caller owns the memory)
      TMatrixD* PopSysErrorMatrix( const std::string& name );

      //! Add this vert error band (with memory ownership) to the MnvH1D 
      bool PushErrorBand( const std::string& name, MnvVertErrorBand2D* errBand );
      //! Add this lat error band (with memory ownership) to the MnvH1D 
      bool PushErrorBand( const std::string& name, MnvLatErrorBand2D* errBand );
      //! Add this lat error band (with memory ownership) to the MnvH1D
      bool PushSysErrorMatrix( const std::string& name, TMatrixD* errBand );


      //! Fill the weights of an MnvVertErrorBand's universes from a vector
      bool FillVertErrorBand( const std::string& name, const double xval, const double yval, const std::vector<double>& weights, const double cvweight  = 1.0, double cvWeightFromMe = 1. );
      //! Fill the weights of a MnvVertErrorBand's universes from array
      bool FillVertErrorBand( const std::string& name, const double xval, const double yval, const double * weights, const double cvweight  = 1.0, double cvWeightFromMe = 1. );
      //! Fill the weights of an MnvVertErrorBand's 2 universes with these 2 weights (must have 2)
      bool FillVertErrorBand( const std::string& name, const double xval, const double yval, const double weightDown, const double weightUp, const double cvweight  = 1.0, double cvWeightFromMe = 1. );

      //! Fill the weights of a MnvLatErrorBand's universes from a vector
      bool FillLatErrorBand( const std::string& name, const double xval, const double yval, const std::vector<double>& xshifts, const std::vector<double>& yshifts, const double cvweight  = 1.0, const bool fillcv = true, const double* weights = 0 );
      //! Fill the weights of an MnvLatErrorBand's universes from array
      bool FillLatErrorBand( const std::string& name, const double xval, const double yval, const double *xshifts, const double *yshifts, const double cvweight  = 1.0, const bool fillcv = true, const double* weights = 0 );
      //! Fill the weights of an MnvLatErrorBand's 2 universes with these 2 weights (must have 2)
      bool FillLatErrorBand( const std::string& name, const double xval, const double yval, const double xshiftDown, const double xshiftUp, const double yshiftDown, const double yshiftUp, const double cvweight = 1.0, const bool fillcv = true );

      //! Get a pointer to this MnvVertErrorBand
      MnvVertErrorBand2D* GetVertErrorBand( const std::string& name );
      //! Get a const pointer to this MnvVertErrorBand
      const MnvVertErrorBand2D* GetVertErrorBand( const std::string& name ) const;

      //! Get a pointer to this MnvLatErrorBand
      MnvLatErrorBand2D* GetLatErrorBand( const std::string& name );
      //! Get a const pointer to this MnvLatErrorBand
      const MnvLatErrorBand2D* GetLatErrorBand( const std::string& name ) const;

      /*! Get a new TH2D which is the central value histogram with statistical error only
        @return copy of a TH2D
       */
      TH2D GetCVHistoWithStatError() const;

      /*! Get a new TH2D which is the central value histogram with errors added in quadrature
        @todo add param errNames to add only certain errors from these error bands.
        @return copy of a TH2D
       */
      TH2D GetCVHistoWithError( bool includeStat = true, bool cov_area_normalize = false ) const;
      //    TH1D GetCVHistoWithError( const std::vector<std::string>& errNames ) const;

      //! Get a TH2D filled with the errors summed in quadrature
      TH2D GetTotalError( bool includeStat = true, bool asFrac = false , bool cov_area_normalize = false ) const;

      //! Get a TH2D filled with the stat error ONLY
      TH2D GetStatError( bool asFrac = false ) const;

      bool PushCovMatrix( const std::string& name, TMatrixD covmx,  bool cov_area_normalize = false );

      //! Modify Covariance Matrix and diagonal stat unc. based on studies by DGR related to unfolding and finite MC sizes. See docDB 28992, 28899
      void ModifyStatisticalUnc(double factor, std::string covmatrixname = "unfoldingCov");

      //! Get the Total Covariance Matrix
      TMatrixD GetTotalErrorMatrix(bool includeStat = true, bool asFrac = false, bool cov_area_normalize = false ) const;
      //! Get Total Correlation Matrix
      TMatrixD GetTotalCorrelationMatrix(bool cov_area_normalize = false ) const;

      //! Get Statistical Error Matrix
      TMatrixD GetStatErrorMatrix( bool asFrac = false ) const;

      //! Get a single Systematical Matrix from the map
      TMatrixD GetSysErrorMatrix(const std::string& name, bool asFrac = false, bool cov_area_normalize = false) const;
  
      //! Allow user to replace a sys error matrix after external manipulation
  
      bool FillSysErrorMatrix(const std::string& name, const TMatrixD&  matrix);

      //! How many MnvLatErrorBands are there?
      size_t GetNLatErrorBands() const { return fLatErrorBandMap.size(); };
      //! How many MnvVertErrorBands are there?
      size_t GetNVertErrorBands() const { return fVertErrorBandMap.size(); };
      //! How many SysErrorMatrices are there?
      size_t GetNSysErrorMatrices() const { return fSysErrorMatrix.size(); };
      //! How many SysErrorMatrices are there?
      size_t GetNRemovedSysErrorMatrices() const { return fRemovedSysErrorMatrix.size(); };
      //! How many error sources are there (exlcudes removed SysErrorMatrices)?
      size_t GetNErrorSources() const { return fLatErrorBandMap.size() + fVertErrorBandMap.size() + fSysErrorMatrix.size(); };

      //! Get a vector of the names of all error bands
      std::vector<std::string> GetVertErrorBandNames() const;

      //! Get a vector of the names of all error bands
      std::vector<std::string> GetLatErrorBandNames() const;

      //! Get a vector of the names of all error bands
      std::vector<std::string> GetErrorBandNames() const;

      std::vector<std::string> GetUncorrErrorNames() const { return std::vector<std::string>(); };

      //! Get a vector of the names of Systematic Error Matrices
      std::vector<std::string> GetSysErrorMatricesNames() const;
  
      //! Get a vector of the names of standalone Systematic Error Matrices
      std::vector<std::string> GetCovMatricesNames() const;


      //! Delete all Systematic Error matrices
      void ClearSysErrorMatrices( );
      //=======================================================================
      // Implementations of ROOT virtual functions
      //=======================================================================
      //! When calling scale we usually need to scale all of the error bands' histos
      virtual void Scale( Double_t c1 = 1., Option_t *option = "", Bool_t allUniv = true );

      /*! Replace this MnvH2D's contents with the result of a division of two other MnvH2Ds
       */
      virtual void Divide( const MnvH2D* h1, const MnvH2D* h2, Double_t c1 = 1, Double_t c2 = 1, Option_t* option="" );

      /*! Replace this MnvH2D's contents with the result of dividing an MnvH2D and all its errors by a TH2
       */
      virtual void DivideSingle( const MnvH2D* h1, const TH2* h2, Double_t c1 = 1, Double_t c2 = 1, Option_t* option="" );

      //! Add Histograms
      virtual void Add( const TH2* h1, const Double_t c1 = 1. );

      //! Replace this MnvH2D's contents with the result of a multiplication of two other MnvH2Ds
      virtual void Multiply( const MnvH2D* h1, const MnvH2D* h2, const Double_t c1 = 1., const Double_t c2 = 1. );

      //! Replace this MnvH2D's contents with the results of multiplying an MnvH2D and all its errors by a TH2
      virtual void MultiplySingle( const MnvH2D* h1, const TH2* h2, const Double_t c1 = 1., const Double_t c2 = 1. );

      //! Reset (empty) the CV and all error histograms
      void Reset( Option_t *option = "" );

      //! Get the default bin width to which bin content/error is normalized
      inline Double_t GetNormBinWidthX() const { return fNormBinWidthX; };
      inline Double_t GetNormBinWidthY() const { return fNormBinWidthY; };

      //! Set the default bin width to which bin content/error is normalized
      inline void     SetNormBinWidthX(const Double_t x ) { fNormBinWidthX = x; };
      inline void     SetNormBinWidthY(const Double_t y ) { fNormBinWidthY = y; };

      //! Does sumw2 and prints no warning
      void SilentSumw2();

      //!Destructor
      virtual ~MnvH2D();
  
    //! Set Strict flag
      inline void SetStrict(Bool_t strict=false){
        fStrict = strict;
      };
  
    //! Get Strict flag
      inline Bool_t GetStrict(){
        return fStrict;
      }
  //! dump a histogram not a method a function
      void MnvH2DToCSV(std::string name, std::string directory="", double scale=1.0, bool fullprecision=true, bool syserrors=true, bool percentage = true, bool binwidth =true);

    protected:
      //! A helper function to check if this string has that ending
      bool HasEnding (std::string const &fullString, std::string const &ending) const;

      //! Stores a map from name of Systematics Error Matrices 
      std::map<std::string, TMatrixD*> fSysErrorMatrix;

      //! Stores a map from name of removed Systematic Error Matrices
      std::map<std::string, TMatrixD*> fRemovedSysErrorMatrix;

      //! Strores a map from name to error band for MnvLatErrorBands
      std::map<std::string, MnvLatErrorBand2D*> fLatErrorBandMap;

      //! Strores a map from name to error band for MnvVertErrorBands
      std::map<std::string, MnvVertErrorBand2D*> fVertErrorBandMap;

      //! Stores the width to which we will normalize bins (e.g. n Events per fNormBinWidth GeV)
      //! If negative, then refuse to normalize to bin width (appropriate for ratios, efficiencies)
      Double_t fNormBinWidthX;
      Double_t fNormBinWidthY;
      Double_t fStrict;

      //!define a class named MnvH2D, at version 1
  // put in some fixes for covariance matrices
      ClassDef(MnvH2D, 2); //MINERvA 2-D histogram

  }; //end of MnvH2D

  // templated methods have to be in the .h file
  template <typename MnvHType>
  bool MnvH2D::AddMissingErrorBandsAndFillWithCV( const MnvHType& ref )
  {
    // Declare container for error band names
    std::vector<std::string> names = ref.GetVertErrorBandNames();

    // Add vertical error bands found in reference MnvH1D
    for( std::vector<std::string>::iterator name=names.begin(); name!=names.end(); ++name )
    {
      // Skip already added vertical error bands
      if( HasVertErrorBand( *name ) ) 
        continue;
      unsigned int nunis = ref.GetVertErrorBand( *name )->GetNHists();
      if( ! this->AddVertErrorBandAndFillWithCV( *name, nunis ) )
        return false;
      if ( ref.GetVertErrorBand(*name)->GetUnivWgts() ) GetVertErrorBand(*name )->SetUnivWgts( *(ref.GetVertErrorBand(*name)->GetUnivWgts()) );
    }

    // Add lateral error bands found in reference MnvH1D
    names = ref.GetLatErrorBandNames();
    for( std::vector<std::string>::iterator name=names.begin(); name!=names.end(); ++name )
    { 
      // Skip already added lateral error bands
      if( HasLatErrorBand( *name ) ) 
        continue;
      unsigned int nunis = ref.GetLatErrorBand( *name )->GetNHists();
      if( ! this->AddLatErrorBandAndFillWithCV( *name, nunis ) ) 
        return false;
      if ( ref.GetLatErrorBand(*name)->GetUnivWgts() ) GetLatErrorBand(*name )->SetUnivWgts( *(ref.GetLatErrorBand(*name)->GetUnivWgts()) );
    }
  // used to be GetSysErrorMatricesNames but meant Covmatrix
  names = ref.GetCovMatricesNames();
  for( std::vector<std::string>::iterator name=names.begin(); name!=names.end(); ++name )
  {
    // Skip already added error matrices
    if( HasErrorMatrix( *name ) )

      std::cout << "MnvH2D::AddMissingErrorBandsAndFillWithCV from " << ref.GetName() << " to " << GetName() << " already has " << *name << std::endl;
       continue;
    TMatrixD* newmatrix = new TMatrixD(ref.GetSysErrorMatrix(*name));
    std::cout << "MnvH2D::AddMissingErrorBandsAndFillWithCV from " << ref.GetName() << " to " << GetName() << " just moved " << *name << std::endl;
    if( ! this->PushSysErrorMatrix( *name, newmatrix ) )
      return false;
  }
  
    

    return true;
  }

}//end PlotUtils


#endif
