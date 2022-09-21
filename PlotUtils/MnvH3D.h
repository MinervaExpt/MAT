#ifndef MNV_MnvH3D_H
#define MNV_MnvH3D_H 1

#include "TObject.h"
#include "TString.h"
#include "TH3D.h"
#include "TVectorD.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvVertErrorBand3D.h"
#include "PlotUtils/MnvLatErrorBand3D.h"
#include <string>
#include <vector>
#include <map>


namespace PlotUtils
{

  class MnvH3D: public TH3D
  {
    public:
      //! Default constructor
      MnvH3D();

      //! Default Constructor with specified normalized bin width
      explicit MnvH3D( Double_t normBinWidthX, Double_t normBinWidthY, Double_t normBinWidthZ);

      //==== Copy Constructors ====//
      //! Copy constructor and Sumw2
      MnvH3D( const TH3D& h2d );

      //==== Copy Constructors with specified normalized bin width ====//
      //! Copy constructor and Sumw2
      MnvH3D( const TH3D& h2d, Double_t fnormBinWidthX, Double_t fnormBinWidthY, Double_t fnormBinWidthZ );

      //==== Constructors with default normalized bin width ====//
      //! Construct with variable bin sizes and Sumw2
      MnvH3D( const char* name, const char* title, Int_t nbinsx, const Float_t* xbins, Int_t nbinsy, const Float_t* ybins, Int_t nbinsz, const Float_t* zbins);

      //! Construct with variable bin sizes and Sumw2
      MnvH3D( const char* name, const char* title, Int_t nbinsx, const Double_t* xbins, Int_t nbinsy, const Double_t* ybins, Int_t nbinsz, const Double_t* zbins);

      //! Construct with constant bin sizes and Sumw2
      MnvH3D( const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, Int_t nbinsz, Double_t zlow, Double_t zup);

      //! Deep copy constructor (the real copy)
      MnvH3D( const MnvH3D& h );

      //! Deep assignment
      MnvH3D& operator=( const MnvH3D& h );

    private:
      //! A helper function which set variables for the deep copy and assignment
      void DeepCopy( const MnvH3D& h );

    public:
      //=== Constructors which specify a default normalization bin width ===//
      //! Construct with variable bin sizes and use this binWidth to normalize bins by default
      MnvH3D( const char* name, const char* title, Int_t nbinsx, const Float_t* xbins, Int_t nbinsy, const Float_t* ybins, Int_t nbinsz, const Float_t* zbins, Double_t normBinWidthX, Double_t normBinWidthY, Double_t normBinWidthZ );

      //! Construct with variable bin sizes and use this binWidth to normalize bins by default
      MnvH3D( const char* name, const char* title, Int_t nbinsx, const Double_t* xbins, Int_t nbinsy, const Double_t* ybins, Int_t nbinsz, const Double_t* zbins, Double_t normBinWidthX, Double_t normBinWidthY, Double_t normBinWidthZ);

      //! Construct with variable bin sizes and use this binWidth to normalize bins by default
      MnvH3D( const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, Int_t nbinsz, Double_t zlow, Double_t zup, Double_t normBinWidthX, Double_t normBinWidthY, Double_t normBinWidthZ);
      //=== end of constructors ===/

      /*! Get an MnvH3D which has its bin content and errors normalized to bin width so it looks smooth
      @param[in] normBinWidthX,normBinWidthY bin width for X and Y axis to normalize to.  normalization = normBinWidthX * normBinWidthY / thisBinArea.
      where thisBinArea = thisbinWidthX*thisbinWidthY (nonpositive means use MnvH3D's default if set)
      @return A copy of this MnvH3D which has its bin content/error normalized to bin width
      */
      //A nonsense function so I can get TransWarpExtractor to compile without figuring out how this was working with the CVS version
      //virtual Int_t GetBin(Int_t binx, Int_t biny){std::cout<<"ERROR: TRYING TO GET 3D BIN USING ONLY X ANY Y. RETURNING 0;"<<std::endl; return 0;  }
      MnvH3D GetBinNormalizedCopy( Double_t normBinWidthX = -1., Double_t normBinWidthY = -1., Double_t normBinWidthZ = -1.) const;
      
      MnvH1D *ProjectionX(const char* name = "_px", Int_t firstybin = 0, Int_t lastybin = -1, Int_t firstzbin = 0, Int_t lastzbin = -1, Option_t* option = "") const;
          
      MnvH1D *ProjectionY(const char* name = "_py", Int_t firstxbin = 0, Int_t lastxbin = -1, Int_t firstzbin = 0, Int_t lastzbin = -1, Option_t* option = "") const;
      
      MnvH1D *ProjectionZ(const char* name = "_pz", Int_t firstxbin = 0, Int_t lastxbin = -1, Int_t firstybin = 0, Int_t lastybin = -1, Option_t* option = "") const;
      
      TH1 *Project3D(Option_t* option = "x") const;

      bool HasVertErrorBand( const std::string& name ) const;
      bool HasLatErrorBand( const std::string& name ) const;

      //! Check for the existence of an error band (of any type)
      bool HasErrorBand( const std::string& name ) const;
      
      bool HasErrorMatrix( const std::string& name ) const;
      
      bool AddVertErrorBand( const std::string& name, const int nhists = -1 );
      //! Add a customed MnvVertErrorBand
      bool AddVertErrorBand( const std::string& name, const std::vector<TH3D*>& base );
      //! Add a new MnvVertErrorBand and fill its universes with the CV 
      bool AddVertErrorBandAndFillWithCV( const std::string& name, const unsigned int nhists );

      bool AddLatErrorBand( const std::string& name, const int nhists = -1 );
      //! Add a customed MnvLatErrorBand
      bool AddLatErrorBand( const std::string& name, const std::vector<TH3D*>& base );

      //! Add a new MnvLattErrorBand and fill its universes with the CV 
      bool AddLatErrorBandAndFillWithCV( const std::string& name, const unsigned int nhists );

      //! Add missing MnvVertErrorBands or MnvLatErrorBands and fill their universes with the CV 
      bool AddMissingErrorBandsAndFillWithCV( const MnvH3D& ref );
      
      //! Delete all Error Bands 
      void ClearAllErrorBands();
      
      //! Fill the weights of an MnvVertErrorBand's universes from a vector
      bool FillVertErrorBand( const std::string& name, const double xval, const double yval, const double zval, const std::vector<double>& weights, const double cvweight  = 1.0, double cvWeightFromMe = 1. );
      //! Fill the weights of a MnvVertErrorBand's universes from array
      bool FillVertErrorBand( const std::string& name, const double xval, const double yval, const double zval, const double * weights, const double cvweight  = 1.0, double cvWeightFromMe = 1. );
      //! Fill the weights of an MnvVertErrorBand's 2 universes with these 2 weights (must have 2)
      bool FillVertErrorBand( const std::string& name, const double xval, const double yval, const double zval, const double weightDown, const double weightUp, const double cvweight  = 1.0, double cvWeightFromMe = 1. );

      //! Fill the weights of a MnvLatErrorBand's universes from a vector
      bool FillLatErrorBand( const std::string& name, const double xval, const double yval, const double zval, const std::vector<double>& xshifts, const std::vector<double>& yshifts, const std::vector<double>& zshifts, const double cvweight  = 1.0, const bool fillcv = true, const double* weights = NULL );
      //! Fill the weights of an MnvLatErrorBand's universes from array
      bool FillLatErrorBand( const std::string& name, const double xval, const double yval, const double zval, const double *xshifts, const double *yshifts, const double *zshifts, const double cvweight  = 1.0, const bool fillcv = true, const double* weights = NULL );
      //! Fill the weights of an MnvLatErrorBand's 2 universes with these 2 weights (must have 2)
      bool FillLatErrorBand( const std::string& name, const double xval, const double yval, const double zval, const double xshiftDown, const double xshiftUp, const double yshiftDown, const double yshiftUp, const double zshiftDown, const double zshiftUp, const double cvweight = 1.0, const bool fillcv = true );

      //! Get a pointer to this MnvVertErrorBand
      MnvVertErrorBand3D* GetVertErrorBand( const std::string& name );
      //! Get a const pointer to this MnvVertErrorBand
      const MnvVertErrorBand3D* GetVertErrorBand( const std::string& name ) const;

      //! Get a pointer to this MnvLatErrorBand
      MnvLatErrorBand3D* GetLatErrorBand( const std::string& name );
      //! Get a const pointer to this MnvLatErrorBand
      const MnvLatErrorBand3D* GetLatErrorBand( const std::string& name ) const;

      /*! Get a new TH3D which is the central value histogram with statistical error only
      @return copy of a TH3D
       */
      TH3D GetCVHistoWithStatError() const;

      /*! Get a new TH3D which is the central value histogram with errors added in quadrature
      @todo add param errNames to add only certain errors from these error bands.
      @return copy of a TH3D
       */
      TH3D GetCVHistoWithError( bool includeStat = true, bool cov_area_normalize = false ) const;
      //    TH1D GetCVHistoWithError( const std::vector<std::string>& errNames ) const;

      //! Get a TH3D filled with the errors summed in quadrature
      TH3D GetTotalError( bool includeStat = true, bool asFrac = false , bool cov_area_normalize = false ) const;

      //! Get a TH3D filled with the stat error ONLY
      TH3D GetStatError( bool asFrac = false ) const;

      //! Get the Total Covariance Matrix
      TMatrixD GetTotalErrorMatrix(bool includeStat = true, bool asFrac = false, bool cov_area_normalize = false ) const;
      //! Get Total Correlation Matrix
      TMatrixD GetTotalCorrelationMatrix(bool cov_area_normalize = false ) const;

      //! Get Statistical Error Matrix
      TMatrixD GetStatErrorMatrix( bool asFrac = false ) const;

      //! Get a single Systematical Matrix from the map
      TMatrixD GetSysErrorMatrix(const std::string& name, bool asFrac = false, bool cov_area_normalize = false) const; 
      //! Get a vector of the names of all error bands
      std::vector<std::string> GetVertErrorBandNames() const;

      //! Get a vector of the names of all error bands
      std::vector<std::string> GetLatErrorBandNames() const;

      //! Get a vector of the names of Systematic Error Matrices
      std::vector<std::string> GetSysErrorMatricesNames() const;

      //! Delete all Systematic Error matrices
      void ClearSysErrorMatrices( );
      //=======================================================================
      // Implementations of ROOT virtual functions
      //=======================================================================
      //! When calling scale we usually need to scale all of the error bands' histos
      virtual void Scale( Double_t c1 = 1., Option_t *option = "", Bool_t allUniv = true );

      /*! Replace this MnvH3D's contents with the result of a division of two other MnvH3Ds
       */
      virtual void Divide( const MnvH3D* h1, const MnvH3D* h2, Double_t c1 = 1, Double_t c2 = 1, Option_t* option="" );

      /*! Replace this MnvH3D's contents with the result of dividing an MnvH3D and all its errors by a TH3
       */
      virtual void DivideSingle( const MnvH3D* h1, const TH3* h2, Double_t c1 = 1, Double_t c2 = 1, Option_t* option="" );
      
      //! Add Histograms
      virtual void Add( const TH3* h1, const Double_t c1 = 1. );
      
      //! Replace this MnvH3D's contents with the result of a multiplication of two other MnvH3Ds
      virtual void Multiply( const MnvH3D* h1, const MnvH3D* h2, const Double_t c1 = 1., const Double_t c2 = 1. );

      //! Get the default bin width to which bin content/error is normalized
      inline Double_t GetNormBinWidthX() const { return fNormBinWidthX; };
      inline Double_t GetNormBinWidthY() const { return fNormBinWidthY; };
      inline Double_t GetNormBinWidthZ() const { return fNormBinWidthZ; };

      //! Set the default bin width to which bin content/error is normalized
      inline void     SetNormBinWidthX(const Double_t x ) { fNormBinWidthX = x; };
      inline void     SetNormBinWidthY(const Double_t y ) { fNormBinWidthY = y; };
      inline void     SetNormBinWidthZ(const Double_t z ) { fNormBinWidthZ = z; };

      //!Destructor (note: root cleans up histograms for us)
      virtual ~MnvH3D() {};

    private:
      //! A helper function to check if this string has that ending
      bool HasEnding (std::string const &fullString, std::string const &ending) const;

      //! Stores a map from name of Systematics Error Matrices 
      std::map<std::string, TMatrixD*> fSysErrorMatrix;

      //! Stores a map from name of removed Systematic Error Matrices
      std::map<std::string, TMatrixD*> fRemovedSysErrorMatrix;

      //! Strores a map from name to error band for MnvLatErrorBands
      std::map<std::string, MnvLatErrorBand3D*> fLatErrorBandMap;

      //! Strores a map from name to error band for MnvVertErrorBands
      std::map<std::string, MnvVertErrorBand3D*> fVertErrorBandMap;

      //! Stores the width to which we will normalize bins (e.g. n Events per fNormBinWidth GeV)
      //! If negative, then refuse to normalize to bin width (appropriate for ratios, efficiencies)
      Double_t fNormBinWidthX;
      Double_t fNormBinWidthY;
      Double_t fNormBinWidthZ;

      //!define a class named MnvH3D, at version 1
      ClassDef(MnvH3D, 1); //MINERvA 2-D histogram

  }; //end of MnvH3D

}//end PlotUtils

#endif
