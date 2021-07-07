#ifndef MNV_MnvRecoShifter_h 
#define MNV_MnvRecoShifter_h 1


#include <vector>
#include <TRandom.h>
#include <TMath.h>


#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/Rotation3D.h"
#include "Math/AxisAngle.h"

#include <TVector3.h>
#include <TLorentzVector.h>

// Do not use the Gaudi version anywhere
//#ifndef PLOTUTILS_STANDALONE
//#include "Kernel/MinervaPhysicalConstants.h"
//#else
#include "PlotUtilsPhysicalConstants.h"
//#endif

namespace PlotUtils
{

  using namespace ROOT::Math;
  using namespace TMath;

  static const double M_p_GeV = 0.938272013;
  static const double M_n_GeV = 0.93956536;
  /*! A class to return the shift due to an nSigma mismeasurement on a reconstructed quantity

    This is a makeshift solution because we do not know our reconstruction systematics well enough yet to implement a Gaudi tool.

    The functions have the following types:
    <ol>
    <li>
    bool GetShift_<on>( double& shift, [necessary input params], double nSigma = 1. ).
    For example: GetShift_muonE( shift, muonE, muonEerr, -.5 )  get shift on muonE if the measurement were -.5 sigma away from its measured value for an event with muonE+/- muonEerr.  Return value is true if everything went well.
    </li>
    <li>
    bool GetShift_<on>_<from>( double& shift, [necessary input params...], double nSigma = 1. ).
    For example: GetShift_x_muonE(muonE, hadronicE, muonEerr) get shift on x due to +1 sigma shift in the measurement of muonE for an event with muonE, hadronicE and with an uncertainty on muonE of muonEerr. That uncertainty now comes from MuonUtils:calculateMomentumCorrection().  Return value is true if everything went well.
    </li>
    </ol>

    All variables are in MeV, ns, mm.  In the code, I use the following abbreviations:
    <table>
    <tr><th> In Code </th><th> Translation </th></tr>
    <tr><td> E </td><td> Neutrino Energy </td></tr>
    <tr><td> muonE </td><td> Muon Energy </td></tr>
    <tr><td> hadronicE </td><td> Hadronic Energy </td></tr>
    <tr><td> q2 </td><td> 4-momentum squared </td></tr>
    <tr><td> x </td><td> Bjorken X </td></tr>
    <tr><td> y </td><td> Inelasticity </td></tr>
    <tr><td> muonTheta </td><td> Theta of the muon wrt beam </td></tr>
    </table>

    @note The GetShift functions are virtual, so each analysis has the chance to use specialize resolution functions.

    @author Brian Tice
   */
  class MnvRecoShifter
  {
    //! vector of double iterators are used a lot
    typedef std::vector<double>  dVec;
    typedef dVec::const_iterator dVecCItr;
    typedef dVec::iterator       dVecItr;

    public:

      //! Default constructor - default to neutrino.
      MnvRecoShifter();

      //! Constructor with nu/antinu declaration.
      explicit MnvRecoShifter( bool isAnti);

      MnvRecoShifter( bool isAnti, int nUniverses);

      //! Default destructor
      virtual ~MnvRecoShifter();

      virtual void setIsAntiNu( bool isAnti )  { m_isAntiNu = isAnti; };
      virtual bool isAntiNu()            const { return m_isAntiNu; };
      virtual int GetNumberOfUniverses() const { return m_nUniverses; };

      //=================================
      // BEGIN - Calculators
      //=================================
      //! @name Calculators
      //!@{

      /*! @brief Calculate Q^2.
       */
      virtual bool Calc_q2( double& q2, const double muonE, const double hadronicE, const double muonTheta) const;
      /*! @brief Calculate W, under the assumption M = (M_proton + M_neutron)/2.
       */
      virtual bool Calc_w( double& w, const double muonE, const double hadronicE, const double muonTheta) const;
      /*! @brief Calculate Bjorken x.
       */
      virtual bool Calc_x( double& x, const double muonE, const double hadronicE, const double muonTheta ) const;
      /*! @brief Calculate the inelasticity.
       */
      virtual bool Calc_y( double& y, const double muonE, const double hadronicE ) const;

      //!@}
      //=================================
      // END - Calculators
      //=================================

      //===========================================
      // BEGIN - CCQE Calculators and Scale Shifts
      //===========================================
    
      virtual bool Calc_ccqe_q2( double& q2, const double lepP, const double lepTheta, const double bindingE, double lepMass=MinervaUnits::M_mu) const;

      virtual bool Calc_ccqe_nuE( double& nuE, const double muonP, const double muonTheta, const double bindingE, double lepMass=MinervaUnits::M_mu) const;

      //===========================================
      // END - CCQE Calculators and Scale Shifts
      //===========================================

      //===========================================
      // BEGIN - Transverse Variables Calculators and Scale Shifts
      //===========================================

      virtual bool Calc_tki_vars( dVec& vars, const TVector3 &nu3P,const TLorentzVector &lepton, const TLorentzVector &nucleon, const double NucleusMass/*GeV*/, const double ISNucleonMass = M_n_GeV, const double bindingE = 0/*GeV*/) const;
      virtual bool Calc_tki_vars( dVec& vars, const XYZVector &nu3P,const XYZTVector &lepton, const XYZTVector &nucleon, const double NucleusMass/*GeV*/, const double ISNucleonMass = M_n_GeV, const double bindingE = 0/*GeV*/) const;

      virtual bool Calc_expected_nucleon( TLorentzVector& expected, TVector3 &nu3P, TLorentzVector &lepton, double ISMass= M_n_GeV, double FSMass= M_p_GeV, double BindingE=0 ) ;
      virtual bool Calc_expected_nucleon( XYZTVector& expected, XYZVector &nu3P, XYZTVector &lepton, double ISMass=M_n_GeV, double FSMass=M_p_GeV, double BindingE=0 ) ;

      virtual bool ComputeNeutronAngularVars( dVec & res, TVector3 &nu, TLorentzVector& lepton, TVector3 &targetVec,double ISMass=M_n_GeV, double FSMass=M_p_GeV, double BindingE=0  );
      virtual bool ComputeNeutronAngularVars( dVec & res, XYZVector &nu, XYZTVector& lepton, XYZVector &targetVec,double ISMass=M_n_GeV, double FSMass=M_p_GeV, double BindingE=0  );


      //virtual bool GetShifts_beam( TVector3& shift, const double error, const int nUniverse=-1, const double nSigma=1 ) const;
      virtual bool GetShifts_beamTheta( TVector3& shifted_beam, const double error, const std::string axis = "x", const int nUniverse=-1, const double nSigma=1 ) const;
      virtual bool GetShifts_muon( TVector3& shift, const double error, const int nUniverse=-1, const double nSigma=1 ) const;
      virtual bool GetShifts_hadron( TVector3& shift, const double error, const int nUniverse=-1, const double nSigma=1 ) const;
      //===========================================
      // END - Transverse Variables Calculators and Scale Shifts
      //===========================================


      //-----------------------------------
      // Muon momentum shifts 
      //-----------------------------------
      virtual bool GetShift_ccqe_muonP( double& shift, const double muonPErr, const int nUniverse = -1, const double nSigma = 1.) const;

      virtual bool GetShift_ccqe_muonE_muonP(double& shift, const double muonP, const double muonPErr, const int nUniverse = -1, const double nSigma = 1.) const;

      virtual bool GetShift_ccqe_muonT_muonP(double& shift, const double muonP, const double muonPErr, const int nUniverse = -1, const double nSigma = 1.) const;

      virtual bool GetShift_ccqe_nuE_muonP(double &shift, const double muonP, const double muonTheta, const double bindingE, const double muonPErr, const int nUniverse = -1, const double nSigma = 1.) const;

      virtual bool GetShift_ccqe_q2_muonP(double &shift, const double muonP, const double muonTheta, const double bindingE, const double muonPErr, const int nUniverse = -1, const double nSigma = 1.) const;

      virtual bool GetShift_ccqe_muonPT_muonP(double &shift, const double muonP, const double muonTheta, const double muonPErr, const int nUniverse = -1, const double nSigma = 1.) const; 
 
      virtual bool GetShift_ccqe_muonPZ_muonP(double &shift, const double muonP, const double muonTheta, const double muonPErr, const int nUniverse = -1, const double nSigma = 1.) const; 

      //-----------------------------------
      // Muon theta shifts 
      //-----------------------------------
      virtual bool GetShift_ccqe_muonTheta(double& shift, const double muonThetaErr, const int nUniverse = -1, const double nSigma = 1.) const;

      virtual bool GetShift_ccqe_muonThetaXY(double& shift, const double muonTheta, const double muonThetaX, const double muonThetaX_Err, const double muonThetaY, const double muonThetaY_Err, const int nUniverse = -1, const double nSigma = 1.) const; 

      virtual bool GetShift_ccqe_nuE_muonTheta(double &shift, const double muonP, const double muonTheta, const double bindingE, const double muonThetaErr, const int nUniverse = -1, const double nSigma = 1.) const;

      virtual bool GetShift_ccqe_q2_muonTheta(double &shift, const double muonP, const double muonTheta, const double bindingE, const double muonThetaErr, const int nUniverse = -1, const double nSigma = 1.) const;

      virtual bool GetShift_ccqe_q2_muonThetaXY(double &shift, const double muonP, const double muonTheta, const double bindingE, const double muonThetaX, const double muonThetaX_Err, const double muonThetaY, const double muonThetaY_Err, const int nUniverse = -1, const double nSigma = 1. ) const;  

      virtual bool GetShift_ccqe_muonPT_muonTheta(double &shift, const double muonP, const double muonTheta, const double muonThetaErr, const int nUniverse = -1, const double nSigma = 1.) const;  

      virtual bool GetShift_ccqe_muonPT_muonThetaXY(double &shift, const double muonP, const double muonTheta, const double muonThetaX, const double muonThetaX_Err, const double muonThetaY, const double muonThetaY_Err, const int nUniverse = -1, const double nSigma = 1.) const; 

      virtual bool GetShift_ccqe_muonPZ_muonTheta(double &shift, const double muonP, const double muonTheta, const double muonThetaErr, const int nUniverse = -1, const double nSigma = 1.) const;  

      virtual bool GetShift_ccqe_muonPZ_muonThetaXY(double &shift, const double muonP, const double muonTheta, const double muonThetaX, const double muonThetaX_Err, const double muonThetaY, const double muonThetaY_Err, const int nUniverse = -1, const double nSigma = 1.) const; 

      //----------------------------------
      // Binding energy shifts 
      //----------------------------------
      virtual bool GetShift_ccqe_q2_bindingEnergy(double &shift, const double muonP, const double muonTheta, const double bindingE, const double bindingE_Err, const int nUniverse = -1, const double nSigma = 1.) const; 

      //-----------------------------------
      // Recoil energy shifts 
      //-----------------------------------
      virtual bool GetShift_ccqe_recoilE( double& shift, const double muonRecoilerr, const int nUniverse = -1, const double nSigma = 1. ) const;

      //===========================================
      // END - CCQE Calculators and Scale Shifts
      //===========================================      
      //================================= 
      // BEGIN - Muon Energy Angle Shifts
      //================================= 
      virtual bool GetShifts_muonTheta( dVec& shifts,
          const double muonTheta, const double muonE, const double hadronicE,
          const double muonThetaErr) const;

      /*! @brief Shift in Muon Energy as a function of itself. 
       */
      virtual bool GetShift_muonTheta( double& shift,
          const double MuonTheta, const double muonE, const double hadronicE,
          const double muonThetaErr, 
          const double nSigma = 1.) const;

      //------------------------------------

      //================================= 
      // BEGIN - Muon Energy Scale Shifts
      //================================= 
      //! @name Shifts in Muon Energy for CC Inclusive
      //!@{

      /*! @brief Get vector of shifts (one for each universe) in  muon energy
        @param[out] shifts Vector of shifts - in MeV
        @param[in] muonE The muon energy in units of MeV
        @param[in] hadronicE The recoil energy to be shifted in units of MeV (used to check limits)
        @param[in] muonEerr The uncertainty on muon energy in units of MeV
       */
      virtual bool GetShifts_muonE( dVec& shifts,
          const double muonE, const double hadronicE,
          const double muonEerr) const;

      /*! @brief Shift in Muon Energy as a function of itself. 
       */
      virtual bool GetShift_muonE( double& shift,
          const double muonE, const double hadronicE,
          const double muonEerr, 
          const double nSigma = 1.) const;

      //------------------------------------
      /*! @brief Get vector of shifts (one for each universe) in neutrino energy
        @param[out] shifts Vector of shifts - in MeV
        @param[in] muonE The muon energy in units of MeV
        @param[in] hadronicE The recoil energy to be shifted in units of MeV
        @param[in] muonEerr The uncertainty on muon energy in units of MeV
       */
      virtual bool GetShifts_E_muonE( dVec& shifts,
          const double muonE, const double hadronicE,
          const double muonEerr) const;

      /*! @brief Shift in Neutrino Energy as a function of Muon Energy.
       */
      virtual bool GetShift_E_muonE( double& shift,
          const double muonE, const double hadronicE, 
          const double muonEerr,
          const double nSigma = 1.) const;

      //------------------------------------
      /*! @brief Get vector of shifts (one for each universe) in Q2
        @param[out] shifts Vector of shifts - in MeV^2
        @param[in] muonE The muon energy in units of MeV
        @param[in] hadronicE The recoil energy to be shifted in units of MeV
        @param muonTheta The muon angle in units of radians.
        @param[in] muonEerr The uncertainty on muon energy in units of MeV
       */
      virtual bool GetShifts_q2_muonE( dVec& shifts,
          const double muonE, const double hadronicE, const double muonTheta,
          const double muonEerr) const;

      virtual bool GetShift_q2_muonE( double& shift,
          const double muonE, const double hadronicE, const double muonTheta, 
          const double muonEerr,
          const double nSigma = 1.) const;

      //------------------------------------
      /*! @brief Get vector of shifts (one for each universe) in W
        @param[out] shifts Vector of shifts - in MeV
        @param[in] muonE The muon energy in units of MeV
        @param[in] hadronicE The recoil energy to be shifted in units of MeV
        @param muonTheta The muon angle in units of radians.
        @param[in] muonEerr The uncertainty on muon energy in units of MeV
       */
      virtual bool GetShifts_w_muonE( dVec& shifts,
          const double muonE, const double hadronicE, const double muonTheta,
          const double muonEerr) const;

      /*! @brief Shift in W as a function of Muon Energy
        @param shift The shifted W - passed back by reference in units of MeV.
        @param muonE The muon energy in units of MeV.
        @param hadronicE The recoil energy to be shifted in units of MeV.
        @param muonTheta The muon angle in units of radians.
        @param muonEerr The uncertainty on the Muon energy in units of MeV.
        @param nSigma Result multiplier: use -1 for the negative shifted value.
        @return false for unphysical W
       */
      virtual bool GetShift_w_muonE( double& shift,
          const double muonE, const double hadronicE, const double muonTheta, 
          const double muonEerr,
          const double nSigma = 1.)const;

      //------------------------------------
      /*! @brief Get vector of shifts (one for each universe) in x
        @param[out] shifts Vector of shifts
        @param[in] muonE The muon energy in units of MeV
        @param[in] hadronicE The recoil energy to be shifted in units of MeV
        @param muonTheta The muon angle in units of radians.
        @param[in] muonEerr The uncertainty on muon energy in units of MeV
       */
      virtual bool GetShifts_x_muonE( dVec& shifts,
          const double muonE, const double hadronicE, const double muonTheta,
          const double muonEerr) const;

      virtual bool GetShift_x_muonE( double& shift,
          const double muonE, const double hadronicE, const double muonTheta, const double muonEerr,
          const double nSigma = 1. ) const;

      //-----------------------------------
      /*! @brief Get vector of shifts (one for each universe) in y
        @param[out] shifts Vector of shifts
        @param[in] muonE The muon energy in units of MeV
        @param[in] hadronicE The recoil energy to be shifted in units of MeV
        @param[in] muonEerr The uncertainty on muon energy in units of MeV
       */
      virtual bool GetShifts_y_muonE( dVec& shifts,
          const double muonE, const double hadronicE,
          const double muonEerr) const;

      virtual bool GetShift_y_muonE( double& shift,
          const double muonE, const double hadronicE,
          const double muonEerr,
          const double nSigma = 1. ) const;

      //!@}
      //================================= 
      // END - Muon Energy Shifts
      //================================= 

      //================================= 
      // BEGIN - Hadronic Energy Shifts
      //================================= 
      //! @name Shifts in Hadronic Energy for CC Inclusive
      //!@{

      /*! @brief Get vector of shifts (one for each universe) in  Hadronic (Recoil) Energy
        @param[out] shifts Vector of shifts - in MeV
        @param[in] muonE The muon energy in units of MeV (used to check limits)
        @param[in] hadronicE The recoil energy to be shifted in units of MeV
       */
      virtual bool GetShifts_hadronicE( dVec& shifts,
          const double muonE, const double hadronicE) const;


      /*! @brief Shift in Hadronic (Recoil) Energy as a function of itself - currently supports 1-sigma shifts only!
        @param shift The shifted energy - passed back by reference in units of MeV.
        @param[in] muonE The muon energy in units of MeV (used to check limits)
        @param hadronicE The recoil energy to be shifted in units of MeV.
        @param nSigma Result multiplier: use -1 for the negative shifted value.
       */
      virtual bool GetShift_hadronicE( double& shift,
          const double muonE, const double hadronicE,
          const double nSigma = 1.) const;

      //---------------------------------
      /*! @brief Get vector of shifts (one for each universe) in neutrino energy cause by shifts in hadronic (recoil) Energy
        @param[out] shifts Vector of shifts - in MeV
        @param[in] muonE The muon energy in units of MeV (used to check limits)
        @param[in] hadronicE The recoil energy to be shifted in units of MeV
       */
      virtual bool GetShifts_E_hadronicE( dVec& shifts,
          const double muonE, const double hadronicE) const;

      /*! @brief Shift in Neutrino Energy as a function of Hadronic (Recoil) Energies.
        @param shift The shift in energy - passed back by reference in units of MeV.
        @param muonE The muon energy in units of MeV.
        @param hadronicE The recoil energy in units of MeV.
        @param nSigma Result multiplier: use -1 for the negative shifted value.
       */
      virtual bool GetShift_E_hadronicE( double& shift,
          const double muonE, const double hadronicE,
          const double nSigma = 1.)const;


      //---------------------------------
      /*! @brief Get vector of shifts (one for each universe) in neutrino energy cause by shifts in hadronic (recoil) Energy
        @param[out] shifts Vector of shifts in energy - in MeV
        @param[in] muonE The muon energy in units of MeV (used to check limits)
        @param[in] hadronicE The recoil energy to be shifted in units of MeV
        @param[in] muonTheta The muon angle in units of radians.
       */
      virtual bool GetShifts_q2_hadronicE( dVec& shifts,
          const double muonE, const double hadroicE, const double muonTheta ) const;

      /*! @brief Shift in Q-squared as a function of Hadronic (Recoil) Energies.
        @param shift The shift in Q-squared - passed back by reference in units of MeV-squared.
        @param muonE The muon energy in units of MeV.
        @param hadronicE The recoil energy to be shifted in units of MeV.
        @param muonTheta The muon angle in units of radians.
        @param nSigma Result multiplier: use -1 for the negative shifted value.
       */
      virtual bool GetShift_q2_hadronicE( double& shift,
          const double muonE, const double hadronicE, const double muonTheta,
          const double nSigma = 1.)const;

      //---------------------------------
      /*! @brief Get vector of shifts (one for each universe) in W caused by shifts in hadronic (recoil) Energy
        @param[out] shifts Vector of shifts in - in MeV
        @param[in] muonE The muon energy in units of MeV (used to check limits)
        @param[in] hadronicE The recoil energy to be shifted in units of MeV
        @param[in] muonTheta The muon angle in units of radians.
       */
      virtual bool GetShifts_w_hadronicE( dVec& shifts,
          const double muonE, const double hadroicE, const double muonTheta ) const;


      /*! @brief Shift in W as a function of Hadronic (Recoil) Energies. Return false for unphysical W.
        @param shift The shift in W - passed back by reference in units of MeV.
        @param muonE The muon energy in units of MeV.
        @param hadronicE The recoil energy to be shifted in units of MeV.
        @param muonTheta The muon angle in units of radians. (Needed to calculate Q^2.)
        @param nSigma Result multiplier: use -1 for the negative shifted value.
       */
      virtual bool GetShift_w_hadronicE( double& shift,
          const double muonE, const double hadronicE, const double muonTheta,
          const double nSigma = 1.)const;

      //---------------------------------
      /*! @brief Get vector of shifts (one for each universe) in x caused by shifts in hadronic (recoil) Energy
        @param[out] shifts Vector of shift amounts - in MeV
        @param[in] muonE The muon energy in units of MeV (used to check limits)
        @param[in] hadronicE The recoil energy to be shifted in units of MeV
        @param[in] muonTheta The muon angle in units of radians.
       */
      virtual bool GetShifts_x_hadronicE( dVec& shifts,
          const double muonE, const double hadroicE, const double muonTheta ) const;

      /*! @brief Shift in Bjorken x. Energy units are MeV.
       */
      virtual bool GetShift_x_hadronicE( double& shift,
          const double muonE, const double hadronicE, const double muonTheta, 
          const double nSigma = 1.)const;


      //---------------------------------
      /*! @brief Get vector of shifts (one for each universe) in y caused by shifts in hadronic (recoil) Energy
        @param[out] shifts Vector of shift amounts
        @param[in] muonE The muon energy in units of MeV (used to check limits)
        @param[in] hadronicE The recoil energy to be shifted in units of MeV
       */
      virtual bool GetShifts_y_hadronicE( dVec& shifts,
          const double muonE, const double hadroicE ) const;

      /*! @brief Shift in Inelasticity. Energy units are MeV.
       */
      virtual bool GetShift_y_hadronicE( double& shift,
          const double muonE, const double hadronicE, 
          const double nSigma = 1.)const;

      //---------------------------------
      /*! @brief Calorimetric systematic uncertainty from J. Devan docdb 7605.
        @param nml Neutrino-minus-lepton energy (recoil). Units are GeV!
        @param anti True if this event is from an antineutrino.
       */
      virtual double CCnumu_calorimetric_err( double nml, bool anti ) const;
      virtual double CCnumu_calorimetric_err( double nml ) const { return CCnumu_calorimetric_err( nml, m_isAntiNu ); };

      //!@}
      //================================= 
      // END - Hadronic Energy Shifts
      //================================= 

      //================================= 
      // BEGIN - Muon Theta Shifts
      //================================= 
      //! @name Shifts in Muon Theta for CC Inclusive
      //!@{

      /*! @brief Get vector of shifts (one for each universe) in Muon Theta
        @param[out] shifts Vector of shifts - in rad
        @param[in] muonTheta The muon angle in rad
       */
      virtual bool GetShifts_muonTheta( dVec& shifts,
          const double muonTheta) const;

      /*! @brief Shift in muon theta as a function of itself
        @param shift The shifted muon theta
        @param[in] muonTheta The muon angle in rad
        @param nSigma Result multiplier: use -1 for the negative shifted value.
       */
      virtual bool GetShift_muonTheta( double& shift,
          const double muonTheta,
          const double nSigma = 1.) const;


      //---------------------------------
      /*! @brief Get vector of shifts (one for each universe) in Q2 caused by shifts in muon theta
        @param[out] shifts Vector of shifts in Q2 - in MeV/c^2
        @param[in] muonE The muon energy in units of MeV (used to check limits)
        @param[in] hadronicE The recoil energy to be shifted in units of MeV
        @param[in] muonTheta The muon angle in units of radians.
       */
      virtual bool GetShifts_q2_muonTheta( dVec& shifts,
          const double muonE, const double hadroicE, const double muonTheta ) const;

      /*! @brief Shift in Q-squared as a function of muon theta.
        @param shift The shift in Q-squared - passed back by reference in units of MeV-squared.
        @param muonE The muon energy in units of MeV.
        @param hadronicE The recoil energy to be shifted in units of MeV.
        @param muonTheta The muon angle in units of radians.
        @param nSigma Result multiplier: use -1 for the negative shifted value.
       */
      virtual bool GetShift_q2_muonTheta( double& shift,
          const double muonE, const double hadronicE, const double muonTheta,
          const double nSigma = 1.)const;

      //---------------------------------
      /*! @brief Get vector of shifts (one for each universe) in W caused by shifts in muon theta
        @param[out] shifts Vector of shifts in - in rad
        @param[in] muonE The muon energy in units of MeV (used to check limits)
        @param[in] hadronicE The recoil energy to be shifted in units of MeV
        @param[in] muonTheta The muon angle in units of radians.
       */
      virtual bool GetShifts_w_muonTheta( dVec& shifts,
          const double muonE, const double hadroicE, const double muonTheta ) const;


      /*! @brief Shift in W as a function of muon theta. Return false for unphysical W.
        @param shift The shift in W - passed back by reference in units of MeV.
        @param muonE The muon energy in units of MeV.
        @param hadronicE The recoil energy to be shifted in units of MeV.
        @param muonTheta The muon angle in units of radians. (Needed to calculate Q^2.)
        @param nSigma Result multiplier: use -1 for the negative shifted value.
       */
      virtual bool GetShift_w_muonTheta( double& shift,
          const double muonE, const double hadronicE, const double muonTheta,
          const double nSigma = 1.)const;

      //---------------------------------
      /*! @brief Get vector of shifts (one for each universe) in x caused by shifts in muon theta
        @param[out] shifts Vector of shifts in x
        @param[in] muonE The muon energy in units of MeV (used to check limits)
        @param[in] hadronicE The recoil energy to be shifted in units of MeV
        @param[in] muonTheta The muon angle in units of radians.
       */
      virtual bool GetShifts_x_muonTheta( dVec& shifts,
          const double muonE, const double hadroicE, const double muonTheta ) const;

      /*! @brief Shift in x as a function of muon theta.
        @param shift The shift in x
        @param muonE The muon energy in units of MeV. 
        @param hadronicE The recoil energy to be shifted in units of MeV.
        @param muonTheta The muon angle in units of radians.
        @param nSigma Result multiplier: use -1 for the negative shifted value.
       */
      virtual bool GetShift_x_muonTheta( double& shift,
          const double muonE, const double hadronicE, const double muonTheta,
          const double nSigma = 1.)const;


      /*! @brief Muon theta uncertainties stores internally
        @note Use a function because muon theta uncertainty might depend on stuff like the angle, hadronic E, if the vertex is fit...
       */
      virtual double GetMuonTheta_err( ) const;

      //!@}
      //================================= 
      // END - Muon Theta Shifts
      //================================= 


      //===========================================================
      // BEGIN - Get a Vector of Gaus Random Initialized variables
      //         for Lateral multiuniverses filling purposes
      //===========================================================
      bool GetRandomShiftVector(std::vector<double>& rShift, const double seed, const int nUniverses = 1000 );
      bool GetRandomShiftVectors(const int nUniverses = 1000 );

      //============================================================
      // END - Get a Vector of Gaus Random Initialized variables
      //============================================================

      //============================================================
      // BEGIN - energy limit get and set functions
      //============================================================

      inline double GetLowerLimit_E() const     { return m_lowerLimit_E; };
      inline void   SetLowerLimit_E( double e ) { m_lowerLimit_E = e;    };
      inline double GetUpperLimit_E() const     { return m_upperLimit_E; };
      inline void   SetUpperLimit_E( double e ) { m_upperLimit_E = e;    };

      inline double GetLowerLimit_hadronicE() const     { return m_lowerLimit_hadronicE; };
      inline void   SetLowerLimit_hadronicE( double e ) { m_lowerLimit_hadronicE = e;    };
      inline double GetUpperLimit_hadronicE() const     { return m_upperLimit_hadronicE; };
      inline void   SetUpperLimit_hadronicE( double e ) { m_upperLimit_hadronicE = e;    };

      inline double GetLowerLimit_muonE() const     { return m_lowerLimit_muonE; };
      inline void   SetLowerLimit_muonE( double e ) { m_lowerLimit_muonE = e;    };
      inline double GetUpperLimit_muonE() const     { return m_upperLimit_muonE; };
      inline void   SetUpperLimit_muonE( double e ) { m_upperLimit_muonE = e;    };


      inline double GetLowerLimit_muonTheta() const     { return m_lowerLimit_muonTheta; };
      inline void   SetLowerLimit_muonTheta( double t ) { m_lowerLimit_muonTheta = t;    };
      inline double GetUpperLimit_muonTheta() const     { return m_upperLimit_muonTheta; };
      inline void   SetUpperLimit_muonTheta( double t ) { m_upperLimit_muonTheta = t;    };

      //============================================================
      // END - energy limit get and set functions
      //============================================================

      //============================================================
      // BEGIN - out of bounds tests
      //============================================================
      inline bool OutOfBounds_E( double e )         const { return ( e < m_lowerLimit_E         || m_upperLimit_E < e ); };
      inline bool OutOfBounds_hadronicE( double e ) const { return ( e < m_lowerLimit_hadronicE || m_upperLimit_hadronicE < e ); };
      inline bool OutOfBounds_muonE( double e )     const { return ( e < m_lowerLimit_muonE     || m_upperLimit_muonE < e ); };
      inline bool OutOfBounds_muonTheta( double t )     const { return ( t < m_lowerLimit_muonTheta     || m_upperLimit_muonTheta < t ); };
      //============================================================
      // END - out of bounds tests
      //============================================================


    protected:

      void SetSeeds(); ///< Set default seeds
      void SetDefaultEnergyLimits(); ///< Set default energy limits which define the not physical region

      bool m_isAntiNu; ///< Are we computing numbers for anti-neutrinos by default?
      UInt_t m_seedMuon;      ///< The seed for random number generator for muon energy shifts
      UInt_t m_seedHadron;    ///< The seed for random number generator for hadron energy shifts
      UInt_t m_seedMuonTheta; ///< The seed for random number generator for muon theta shifts

      //UInt_t m_seedBeamX;
      //UInt_t m_seedBeamY;
      //UInt_t m_seedBeamZ;
      UInt_t m_seedBeamThetaX;
      UInt_t m_seedBeamThetaY;

      UInt_t m_seedMuonX;
      UInt_t m_seedMuonY;
      UInt_t m_seedMuonZ;

      UInt_t m_seedHadronX;
      UInt_t m_seedHadronY;
      UInt_t m_seedHadronZ;

      int m_nUniverses;    ///< How many universes of shifts will I generate?

      double m_lowerLimit_E; ///< A event with neutrino energy lower than this is "not physical" (MeV)
      double m_upperLimit_E; ///< A event with neutrino energy greater than this is "not physical" (MeV)

      double m_lowerLimit_hadronicE; ///< A event with hadronic energy lower than this is "not physical" (MeV)
      double m_upperLimit_hadronicE; ///< A event with hadronic energy greater than this is "not physical" (MeV)

      double m_lowerLimit_muonE; ///< A event with muon energy lower than this is "not physical" (MeV)
      double m_upperLimit_muonE; ///< A event with muon energy greater than this is "not physical" (MeV)

      double m_lowerLimit_muonTheta; ///< A event with muon theta lower than this is "not physical" (rad)
      double m_upperLimit_muonTheta; ///< A event with muon theta greater than this is "not physical" (rad)

      dVec m_rShiftMuon;      ///< Vector of random shifts for muon in units of nSigmas
      dVec m_rShiftHadron;    ///< Vector of random shifts for hadron in units of nSigmas
      dVec m_rShiftMuonTheta; ///< Vector of random shifts for muon theta in units of nSigmas
      //dVec m_rShiftBeamDirX;
      //dVec m_rShiftBeamDirY;
      //dVec m_rShiftBeamDirZ;
      dVec m_rShiftBeamThetaX;
      dVec m_rShiftBeamThetaY;

      dVec m_rShiftMuonDirX;
      dVec m_rShiftMuonDirY;
      dVec m_rShiftMuonDirZ;
      dVec m_rShiftHadronDirX;
      dVec m_rShiftHadronDirY;
      dVec m_rShiftHadronDirZ;

  }; //end of MnvRecoShifter

} //end of PlotUtils

#endif
