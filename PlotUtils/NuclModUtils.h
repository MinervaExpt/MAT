#ifndef MNV_NuclModUtils_h
#define MNV_NuclModUtils_h 1

namespace PlotUtils
{

  //==========================================
  // NuclModUtils class
  //==========================================
  //! Class to calculate nuclear modification of structure functions
  class NuclModUtils
  {
    public:
      //! Default constructor
      NuclModUtils() {};

      //! Default destructor
      virtual ~NuclModUtils() {};

      //! singleton gettor
      static NuclModUtils& Get();

      //! Bodek-Yang scaling variable xi_tm used for 2013 nuclei
      double BY_ScalingVar_Xi_TM( double x, double q2 ) const;
      double BY_ScalingVar_Xi_w( double x, double q2 ) const;

      //! Bodek-Yang 03/13 equation for F2_D/F2_(n+p)
      double BY_Free_to_D( double x ) const;
      //! Bodek-Yang 03 equation for F2_Fe/F2_D
      double BY_D_to_Fe( double x ) const;

      //! Bodek-Yang 13 equation for F2_Fe/F2_D
      double BY13_D_to_Fe( double x) const;
      //! Bodek-Yang 13 equation for F2_Fe/F2_C
      double BY13_C_to_Fe( double x) const;
      //! Bodek-Yang 13 equation for F2_Pb/F2_Fe
      double BY13_Fe_to_Pb( double x) const;

      //! E139 fit for F2_A/F2_D
      double E139( double x, double A ) const;
      //! E139 fit for F2_A/F2_D
      double E139( double x, double A, bool &goodfit ) const;

      //! GENIE correction for BY03
      double GENIE_BY( double x, int A ) const;
      //! GENIE correction for BY03 with the bug
      double GENIE_BY_bug( double x, double q2, int A ) const;
      //! GENIE correction for BY13
      double GENIE_BY13( double x, double q2, int A ) const;
      //! GENIE correction for E139
      double GENIE_E139( double x, int A ) const;

  };//end of NuclModUtils

}//end PlotUtils namespace

#endif //Mnv_NuclModUtils
