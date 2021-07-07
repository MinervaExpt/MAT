
#ifndef BASEUNIVERSE_cxx
#define BASEUNIVERSE_cxx

#include "BaseUniverse.h"

using namespace PlotUtils;

bool BaseUniverse::m_is_truth                  = false;

// CTOR
  BaseUniverse::BaseUniverse(PlotUtils::TreeWrapper* chw,
                                       double nsigma )
    : m_chw(chw), m_nsigma(nsigma), m_entry(-1)
  {}

//! Generic Branch Getters
  int BaseUniverse::GetInt(const char* name) const {
    return m_chw->GetInt(name, m_entry); 
  }


  bool BaseUniverse::GetBool(const char* name) const {
    return m_chw->GetValue(name, m_entry);
  }


  double BaseUniverse::GetDouble(const char* name) const {
    return m_chw->GetValue(name, m_entry); 
  }


  int BaseUniverse::GetVecElemInt(const char* name, const int i) const {
    return m_chw->GetValue(name, m_entry, i);
  }


  double BaseUniverse::GetVecElem(const char* name, const int i) const {
    return m_chw->GetValue(name, m_entry, i); 
  }

  double BaseUniverse::GetVecElem(const char* name,
                                       const int i,
                                       const int j) const {
    std::vector<std::vector<double>> vec = GetVecOfVecDouble(name);
    return vec[i][j];
  }

  std::vector<int> BaseUniverse::GetVecInt(const char* name) const {
    return m_chw->GetValueVector<int>(name, m_entry);
  }

  std::vector<double> BaseUniverse::GetVecDouble(const char* name) const {
    return m_chw->GetValueVector<double>(name, m_entry);
  }

  std::vector<std::vector<double> > BaseUniverse::GetVecOfVecDouble(const char* name) const {
    return m_chw->GetValueNDVector<std::vector<std::vector<double>>> (name,m_entry);
    //return GetVec<std::vector<double>> (name);
  }

// helper functions
  double BaseUniverse::phi3D(double thetaX,double thetaY) const {
    return std::atan2(std::tan(thetaY),std::tan(thetaX));
  }

  double BaseUniverse::theta3D(double thetaX, double thetaY) const {
    double sec_thetaX = 1./std::cos(thetaX);
    double sec_thetaY = 1./std::cos(thetaY);
    double inter = std::sqrt(1./(sec_thetaX*sec_thetaX + sec_thetaY*sec_thetaY -1.));
    double angle = std::acos(inter);

    if(thetaX > 3.14159265/2 || thetaY > 3.14159265/2 ) return 3.14159265 - angle;
    else return angle;
  }

  double BaseUniverse::calcq0(const double Enu, const double Elep) const {
    return Enu-Elep;
  }

  double BaseUniverse::calcq3(const double Q2, const double Enu, const double Elep) const {
    return sqrt(Q2 + pow(Enu - Elep,2.0)); 
  }

#endif // BaseUniverse_cxx
