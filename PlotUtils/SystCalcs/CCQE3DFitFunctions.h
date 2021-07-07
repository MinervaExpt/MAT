#ifndef CCQE3DFITFUNCTIONS_H
#define CCQE3DFITFUNCTIONS_H

double GetCCQE3DFitsWeight( ) const {
  return 1.0;
}

double GetPmuTransverseTrue() const { /* GeV */

  const double px=GetVecElem("mc_primFSLepton", 0);
  const double py=GetVecElem("mc_primFSLepton", 1);
  const double pz=GetVecElem("mc_primFSLepton", 2);
  // Copied (and modified) from MinervaCoordSysTool::thetaWRTBeam
  const double numi_beam_angle_rad = -0.05887; // Taken from MinervaPhysicalConstants.h
  double pyp = -1.0 *sin( numi_beam_angle_rad )*pz + cos( numi_beam_angle_rad )*py;
  double pt   = sqrt( pow(px,2) +pow(pyp,2) );
  
  return 1.0e-3*pt; //GeV

}

double GetPmuLongitudinalTrue() const { /* GeV */

  const double pylep=GetVecElem("mc_primFSLepton", 1);
  const double pzlep=GetVecElem("mc_primFSLepton", 2);
  // Copied (and modified) from MinervaCoordSysTool::thetaWRTBeam
  const double numi_beam_angle_rad = -0.05887; // Taken from MinervaPhysicalConstants.h
  double pzp = cos( numi_beam_angle_rad )*pzlep + sin( numi_beam_angle_rad )*pylep;
  
  return 1.0e-3*pzp; //GeV

}

double GetEAvailableTrue() const { /* MeV */

  double recoil = 0;
  int n_parts = GetInt("mc_nFSPart");
  double mass_pion = 135;
  double mass_proton = 938.27;
  for(int i=0;i<n_parts;i++){
    int pdg = GetVecElemInt("mc_FSPartPDG",i);
    if(pdg == 22) recoil+=GetVecElem("mc_FSPartE",i);//total energy
    if(pdg == 211 || pdg == -211) recoil+=GetVecElem("mc_FSPartE",i)-mass_pion;//kinetic
    if(pdg == 111) recoil+=GetVecElem("mc_FSPartE",i);//total energy
    if(pdg == 2212) recoil+=GetVecElem("mc_FSPartE",i)-mass_proton;//kinetic
  } 

  return recoil;

}

#endif  // CCQE3DFITFUNCTIONS
