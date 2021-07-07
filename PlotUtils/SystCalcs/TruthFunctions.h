#ifndef TRUTHFUNCTIONS_H
#define TRUTHFUNCTIONS_H

//! Get true event quantities -- these are needed to calculate several
//! event weights
double GetEnuTrue() const { /* MeV */
  return GetDouble("mc_incomingE");
}

double GetElepTrue() const { /* MeV */
  return GetVecElem("mc_primFSLepton", 3);
}

double GetPlepTrue() const { /* MeV */
  TVector3 p3lep(GetVecElem("mc_primFSLepton", 0),
                 GetVecElem("mc_primFSLepton", 1),
                 GetVecElem("mc_primFSLepton", 2));
  return p3lep.Mag();
}

double GetThetalepTrue() const { /* radians w.r.t. incident nu dirn */
  TVector3 p3lep(GetVecElem("mc_primFSLepton", 0),
                 GetVecElem("mc_primFSLepton", 1),
                 GetVecElem("mc_primFSLepton", 2));
  p3lep.RotateX(MinervaUnits::numi_beam_angle_rad);
  return p3lep.Theta();
}

double GetPhilepTrue() const { /* radians w.r.t. incident nu dirn */
  TVector3 p3lep(GetVecElem("mc_primFSLepton", 0),
                 GetVecElem("mc_primFSLepton", 1),
                 GetVecElem("mc_primFSLepton", 2));
  p3lep.RotateX(MinervaUnits::numi_beam_angle_rad);
  return p3lep.Phi();
}

double GetQ2True() const { /* MeV^2 */
  return GetDouble("mc_Q2");
}

double Getq0True() const { /* MeV */
  return calcq0(GetEnuTrue(), GetElepTrue());
}

double Getq3True() const { /* MeV */
  return calcq3(GetQ2True(), GetEnuTrue(), GetElepTrue());
}

int GetTargetZTrue() const { /* atomic number (Z) of struck nucleus */
  return GetInt("mc_targetZ");
}

double GetVertexZTrue() const { /* longitudinal coord of interaction vertex */
  return GetVecElem("mc_vtx",2);
}

#endif  // TRUTHFUNCTIONS
