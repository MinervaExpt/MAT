#ifndef genie_particle_cxx
#define genie_particle_cxx

#include "genie_particle.h"
using namespace MAT;

TLorentzVector genie_particle::GetMomentumEnergy() const {
    return __p4;
}

TLorentzVector genie_particle::p4() const {
    return __p4;
}

void genie_particle::SetMomentumEnergy(TLorentzVector p4) {
    __p4 = p4;

    double px   = p4.X();
    double py   = p4.Y();
    double pz   = p4.Z();
    double etot = p4.E();

    double p = sqrt(px*px + py*py + pz*pz);
    double m = sqrt(etot*etot - p*p);

    __kinetic_energy = etot - m;
    
}

double genie_particle::GetKineticEnergy() const {
    return __kinetic_energy;
}

double genie_particle::T() const {
    return __kinetic_energy;
}

bool genie_particle::default_particle() const {
    return __id == -1 && __pdg == -1 && __status == -1 && __mother == -1;
    
}

void genie_particle::Print()
{
    std::cout.setf(std::ios_base::fixed);
    std::cout.precision(4);
    std::cout << "Particle PDG,ID,MotherID: "<< __pdg << "," << __id <<","<< __mother << std::endl;
    std::cout << "Particle 4-momentum: (" << __p4.Px() << "," << __p4.Py() << "," << __p4.Pz() << "," << __p4.E()
              << ")" << std::endl;
}



#endif
