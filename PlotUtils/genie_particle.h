#ifndef genie_particle_h
#define genie_particle_h

#include <iostream>

#include <cmath>
#include <iomanip>
#include <TLorentzVector.h>


using std::abs;
using std::setw;
using std::sqrt;

namespace MAT{

class genie_particle {
public:
    genie_particle ()
        : __id(-1),
        __pdg(-1),
        __status(-1),
        __first(-1),
        __last(-1),
        __mother(-1) {}
    
    genie_particle(int id, int pdg, int status, int first, int last, int mother)
        : __id(id),
        __pdg(pdg),
        __status(status),
        __first(first),
        __last(last),
        __mother(mother) {}

    int GetId()      const { return __id;     }
    int GetPDGCode() const { return __pdg;    }
    int GetStatus()  const { return __status; }
    int GetFD()      const { return __first;  }
    int GetLD()      const { return __last;   }
    int GetMother()  const { return __mother; }

    TLorentzVector GetMomentumEnergy() const;
    TLorentzVector p4() const;
    void SetMomentumEnergy(TLorentzVector);

    double GetKineticEnergy() const;
    double T() const;

    bool default_particle() const;
    
    void Print();

    bool operator==(const genie_particle& rhs)
    {
        return this->GetId() == rhs.GetId()
        && this->p4().Px() == rhs.p4().Px()
        && this->p4().Py() == rhs.p4().Py()
        && this->p4().Pz() == rhs.p4().Pz();
    }

    bool operator!=( const genie_particle& rhs)
    {
        return this->GetId() != rhs.GetId()
        || this->p4().Px() != rhs.p4().Px()
        || this->p4().Py() != rhs.p4().Py()
        || this->p4().Pz() != rhs.p4().Pz();
    }   
private:

    int __id;
    int __pdg;
    int __status;
    int __first;
    int __last;
    int __mother;

    TLorentzVector __p4;
    double __kinetic_energy;
};

/*
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

bool operator==(const genie_particle& lhs, const genie_particle& rhs)
{
    return lhs.GetId() == rhs.GetId()
    && lhs.p4().Px() == rhs.p4().Px()
    && lhs.p4().Py() == rhs.p4().Py()
    && lhs.p4().Pz() == rhs.p4().Pz();
}

bool operator!=(const genie_particle& lhs, const genie_particle& rhs)
{
    return lhs.GetId() != rhs.GetId()
    || lhs.p4().Px() != rhs.p4().Px()
    || lhs.p4().Py() != rhs.p4().Py()
    || lhs.p4().Pz() != rhs.p4().Pz();
}
*/

}
#endif
