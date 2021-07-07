#include <iostream>
#include <vector>

#include <TMath.h>

#include "ChainWrapper.h"

#include "weight_fsi_absorption.h"


double PlotUtils::gethAFSIAbsorptionWeight(int mc_incoming,
                                           int mc_primaryLepton,
                                           int mc_charm,
                                           int mc_intType,
                                           int mc_targetA,
                                           int mc_targetZ,
                                           int mc_er_nPart,
                                           const int* mc_er_ID,
                                           const int* mc_er_status,
                                           const int* mc_er_FD,
                                           const int* mc_er_LD,
                                           const int* mc_er_mother,
                                           const double* mc_er_Px,
                                           const double* mc_er_Py,
                                           const double* mc_er_Pz,
                                           const double* mc_er_E,
                                           bool verbose)
{
    const Double_t protonmass = 938.272046;
    const Double_t neutronmass = 939.565379;
    const Double_t pionmass = 139.57061;
    const Double_t pizeromass = 134.9770;
    
    Int_t proton14N = 0;  
    Int_t neutron14N = 0;  
    Int_t piplus14N = 0; 
    Int_t piminus14N = 0;  
    Int_t pizero14N = 0;
    
    Double_t proton14KE[4]  = {0.0, 0.0, 0.0, 0.0}; 
    Double_t neutron14KE[4] = {0.0, 0.0, 0.0, 0.0}; 
    Double_t piplus14KE[4]  = {0.0, 0.0, 0.0, 0.0}; 
    Double_t piminus14KE[4] = {0.0, 0.0, 0.0, 0.0}; 
    Double_t pizero14KE[4]  = {0.0, 0.0, 0.0, 0.0}; 
    
    Int_t proton14fate[4]  = {0, 0, 0, 0}; 
    Int_t neutron14fate[4] = {0, 0, 0, 0}; 
    Int_t piplus14fate[4]  = {0, 0, 0, 0}; 
    Int_t piminus14fate[4] = {0, 0, 0, 0}; 
    Int_t pizero14fate[4]  = {0, 0, 0, 0}; 
   
    double weightPionAbsorptionBug = 1.0;
        //double weightElasticFSIforQE = 1.0;

    
      
    for(Int_t i=0; i < mc_er_nPart; ++i){
	if(mc_er_status[i] == 14){
            Int_t tempfate = -1;
            
                // Only track fates for proton, neutron, pizero, piplus, piminus.
            if(mc_er_ID[i] == 2212 || mc_er_ID[i] == 2112 || mc_er_ID[i] == 111 || TMath::Abs(mc_er_ID[i])==211){
                
                    // The distinguishing feature is first daughter last daughter
                    // get an easy handle to these for use later
                Int_t fd = mc_er_FD[i];
                Int_t ld = mc_er_LD[i];

                    // is one of these a pion
                Int_t isaPion = 0;
                if(mc_er_ID[i] == 111 || TMath::Abs(mc_er_ID[i] == 211))isaPion = 1;
                
                    // how many daughters are pions
                Int_t FSpions = 0;
                for(int jj = mc_er_FD[i]; jj <= mc_er_LD[i]; ++jj){
                    if(TMath::Abs(mc_er_ID[jj]) == 211 || mc_er_ID[jj] == 111)FSpions++;
                }
                
                if(mc_er_FD[i] == mc_er_LD[i]){
                        //fate is either 1 = no interaction or 3 = elastic
                        // or 5 = absorption with a 3xxxxxxxx hadblob 
                    
                        // This test of 25.00000 MeV is to separate elastic
                        //TODO replace this with a nucleus specific offset
                        // unknown nuclei get an offset like 8 MeV
                    Double_t offset = 25.000000;
                    if(mc_intType != 1)offset = 0.0;
                    if(mc_er_ID[fd] > 2000000000){
                        tempfate = 5;   // absorption on 3 nucleons
                    } else if(TMath::Abs(mc_er_E[i] - mc_er_E[fd] - offset) < 0.00001){
                        tempfate = 1;   // no scattering
                    } else if(mc_er_ID[fd] != mc_er_ID[i]){
                        tempfate = 2;   // charge exchange ?
                    } else {
                        tempfate = 3;   // elastic fate
                    }
                } else if(FSpions > 0 && !isaPion){
                    tempfate = 8;   // nucleon goes to pions
                } else if(isaPion && FSpions == 0){
                    tempfate = 5;
                        // was this two body absorption?
                    if(ld-fd == 1){
                        int p = 0;
                        int n = 0;
                        if(mc_er_ID[fd] == 2212)p++;
                        if(mc_er_ID[fd] == 2112)n++;
                        if(mc_er_ID[ld] == 2212)p++;
                        if(mc_er_ID[ld] == 2112)n++;
                            // use special codes for producing nn, pn, pp
                        if(p==1 && n==1)tempfate = 51;
                        if(p==2)tempfate = 52;
                        if(n==2)tempfate = 50;
                        
                            //std::cout << "Pion Absorption2? " << mc_er_ID[fd] << " " << mc_er_ID[ld] << " " << fd << " " << ld;
                            //if(tempfate==51)std::cout << " got pn ";
                            //if(tempfate==52)std::cout << " got pp ";
                            //if(tempfate==50)std::cout << " got nn ";
                            //std::cout << std::endl;
                    }
                } else if(isaPion && FSpions >= 2){
                    tempfate = 8;
                } else if(isaPion && FSpions == 1){
                    tempfate = 4;
                    for(int jj = mc_er_FD[i]; jj <= mc_er_LD[i]; ++jj){
                        if((mc_er_ID[jj] == 111 || TMath::Abs(mc_er_ID[jj]) == 211) && mc_er_ID[jj] != mc_er_ID[i]) tempfate = 2;
                    }
                } else {
                    if (verbose) std::cout << " Got unknown fate " << tempfate << std::endl;
                    tempfate = 4;  // or four, can't tell.
                }

                    //std::cout << "fate inttype " << mc_intType << " fd " << mc_er_ID[fd] << " "  << mc_er_E[fd] << " ld " << mc_er_ID[ld] << " " << mc_er_E[ld] << " orig " << mc_er_ID[i] << " " << mc_er_E[i] << " Ediff " << mc_er_E[i]-mc_er_E[fd] << " fate " << tempfate << std::endl;
                    //std::cout << "fate inttype " << mc_intType << " nD " << ld - fd + 1 << " fate " << tempfate <<  "      fd " << mc_er_ID[fd] << " "  << mc_er_E[fd] << " ld " << mc_er_ID[ld] << " " << mc_er_E[ld] << " orig " << mc_er_ID[i] << " " << mc_er_E[i] << " Ediff " << mc_er_E[i]-mc_er_E[fd] << std::endl;
                
                
            }
            
            
                // The next block of code reorders the list of hadrons
                // So that the first four are in decreasing KE order
                // (and more than four drop out of the array, because I'm low-recoiling).
            
            if(mc_er_ID[i] == 2212){
                Double_t myKE = mc_er_E[i] - protonmass;
                    // What was this protons fate ?
                
                
                    // can I get back the FSI code.
                if(proton14N == 0){
                    proton14KE[0] = myKE;
                    proton14fate[0] = tempfate;
                    proton14N++;
                } else if(proton14KE[0] < myKE){
                    proton14KE[3] = proton14KE[2]; proton14KE[2] = proton14KE[1]; proton14KE[1] = proton14KE[0];
                    proton14fate[3] = proton14fate[2]; proton14fate[2] = proton14fate[1]; proton14fate[1] = proton14fate[0];
                    proton14KE[0] = myKE;
                    proton14fate[0] = tempfate;
                    proton14N++;
                } else if(proton14KE[1] < myKE){
                    proton14KE[3] = proton14KE[2]; proton14KE[2] = proton14KE[1];
                    proton14fate[3] = proton14fate[2]; proton14fate[2] = proton14fate[1];
                    proton14KE[1] = myKE;
                    proton14fate[1] = tempfate;
                    proton14N++;
                } else if(proton14KE[2] < myKE){
                    proton14KE[3] = proton14KE[2];
                    proton14fate[3] = proton14fate[2]; 
                    proton14KE[2] = myKE;
                    proton14fate[2] = tempfate;
                    proton14N++;
                } else if(proton14KE[3] < myKE){
                    proton14KE[3] = myKE;
                    proton14fate[3] = tempfate;
                    proton14N++;
                }
                    //std::cout << "proton14 N " << proton14N << " E " << proton14KE[0] << " " << proton14KE[1] << " " << proton14KE[2] << " " << proton14KE[3] << std::endl;
                    // end proton	      
            } else if(mc_er_ID[i] == 2112){
                Double_t myKE = mc_er_E[i] - neutronmass;
                    // can I get back the FSI code.
                if(neutron14N == 0){
                    neutron14KE[0] = myKE;
                    neutron14fate[0] = tempfate;
                    neutron14N++;
                } else if(neutron14KE[0] < myKE){
                    neutron14KE[3] = neutron14KE[2]; neutron14KE[2] = neutron14KE[1]; neutron14KE[1] = neutron14KE[0];
                    neutron14fate[3] = neutron14fate[2]; neutron14fate[2] = neutron14fate[1]; neutron14fate[1] = neutron14fate[0];
                    neutron14KE[0] = myKE;
                    neutron14fate[0] = tempfate;
                    neutron14N++;
                } else if(neutron14KE[1] < myKE){
                    neutron14KE[3] = neutron14KE[2]; neutron14KE[2] = neutron14KE[1];
                    neutron14fate[3] = neutron14fate[2]; neutron14fate[2] = neutron14fate[1];
                    neutron14KE[1] = myKE;
                    neutron14fate[1] = tempfate;
                    neutron14N++;
                } else if(neutron14KE[2] < myKE){
                    neutron14KE[3] = neutron14KE[2];
                    neutron14fate[3] = neutron14fate[2]; 
                    neutron14KE[2] = myKE;
                    neutron14fate[2] = tempfate;
                    neutron14N++;
                } else if(neutron14KE[3] < myKE){
                    neutron14KE[3] = myKE;
                    neutron14fate[3] = tempfate;
                    neutron14N++;
                }
                    //std::cout << "neutron14 N " << neutron14N << " E " << neutron14KE[0] << " " << neutron14KE[1] << " " << neutron14KE[2] << " " << neutron14KE[3] << std::endl;
	  
            } else if(mc_er_ID[i] == 111){
                Double_t myKE = mc_er_E[i] - pizeromass;
                    // can I get back the FSI code.
                if(pizero14N == 0){
                    pizero14KE[0] = myKE;
                    pizero14fate[0] = tempfate;
                    pizero14N++;
                } else if(pizero14KE[0] < myKE){
                    pizero14KE[3] = pizero14KE[2]; pizero14KE[2] = pizero14KE[1]; pizero14KE[1] = pizero14KE[0];
                    pizero14fate[3] = pizero14fate[2]; pizero14fate[2] = pizero14fate[1]; pizero14fate[1] = pizero14fate[0];
                    pizero14KE[0] = myKE;
                    pizero14fate[0] = tempfate;
                    pizero14N++;
                } else if(pizero14KE[1] < myKE){
                    pizero14KE[3] = pizero14KE[2]; pizero14KE[2] = pizero14KE[1];
                    pizero14fate[3] = pizero14fate[2]; pizero14fate[2] = pizero14fate[1];
                    pizero14KE[1] = myKE;
                    pizero14fate[1] = tempfate;
                    pizero14N++;
                } else if(pizero14KE[2] < myKE){
                    pizero14KE[3] = pizero14KE[2];
                    pizero14fate[3] = pizero14fate[2]; 
                    pizero14KE[2] = myKE;
                    pizero14fate[2] = tempfate;
                    pizero14N++;
                } else if(pizero14KE[3] < myKE){
                    pizero14KE[3] = myKE;
                    pizero14fate[3] = tempfate;
                    pizero14N++;
                }
                    //std::cout << "pizero14 N " << pizero14N << " E " << pizero14KE[0] << " " << pizero14KE[1] << " " << pizero14KE[2] << " " << pizero14KE[3] << std::endl;
	  
            } else if(mc_er_ID[i] == 211){
                Double_t myKE = mc_er_E[i] - pionmass;
                    // can I get back the FSI code.
                if(piplus14N == 0){
                    piplus14KE[0] = myKE;
                    piplus14fate[0] = tempfate;
                    piplus14N++;
                } else if(piplus14KE[0] < myKE){
                    piplus14KE[3] = piplus14KE[2]; piplus14KE[2] = piplus14KE[1]; piplus14KE[1] = piplus14KE[0];
                    piplus14fate[3] = piplus14fate[2]; piplus14fate[2] = piplus14fate[1]; piplus14fate[1] = piplus14fate[0];
                    piplus14KE[0] = myKE;
                    piplus14fate[0] = tempfate;
                    piplus14N++;
                } else if(piplus14KE[1] < myKE){
                    piplus14KE[3] = piplus14KE[2]; piplus14KE[2] = piplus14KE[1];
                    piplus14fate[3] = piplus14fate[2]; piplus14fate[2] = piplus14fate[1];
                    piplus14KE[1] = myKE;
                    piplus14fate[1] = tempfate;
                    piplus14N++;
                } else if(piplus14KE[2] < myKE){
                    piplus14KE[3] = piplus14KE[2];
                    piplus14fate[3] = piplus14fate[2]; 
                    piplus14KE[2] = myKE;
                    piplus14fate[2] = tempfate;
                    piplus14N++;
                } else if(piplus14KE[3] < myKE){
                    piplus14KE[3] = myKE;
                    piplus14fate[3] = tempfate;
                    piplus14N++;
                }
                    //std::cout << "pion14 N " << pion14N << " E " << pion14KE[0] << " " << pion14KE[1] << " " << pion14KE[2] << " " << pion14KE[3] << std::endl;
	  
            } else if(mc_er_ID[i] == -211){
                Double_t myKE = mc_er_E[i] - pionmass;
                    // can I get back the FSI code.
                if(piminus14N == 0){
                    piminus14KE[0] = myKE;
                    piminus14fate[0] = tempfate;
                    piminus14N++;
                } else if(piminus14KE[0] < myKE){
                    piminus14KE[3] = piminus14KE[2]; piminus14KE[2] = piminus14KE[1]; piminus14KE[1] = piminus14KE[0];
                    piminus14fate[3] = piminus14fate[2]; piminus14fate[2] = piminus14fate[1]; piminus14fate[1] = piminus14fate[0];
                    piminus14KE[0] = myKE;
                    piminus14fate[0] = tempfate;
                    piminus14N++;
                } else if(piminus14KE[1] < myKE){
                    piminus14KE[3] = piminus14KE[2]; piminus14KE[2] = piminus14KE[1];
                    piminus14fate[3] = piminus14fate[2]; piminus14fate[2] = piminus14fate[1];
                    piminus14KE[1] = myKE;
                    piminus14fate[1] = tempfate;
                    piminus14N++;
                } else if(piminus14KE[2] < myKE){
                    piminus14KE[3] = piminus14KE[2];
                    piminus14fate[3] = piminus14fate[2]; 
                    piminus14KE[2] = myKE;
                    piminus14fate[2] = tempfate;
                    piminus14N++;
                } else if(piminus14KE[3] < myKE){
                    piminus14KE[3] = myKE;
                    piminus14fate[3] = tempfate;
                    piminus14N++;
                }
	    //std::cout << "pion14 N " << pion14N << " E " << pion14KE[0] << " " << pion14KE[1] << " " << pion14KE[2] << " " << pion14KE[3] << std::endl;
                
            }
	}   // end status14
    }
    
        // Now identify what weights we need to use.
    if(mc_intType==1 && !mc_charm){
            // Get the QE elastic fate weight
	bool isanu = false;
	int myfate = -1;
	double myKE = -999.0;
	
	if( mc_incoming == 14 ||  mc_incoming == 12 ||  mc_incoming == 16 ){
            myfate = proton14fate[0];
            myKE = proton14KE[0];
	  
	} else if ( mc_incoming == -14 ||  mc_incoming == -12 ||  mc_incoming == -16 ){
            myfate = neutron14fate[0];
            myKE = neutron14KE[0];
            isanu = true;
            
	} else {
            if (verbose) std::cout << "unexpected lepton " << mc_primaryLepton << " " << mc_targetA << std::endl;
                // unexpected situation.
	}
	double myKEgev = myKE * 0.001;
            // now get the fate from the parameterization code
            // config is 0 (no weight) 1 (weight to noFSI) or 2 (weight to other FSI)
            // nuantinu is 0 = neutrino 1 = antineutrino
	if (verbose) std::cout << " Found QE with fate " << myfate << " " << myKEgev << " A " << mc_targetA << " nu/anu " << isanu << std::endl;
            //weightElasticFSIforQE = GetQEElasticWeight( myfate, myKEgev, mc_targetA, isanu, config);
    }
    
        //if(isaPiPlus && finalpions==0){
    for(int i=0; i<4; i++){
	if(i==piplus14N)break;
	if(piplus14KE[i] <= 0.0)break;
	
            // weight 4.0 nn, 74.0 pn, 22.0 pp to 0 nn 5 pn 95 pp.
            // These are for isoscalar nuclei, replace for different Z/A like Fe, Pb
	if(piplus14fate[i]==51)weightPionAbsorptionBug *= 3.98464/75.8621;
	else if(piplus14fate[i]==50)weightPionAbsorptionBug *= 0.0;
	else if(piplus14fate[i]==52)weightPionAbsorptionBug *= 96.01536/21.2247;
	
            // else found something other than nn, pn, pp
	if (verbose) std::cout << i << " fatepi+ " << piplus14fate[i] << " weight " << weightPionAbsorptionBug << " N " << piplus14N << std::endl;
	
            // alternative, weigh this to some number inspired by LADS
            // D. Rowntree et al PRC 60 (1999) 054610
            // Original GENIE values were from Engel, Mosel arxiv.org nucl-th/9307008
            // A. Engel, W. Cassing, U. Mosel, M. Schafer, Gy. Wolf Nucl.Phys. A572 (1994) p.657
            // another alternative, choose totally different fractions than GENIE
            // like 90% and 100% instead of 96% and treat them as +/- 1 sigma
	
    }
        // do it again for pizero
    for(int i=0; i<4; i++){
	if(i==pizero14N)break;
	if(pizero14KE[i] <= 0.0)break;
	
            // weight 4.0 nn, 74.0 pn, 22.0 pp to 0 nn 5 pn 95 pp.
            // These are for isoscalar nuclei, replace for different Z/A like Fe, Pb
	if(pizero14fate[i]==51)weightPionAbsorptionBug *= 1.0;
	else if(pizero14fate[i]==50)weightPionAbsorptionBug *= 0.12069/0.029132;
	else if(pizero14fate[i]==52)weightPionAbsorptionBug *= 0.12069/0.212247;
            // else found something other than nn, pn, pp
	if (verbose) std::cout << i << " fatepi0 " << pizero14fate[i] << " weight " << weightPionAbsorptionBug << " N " << pizero14N << std::endl;
	
            // alternative, weight 4.0 nn, 74.0 pn, 22.0 pp to LADS
            // another alternative, choose totally different fractions than GENIE
	
    }


        //truth_pionabsorptionbugweight = weightPionAbsorptionBug;	    
    if (verbose) std::cout << "weightPionAbsorptionBug = " << weightPionAbsorptionBug << std::endl;
      
	// print out the entire record to debug.
	//for(Int_t i=0; i < mc_er_nPart; ++i){
	//  std::cout << i << " " << mc_er_ID[i] << " " << mc_er_status[i] << " " << mc_er_FD[i] << " " << mc_er_LD[i] << " " << mc_er_E[i] << std::endl;
	//}
      
    return weightPionAbsorptionBug;

}


double PlotUtils::gethAFSIAbsorptionWeight(PlotUtils::ChainWrapper* chw,
                                           Long64_t entry,
                                           bool verbose)
{

    
    int mc_incoming = chw->GetInt("mc_incoming", entry);
    int mc_primaryLepton = chw->GetInt("mc_primaryLepton", entry);
    int mc_charm = chw->GetInt("mc_charm", entry);
    int mc_intType = chw->GetInt("mc_intType", entry);
    int mc_targetA = chw->GetInt("mc_targetA", entry);
        //int mc_targetZ = chw->GetInt("mc_targetZ", entry);
    int mc_er_nPart = chw->GetInt("mc_er_nPart", entry);
    std::vector<int> mc_er_ID     = chw->GetValueVector<int>("mc_er_ID", entry);
    std::vector<int> mc_er_status = chw->GetValueVector<int>("mc_er_status", entry);
    std::vector<int> mc_er_FD     = chw->GetValueVector<int>("mc_er_FD", entry);
    std::vector<int> mc_er_LD     = chw->GetValueVector<int>("mc_er_LD", entry);
        //std::vector<int> mc_er_mother = chw->GetValueVector<int>("mc_er_mother", entry);
        //std::vector<double> mc_er_Px = chw->GetValueVector<double>("mc_er_Px", entry);
        //std::vector<double> mc_er_Py = chw->GetValueVector<double>("mc_er_Py", entry);
        //std::vector<double> mc_er_Pz = chw->GetValueVector<double>("mc_er_Pz", entry);
    std::vector<double> mc_er_E  = chw->GetValueVector<double>("mc_er_E",  entry);


    
    const Double_t protonmass = 938.272046;
    const Double_t neutronmass = 939.565379;
    const Double_t pionmass = 139.57061;
    const Double_t pizeromass = 134.9770;
    
    Int_t proton14N = 0;  
    Int_t neutron14N = 0;  
    Int_t piplus14N = 0; 
    Int_t piminus14N = 0;  
    Int_t pizero14N = 0;
    
    Double_t proton14KE[4]  = {0.0, 0.0, 0.0, 0.0}; 
    Double_t neutron14KE[4] = {0.0, 0.0, 0.0, 0.0}; 
    Double_t piplus14KE[4]  = {0.0, 0.0, 0.0, 0.0}; 
    Double_t piminus14KE[4] = {0.0, 0.0, 0.0, 0.0}; 
    Double_t pizero14KE[4]  = {0.0, 0.0, 0.0, 0.0}; 
    
    Int_t proton14fate[4]  = {0, 0, 0, 0}; 
    Int_t neutron14fate[4] = {0, 0, 0, 0}; 
    Int_t piplus14fate[4]  = {0, 0, 0, 0}; 
    Int_t piminus14fate[4] = {0, 0, 0, 0}; 
    Int_t pizero14fate[4]  = {0, 0, 0, 0}; 
   
    double weightPionAbsorptionBug = 1.0;
        //double weightElasticFSIforQE = 1.0;

    
      
    for(Int_t i=0; i < mc_er_nPart; ++i){
	if(mc_er_status[i] == 14){
            Int_t tempfate = -1;
            
                // Only track fates for proton, neutron, pizero, piplus, piminus.
            if(mc_er_ID[i] == 2212 || mc_er_ID[i] == 2112 || mc_er_ID[i] == 111 || TMath::Abs(mc_er_ID[i])==211){
                
                    // The distinguishing feature is first daughter last daughter
                    // get an easy handle to these for use later
                Int_t fd = mc_er_FD[i];
                Int_t ld = mc_er_LD[i];

                    // is one of these a pion
                Int_t isaPion = 0;
                if(mc_er_ID[i] == 111 || TMath::Abs(mc_er_ID[i] == 211))isaPion = 1;
                
                    // how many daughters are pions
                Int_t FSpions = 0;
                for(int jj = mc_er_FD[i]; jj <= mc_er_LD[i]; ++jj){
                    if(TMath::Abs(mc_er_ID[jj]) == 211 || mc_er_ID[jj] == 111)FSpions++;
                }
                
                if(mc_er_FD[i] == mc_er_LD[i]){
                        //fate is either 1 = no interaction or 3 = elastic
                        // or 5 = absorption with a 3xxxxxxxx hadblob 
                    
                        // This test of 25.00000 MeV is to separate elastic
                        //TODO replace this with a nucleus specific offset
                        // unknown nuclei get an offset like 8 MeV
                    Double_t offset = 25.000000;
                    if(mc_intType != 1)offset = 0.0;
                    if(mc_er_ID[fd] > 2000000000){
                        tempfate = 5;   // absorption on 3 nucleons
                    } else if(TMath::Abs(mc_er_E[i] - mc_er_E[fd] - offset) < 0.00001){
                        tempfate = 1;   // no scattering
                    } else if(mc_er_ID[fd] != mc_er_ID[i]){
                        tempfate = 2;   // charge exchange ?
                    } else {
                        tempfate = 3;   // elastic fate
                    }
                } else if(FSpions > 0 && !isaPion){
                    tempfate = 8;   // nucleon goes to pions
                } else if(isaPion && FSpions == 0){
                    tempfate = 5;
                        // was this two body absorption?
                    if(ld-fd == 1){
                        int p = 0;
                        int n = 0;
                        if(mc_er_ID[fd] == 2212)p++;
                        if(mc_er_ID[fd] == 2112)n++;
                        if(mc_er_ID[ld] == 2212)p++;
                        if(mc_er_ID[ld] == 2112)n++;
                            // use special codes for producing nn, pn, pp
                        if(p==1 && n==1)tempfate = 51;
                        if(p==2)tempfate = 52;
                        if(n==2)tempfate = 50;
                        
                            //std::cout << "Pion Absorption2? " << mc_er_ID[fd] << " " << mc_er_ID[ld] << " " << fd << " " << ld;
                            //if(tempfate==51)std::cout << " got pn ";
                            //if(tempfate==52)std::cout << " got pp ";
                            //if(tempfate==50)std::cout << " got nn ";
                            //std::cout << std::endl;
                    }
                } else if(isaPion && FSpions >= 2){
                    tempfate = 8;
                } else if(isaPion && FSpions == 1){
                    tempfate = 4;
                    for(int jj = mc_er_FD[i]; jj <= mc_er_LD[i]; ++jj){
                        if((mc_er_ID[jj] == 111 || TMath::Abs(mc_er_ID[jj]) == 211) && mc_er_ID[jj] != mc_er_ID[i]) tempfate = 2;
                    }
                } else {
                    if (verbose) std::cout << " Got unknown fate " << tempfate << std::endl;
                    tempfate = 4;  // or four, can't tell.
                }

                    //std::cout << "fate inttype " << mc_intType << " fd " << mc_er_ID[fd] << " "  << mc_er_E[fd] << " ld " << mc_er_ID[ld] << " " << mc_er_E[ld] << " orig " << mc_er_ID[i] << " " << mc_er_E[i] << " Ediff " << mc_er_E[i]-mc_er_E[fd] << " fate " << tempfate << std::endl;
                    //std::cout << "fate inttype " << mc_intType << " nD " << ld - fd + 1 << " fate " << tempfate <<  "      fd " << mc_er_ID[fd] << " "  << mc_er_E[fd] << " ld " << mc_er_ID[ld] << " " << mc_er_E[ld] << " orig " << mc_er_ID[i] << " " << mc_er_E[i] << " Ediff " << mc_er_E[i]-mc_er_E[fd] << std::endl;
                
                
            }
            
            
                // The next block of code reorders the list of hadrons
                // So that the first four are in decreasing KE order
                // (and more than four drop out of the array, because I'm low-recoiling).
            
            if(mc_er_ID[i] == 2212){
                Double_t myKE = mc_er_E[i] - protonmass;
                    // What was this protons fate ?
                
                
                    // can I get back the FSI code.
                if(proton14N == 0){
                    proton14KE[0] = myKE;
                    proton14fate[0] = tempfate;
                    proton14N++;
                } else if(proton14KE[0] < myKE){
                    proton14KE[3] = proton14KE[2]; proton14KE[2] = proton14KE[1]; proton14KE[1] = proton14KE[0];
                    proton14fate[3] = proton14fate[2]; proton14fate[2] = proton14fate[1]; proton14fate[1] = proton14fate[0];
                    proton14KE[0] = myKE;
                    proton14fate[0] = tempfate;
                    proton14N++;
                } else if(proton14KE[1] < myKE){
                    proton14KE[3] = proton14KE[2]; proton14KE[2] = proton14KE[1];
                    proton14fate[3] = proton14fate[2]; proton14fate[2] = proton14fate[1];
                    proton14KE[1] = myKE;
                    proton14fate[1] = tempfate;
                    proton14N++;
                } else if(proton14KE[2] < myKE){
                    proton14KE[3] = proton14KE[2];
                    proton14fate[3] = proton14fate[2]; 
                    proton14KE[2] = myKE;
                    proton14fate[2] = tempfate;
                    proton14N++;
                } else if(proton14KE[3] < myKE){
                    proton14KE[3] = myKE;
                    proton14fate[3] = tempfate;
                    proton14N++;
                }
                    //std::cout << "proton14 N " << proton14N << " E " << proton14KE[0] << " " << proton14KE[1] << " " << proton14KE[2] << " " << proton14KE[3] << std::endl;
                    // end proton	      
            } else if(mc_er_ID[i] == 2112){
                Double_t myKE = mc_er_E[i] - neutronmass;
                    // can I get back the FSI code.
                if(neutron14N == 0){
                    neutron14KE[0] = myKE;
                    neutron14fate[0] = tempfate;
                    neutron14N++;
                } else if(neutron14KE[0] < myKE){
                    neutron14KE[3] = neutron14KE[2]; neutron14KE[2] = neutron14KE[1]; neutron14KE[1] = neutron14KE[0];
                    neutron14fate[3] = neutron14fate[2]; neutron14fate[2] = neutron14fate[1]; neutron14fate[1] = neutron14fate[0];
                    neutron14KE[0] = myKE;
                    neutron14fate[0] = tempfate;
                    neutron14N++;
                } else if(neutron14KE[1] < myKE){
                    neutron14KE[3] = neutron14KE[2]; neutron14KE[2] = neutron14KE[1];
                    neutron14fate[3] = neutron14fate[2]; neutron14fate[2] = neutron14fate[1];
                    neutron14KE[1] = myKE;
                    neutron14fate[1] = tempfate;
                    neutron14N++;
                } else if(neutron14KE[2] < myKE){
                    neutron14KE[3] = neutron14KE[2];
                    neutron14fate[3] = neutron14fate[2]; 
                    neutron14KE[2] = myKE;
                    neutron14fate[2] = tempfate;
                    neutron14N++;
                } else if(neutron14KE[3] < myKE){
                    neutron14KE[3] = myKE;
                    neutron14fate[3] = tempfate;
                    neutron14N++;
                }
                    //std::cout << "neutron14 N " << neutron14N << " E " << neutron14KE[0] << " " << neutron14KE[1] << " " << neutron14KE[2] << " " << neutron14KE[3] << std::endl;
	  
            } else if(mc_er_ID[i] == 111){
                Double_t myKE = mc_er_E[i] - pizeromass;
                    // can I get back the FSI code.
                if(pizero14N == 0){
                    pizero14KE[0] = myKE;
                    pizero14fate[0] = tempfate;
                    pizero14N++;
                } else if(pizero14KE[0] < myKE){
                    pizero14KE[3] = pizero14KE[2]; pizero14KE[2] = pizero14KE[1]; pizero14KE[1] = pizero14KE[0];
                    pizero14fate[3] = pizero14fate[2]; pizero14fate[2] = pizero14fate[1]; pizero14fate[1] = pizero14fate[0];
                    pizero14KE[0] = myKE;
                    pizero14fate[0] = tempfate;
                    pizero14N++;
                } else if(pizero14KE[1] < myKE){
                    pizero14KE[3] = pizero14KE[2]; pizero14KE[2] = pizero14KE[1];
                    pizero14fate[3] = pizero14fate[2]; pizero14fate[2] = pizero14fate[1];
                    pizero14KE[1] = myKE;
                    pizero14fate[1] = tempfate;
                    pizero14N++;
                } else if(pizero14KE[2] < myKE){
                    pizero14KE[3] = pizero14KE[2];
                    pizero14fate[3] = pizero14fate[2]; 
                    pizero14KE[2] = myKE;
                    pizero14fate[2] = tempfate;
                    pizero14N++;
                } else if(pizero14KE[3] < myKE){
                    pizero14KE[3] = myKE;
                    pizero14fate[3] = tempfate;
                    pizero14N++;
                }
                    //std::cout << "pizero14 N " << pizero14N << " E " << pizero14KE[0] << " " << pizero14KE[1] << " " << pizero14KE[2] << " " << pizero14KE[3] << std::endl;
	  
            } else if(mc_er_ID[i] == 211){
                Double_t myKE = mc_er_E[i] - pionmass;
                    // can I get back the FSI code.
                if(piplus14N == 0){
                    piplus14KE[0] = myKE;
                    piplus14fate[0] = tempfate;
                    piplus14N++;
                } else if(piplus14KE[0] < myKE){
                    piplus14KE[3] = piplus14KE[2]; piplus14KE[2] = piplus14KE[1]; piplus14KE[1] = piplus14KE[0];
                    piplus14fate[3] = piplus14fate[2]; piplus14fate[2] = piplus14fate[1]; piplus14fate[1] = piplus14fate[0];
                    piplus14KE[0] = myKE;
                    piplus14fate[0] = tempfate;
                    piplus14N++;
                } else if(piplus14KE[1] < myKE){
                    piplus14KE[3] = piplus14KE[2]; piplus14KE[2] = piplus14KE[1];
                    piplus14fate[3] = piplus14fate[2]; piplus14fate[2] = piplus14fate[1];
                    piplus14KE[1] = myKE;
                    piplus14fate[1] = tempfate;
                    piplus14N++;
                } else if(piplus14KE[2] < myKE){
                    piplus14KE[3] = piplus14KE[2];
                    piplus14fate[3] = piplus14fate[2]; 
                    piplus14KE[2] = myKE;
                    piplus14fate[2] = tempfate;
                    piplus14N++;
                } else if(piplus14KE[3] < myKE){
                    piplus14KE[3] = myKE;
                    piplus14fate[3] = tempfate;
                    piplus14N++;
                }
                    //std::cout << "pion14 N " << pion14N << " E " << pion14KE[0] << " " << pion14KE[1] << " " << pion14KE[2] << " " << pion14KE[3] << std::endl;
	  
            } else if(mc_er_ID[i] == -211){
                Double_t myKE = mc_er_E[i] - pionmass;
                    // can I get back the FSI code.
                if(piminus14N == 0){
                    piminus14KE[0] = myKE;
                    piminus14fate[0] = tempfate;
                    piminus14N++;
                } else if(piminus14KE[0] < myKE){
                    piminus14KE[3] = piminus14KE[2]; piminus14KE[2] = piminus14KE[1]; piminus14KE[1] = piminus14KE[0];
                    piminus14fate[3] = piminus14fate[2]; piminus14fate[2] = piminus14fate[1]; piminus14fate[1] = piminus14fate[0];
                    piminus14KE[0] = myKE;
                    piminus14fate[0] = tempfate;
                    piminus14N++;
                } else if(piminus14KE[1] < myKE){
                    piminus14KE[3] = piminus14KE[2]; piminus14KE[2] = piminus14KE[1];
                    piminus14fate[3] = piminus14fate[2]; piminus14fate[2] = piminus14fate[1];
                    piminus14KE[1] = myKE;
                    piminus14fate[1] = tempfate;
                    piminus14N++;
                } else if(piminus14KE[2] < myKE){
                    piminus14KE[3] = piminus14KE[2];
                    piminus14fate[3] = piminus14fate[2]; 
                    piminus14KE[2] = myKE;
                    piminus14fate[2] = tempfate;
                    piminus14N++;
                } else if(piminus14KE[3] < myKE){
                    piminus14KE[3] = myKE;
                    piminus14fate[3] = tempfate;
                    piminus14N++;
                }
	    //std::cout << "pion14 N " << pion14N << " E " << pion14KE[0] << " " << pion14KE[1] << " " << pion14KE[2] << " " << pion14KE[3] << std::endl;
                
            }
	}   // end status14
    }
    
        // Now identify what weights we need to use.
    if(mc_intType==1 && !mc_charm){
            // Get the QE elastic fate weight
	bool isanu = false;
	int myfate = -1;
	double myKE = -999.0;
	
	if( mc_incoming == 14 ||  mc_incoming == 12 ||  mc_incoming == 16 ){
            myfate = proton14fate[0];
            myKE = proton14KE[0];
	  
	} else if ( mc_incoming == -14 ||  mc_incoming == -12 ||  mc_incoming == -16 ){
            myfate = neutron14fate[0];
            myKE = neutron14KE[0];
            isanu = true;
            
	} else {
            if (verbose) std::cout << "unexpected lepton " << mc_primaryLepton << " " << mc_targetA << std::endl;
                // unexpected situation.
	}
	double myKEgev = myKE * 0.001;
            // now get the fate from the parameterization code
            // config is 0 (no weight) 1 (weight to noFSI) or 2 (weight to other FSI)
            // nuantinu is 0 = neutrino 1 = antineutrino
	if (verbose) std::cout << " Found QE with fate " << myfate << " " << myKEgev << " A " << mc_targetA << " nu/anu " << isanu << std::endl;
            //weightElasticFSIforQE = GetQEElasticWeight( myfate, myKEgev, mc_targetA, isanu, config);
    }
    
        //if(isaPiPlus && finalpions==0){
    for(int i=0; i<4; i++){
	if(i==piplus14N)break;
	if(piplus14KE[i] <= 0.0)break;
	
            // weight 4.0 nn, 74.0 pn, 22.0 pp to 0 nn 5 pn 95 pp.
            // These are for isoscalar nuclei, replace for different Z/A like Fe, Pb
	if(piplus14fate[i]==51)weightPionAbsorptionBug *= 3.98464/75.8621;
	else if(piplus14fate[i]==50)weightPionAbsorptionBug *= 0.0;
	else if(piplus14fate[i]==52)weightPionAbsorptionBug *= 96.01536/21.2247;
	
            // else found something other than nn, pn, pp
	if (verbose) std::cout << i << " fatepi+ " << piplus14fate[i] << " weight " << weightPionAbsorptionBug << " N " << piplus14N << std::endl;
	
            // alternative, weigh this to some number inspired by LADS
            // D. Rowntree et al PRC 60 (1999) 054610
            // Original GENIE values were from Engel, Mosel arxiv.org nucl-th/9307008
            // A. Engel, W. Cassing, U. Mosel, M. Schafer, Gy. Wolf Nucl.Phys. A572 (1994) p.657
            // another alternative, choose totally different fractions than GENIE
            // like 90% and 100% instead of 96% and treat them as +/- 1 sigma
	
    }
        // do it again for pizero
    for(int i=0; i<4; i++){
	if(i==pizero14N)break;
	if(pizero14KE[i] <= 0.0)break;
	
            // weight 4.0 nn, 74.0 pn, 22.0 pp to 0 nn 5 pn 95 pp.
            // These are for isoscalar nuclei, replace for different Z/A like Fe, Pb
	if(pizero14fate[i]==51)weightPionAbsorptionBug *= 1.0;
	else if(pizero14fate[i]==50)weightPionAbsorptionBug *= 0.12069/0.029132;
	else if(pizero14fate[i]==52)weightPionAbsorptionBug *= 0.12069/0.212247;
            // else found something other than nn, pn, pp
	if (verbose) std::cout << i << " fatepi0 " << pizero14fate[i] << " weight " << weightPionAbsorptionBug << " N " << pizero14N << std::endl;
	
            // alternative, weight 4.0 nn, 74.0 pn, 22.0 pp to LADS
            // another alternative, choose totally different fractions than GENIE
	
    }


        //truth_pionabsorptionbugweight = weightPionAbsorptionBug;	    
    if (verbose) std::cout << "weightPionAbsorptionBug = " << weightPionAbsorptionBug << std::endl;
      
	// print out the entire record to debug.
	//for(Int_t i=0; i < mc_er_nPart; ++i){
	//  std::cout << i << " " << mc_er_ID[i] << " " << mc_er_status[i] << " " << mc_er_FD[i] << " " << mc_er_LD[i] << " " << mc_er_E[i] << std::endl;
	//}
      
    return weightPionAbsorptionBug;


}
