 
#include "weight_fsi_cai.h"

using namespace PlotUtils;



void weight_fsi_cai::read(const TString  f)
//Read in the params doubles from a file
//argument: valid filename
{
  if(f!="")
  {
    fFSIWeight = TFile::Open(f,"READONLY");
    if (fFSIWeight){
      hElaFSIWeight = (TH1D*) fFSIWeight->Get("weight_elastic_fsi");
      hNoFSIWeight  = (TH1D*) fFSIWeight->Get("weight_no_fsi");
      std::cout << "have read in weights from file " << f <<std::endl;
    }
    else{
      //File could not be read
      std::cout << "File could not be read" << std::endl;
    }
  }
}

double weight_fsi_cai::getGenieBE(const int A, const int Z)
{
  if (A == 1) return 0; // hydrogen
  if (A == 6) return 17.0; // lithium
  if (A == 12) return 25.0; // carbon
  if (A == 16) return 27.0; // oxygen
  if (A == 24) return 32.0; // magnesium
  if (A == 56) return 36.0; // 56 iron
  if (A == 58) return 36.0; // 58 nickel
  if (A == 208) return 44.0; // 208 lead
  return 0;
}

double weight_fsi_cai::getWeightInternal(double Ti, int mode )
{
  if (mode == 1) 
  {
    if( Ti > hNoFSIWeight->GetXaxis()->GetXmax() ) Ti = hNoFSIWeight->GetXaxis()->GetXmax();
    else if( Ti < hNoFSIWeight->GetXaxis()->GetXmin() ) Ti = hNoFSIWeight->GetXaxis()->GetXmin();
    return hNoFSIWeight->Interpolate(Ti);
  }
  if (mode == 3)
  { 
    if( Ti > hElaFSIWeight->GetXaxis()->GetXmax() ) Ti = hElaFSIWeight->GetXaxis()->GetXmax();
    if( Ti < hElaFSIWeight->GetXaxis()->GetXmin() ) Ti = hElaFSIWeight->GetXaxis()->GetXmin();
    return hElaFSIWeight->Interpolate(Ti);
  }
  return 1;
}

int weight_fsi_cai::getQEFSIMode( int mc_intType, int mc_targetA, int mc_targetZ, int mc_er_nPart, const int* mc_er_ID, const int* mc_er_status, const int* mc_er_FD, const int* mc_er_LD, const int* mc_er_mother, const double* mc_er_Px, const double* mc_er_Py, const double* mc_er_Pz, const double* mc_er_E, double &Ti)
{
  if (mc_intType != 1) return -999; //not QE
  if (mc_targetA == 1) return -999; // no FSI

  genie_particle fsi_part, prefsi_part;
  int fsi_npart = 0;
  int prefsi_npart = 0;

  for ( int i = 0; i< mc_er_nPart; ++ i )
  {
    //cout<<"iPart = "<<i<<endl;
    const int pdg         = mc_er_ID[i];
    const int status      = mc_er_status[i];
    const int fd          = mc_er_FD[i];
    const int ld          = mc_er_LD[i];
    const int fm          = mc_er_mother[i];
    const double px       = mc_er_Px[i];
    const double py       = mc_er_Py[i];
    const double pz       = mc_er_Pz[i];
    const double E        = mc_er_E[i];


    TLorentzVector p4(px,py,pz,E);
    //double m = p4.M();
    //double p = p4.P();
    //cout<<"Defined m, p"<<endl;

    //if (abs(pdg) == 14 || pdg == 22) m = 0;

    genie_particle pobj(i,pdg,status,fd,ld,fm);
    pobj.SetMomentumEnergy( p4 );
    //cout<<"pobj defined"<<endl;
    
    if ( pdg == m_FSPDG && status == 1) {fsi_part = pobj; ++fsi_npart;}
    if ( pdg == m_FSPDG && status == 14) {prefsi_part = pobj; ++prefsi_npart;}

    //cout<<"inserting "<<i<<"th part"<<endl;
    //particle_map.insert( std::make_pair( i, *pobj ) );
  }
  //cout<<"end loop part"<<endl;
  //assert( !fsi_part.default_particle() );
  Ti = prefsi_part.T();
  if (prefsi_part.default_particle() ) return 11; // expected proton
  if (fsi_part.default_particle() ) return 12; // expected proton
  if (fsi_part.GetMother() != prefsi_part.GetId() ) return 13;
  if (prefsi_part.GetFD() != prefsi_part.GetLD() ) return 14;

  if (fsi_npart !=1 && prefsi_npart !=1)
  {
    return -1; //shouldn't exist
  }

  double Ei = prefsi_part.p4().E();
  double Ep = Ei - getGenieBE( mc_targetA, mc_targetZ );

  double pprime = sqrt( Ep*Ep - prefsi_part.p4().M2() );

  double pold = prefsi_part.p4().P();

  if (pold == 0 ) return 15;

  double scale = pprime/pold;

  double pxp = scale* prefsi_part.p4().Px();
  double pyp = scale* prefsi_part.p4().Py();
  double pzp = scale* prefsi_part.p4().Pz();

  TLorentzVector p4p( pxp, pyp,pzp,Ep);
  
  genie_particle be_part = prefsi_part;
  be_part.SetMomentumEnergy( p4p );

  //part = fsi_part.p4();
  //bePart = be_part.p4();

  bool noint = abs( be_part.p4().E() - fsi_part.p4().E() )< 0.0001 && 
                fsi_part.p4().Angle( prefsi_part.p4().Vect() )<0.01;
  bool elastic = abs( be_part.p4().E() - fsi_part.p4().E() )< 5 ;
  if (noint) return 1;
  if (elastic) return 3;
  return 5;
}


int weight_fsi_cai::getResFSIMode( int mc_intType, int mc_targetA, int mc_targetZ, int mc_er_nPart, const int* mc_er_ID, const int* mc_er_status, const int* mc_er_FD, const int* mc_er_LD, const int* mc_er_mother, const double* mc_er_Px, const double* mc_er_Py, const double* mc_er_Pz, const double* mc_er_E, double &Ti)
{
  if (mc_intType ==4 ) return 0;
  if (mc_targetA == 1) return 0;

  std::map<int, genie_particle> particle_map;
  genie_particle nofsi_pi;
  genie_particle fsi_pi;
  int fsi_npi = 0;
  for (int i = 0; i < mc_er_nPart; ++i) 
  {
    const int pdg       = mc_er_ID[i];
    const int status    = mc_er_status[i];
    const int fd        = mc_er_FD[i];
    const int ld        = mc_er_LD[i];
    const int fm        = mc_er_mother[i];
    const double px = mc_er_Px[i];
    const double py = mc_er_Py[i];
    const double pz = mc_er_Pz[i];
    const double e  = mc_er_E[i];
    
    //const double p = sqrt(px * px + py * py + pz * pz);
    //double m = sqrt(e * e - p * p);

    //if (pdg == -14 || pdg == 22) m = 0.0;
    
    genie_particle pobj(i,pdg,status,fd,ld,fm);
    TLorentzVector p4(px,py,pz,e);
    p4.RotateX(-theta_b);
    pobj.SetMomentumEnergy(p4);
    
    if (pdg == -211 && status ==  1) { fsi_pi = pobj; ++fsi_npi; }
    if (pdg == -211 && status == 14) nofsi_pi = pobj;

    
    particle_map.insert(std::make_pair(i,pobj));

        //bool __verbose = false;
        //if (__verbose) {
        //    std::cout.setf(std::ios_base::fixed);
        //    std::cout.precision(2);
        //    std::cout << setw(8) << i << setw(15) << pdg << setw(8) << status
        //              << setw(8) << fm << "   ("
        //              << setw(10) << px << setw(10) << py
        //              << setw(10) << pz << setw(10) << e-m << "  )"
        //              << std::endl;;
        //}
  }

  if (fsi_npi < 1) { std::cout << "warning: no final state pi-" << std::endl; return 0; }

  assert(!fsi_pi.default_particle());

  if (nofsi_pi.default_particle()) return 0; // pi- not come from pi- as in cex or pi prod

  std::vector<genie_particle> pim_ancestry;
  if (!fsi_pi.default_particle()) 
  {
    genie_particle last_particle = fsi_pi;
    pim_ancestry.push_back(last_particle);
    for (;;) 
    {
      int mother_id = last_particle.GetMother();
      if (mother_id < 1) break;
      genie_particle mother = particle_map[mother_id];
      pim_ancestry.push_back(mother);
      last_particle = mother;
    }
  }

  if (pim_ancestry[1] != nofsi_pi) return 0;
    
  const double delta_px = fsi_pi.p4().Px() - nofsi_pi.p4().Px();
  const double delta_py = fsi_pi.p4().Px() - nofsi_pi.p4().Px();
  const double delta_pz = fsi_pi.p4().Px() - nofsi_pi.p4().Px();
  const double delta_k  = fsi_pi.T() - nofsi_pi.T();

  Ti =  nofsi_pi.T();
 
  bool noint = delta_px == 0.0 && delta_py == 0.0 && delta_pz == 0.0;
  bool elastic   = (!noint) && abs(delta_k)  < 10.0;
  bool inelastic = (!noint) && abs(delta_k) >= 10.0; 
    
  if (noint) return 1;
  else if (elastic) return 3;
  else if (inelastic) return 4;
  else return 10;
}

double weight_fsi_cai::getWeight( int mc_intType, int mc_targetA, int mc_targetZ, 
          int mc_er_nPart, const int* mc_er_ID, const int* mc_er_status, 
          const int* mc_er_FD, const int* mc_er_LD, const int* mc_er_mother, 
          const double* mc_er_Px, const double* mc_er_Py, const double* mc_er_Pz,
          const double* mc_er_E)
{
  double Ti = 0;
  int Mode = (m_reweighQE)? 
    getQEFSIMode( mc_intType, mc_targetA, mc_targetZ, mc_er_nPart,
      mc_er_ID, mc_er_status, mc_er_FD, mc_er_LD, mc_er_mother,
      mc_er_Px, mc_er_Py, mc_er_Pz, mc_er_E, Ti ): 
    getResFSIMode( mc_intType, mc_targetA, mc_targetZ, mc_er_nPart,
      mc_er_ID, mc_er_status, mc_er_FD, mc_er_LD, mc_er_mother,
      mc_er_Px, mc_er_Py, mc_er_Pz, mc_er_E, Ti );
  Ti/=1000; // convert to GeV
  return getWeightInternal( Ti, Mode );
}




double weight_fsi_cai::getWeightPoly( int mc_intType, int mc_targetA, int mc_targetZ, 
          int mc_er_nPart, int* mc_er_ID, int* mc_er_status, 
          int* mc_er_FD, int* mc_er_LD, int* mc_er_mother, 
          double* mc_er_Px, double* mc_er_Py, double* mc_er_Pz,
          double* mc_er_E)
{
  if (mc_intType != 1 || mc_targetA==1 ) return 1.;
  double Ti = 0;
  int Mode = (m_reweighQE)? 
    getQEFSIMode( mc_intType, mc_targetA, mc_targetZ, mc_er_nPart,
      mc_er_ID, mc_er_status, mc_er_FD, mc_er_LD, mc_er_mother,
      mc_er_Px, mc_er_Py, mc_er_Pz, mc_er_E, Ti ): 
    getResFSIMode( mc_intType, mc_targetA, mc_targetZ, mc_er_nPart,
      mc_er_ID, mc_er_status, mc_er_FD, mc_er_LD, mc_er_mother,
      mc_er_Px, mc_er_Py, mc_er_Pz, mc_er_E, Ti );
  Ti/=1000; // convert to GeV
  return getWeightInternalRik( Ti,mc_targetA, Mode, m_config );
}





double weight_fsi_cai::getWeightInternalRik( double inputKEgev, int inputA, int mode, int config )
{
  //if (config == 0 or config == 3 or config>4) return 1.;
  if (config == 0) return 1;
  if (mode <= 0) return 1;
  if (config == 1 || config == 2)
  {

    if(inputKEgev < getGenieBE(inputA)/1000.+0.010 )return 1.0;

    if ( mode == 3 ) return 0.;

    //config == 2
    //double lowEpolyC12[10] = {17.8228, 1551.6, -170931, 6.48098e6, -1.34102e8, 1.6961e9, -1.35026e10, 6.6204e10, -1.82973e11, 2.18416e11};
    //double midEpolyC12[10] = {4.17489, -45.544, 345.733, -1485.0, 3875.45, -6356.52, 6588.9, -4194.6, 1499.09, -230.484};
    //double highEpolyC12[10] = {1.23147, -0.0527092, 0.0081572, 0.0, 0.0,   0.0, 0.0, 0.0, 0.0, 0.0};

    //double lowEpolyO16[10] = {-11.4283, 4987.76, -340194., 1.1199e7, -2.16573e8, 2.63567e9, -2.04893e10, 9.89449e10, -2.70812e11, 3.21319e11};
    //double midEpolyO16[10] = {4.85878, -55.0283, 418.162, -1798.99, 4700.67, -7716.98, 8001.46, -5089.56, 1814.09, -277.456};
    //double highEpolyO16[10] = {1.26765, -0.0583833, 0.00853185, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0};

    //double lowEpolyAr40[10] = {-12.4288, 46744.8, -2.47696e6, 5.40691e7, -5.84637e8, 2.5746e9, 7.87892e9, -1.43636e11, 6.19253e11, -9.44823e11};
    //double midEpolyAr40[10] = {31.0973, -517.199, 4089.58, -18404.7, 51320.3, -91701.5, 105134., -74779.8, 3032.9, -5203.64};
    //double highEpolyAr40[10] = {1.47462, -0.111806, 0.0166705, 0.0, 0.0,   0.0, 0.0, 0.0, 0.0, 0.0};

    //double lowEpolyFe56[10] = {2533.41, -140121., 3.27464e6, -4.09799e7, 2.88262e8, -1.07774e9, 1.67053e9, 0.0, 0.0, 0.0};
    //double midEpolyFe56[10] = {31.9074, -525.4, 4131.33, -18520.4, 51485.8, -91755.4, 104940., -74465.9, 29837.6, -5157.84};
    //double highEpolyFe56[10] = {1.52685, -0.127437, 0.0198632, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0};

    //double lowEpolyPb208[10] = {6301.34, -347755., 8.07485e6, -1.00173e8, 6.97538e8, -2.57953e9, 3.95315e9, 0.0, 0.0, 0.0};
    //double midEpolyPb208[10] = {61.6161, -1022.96, 7948.2, -35212.3, 96827.2, -170874., 193711., -136368., 54248.1, -9316.09};
    //double highEpolyPb208[10] = {1.90105, -0.195386, 0.0286129, 0.0, 0.0,   0.0, 0.0, 0.0, 0.0, 0.0};

    ////======================= antinu
    //double lowEpolyC12Anu[10] = {91.8662, -3841.37, 17803.8, 2.46643e6, -7.77780e7, 1.16001e9, -1.00604e10, 5.18932e10, -1.48160e11, 1.80790e11};
    //double midEpolyC12Anu[10] = {6.87059, -87.0001, 658.091, -2865.93, 7713.87, -13257.9, 14592.2, -9963.09, 3845.83, -641.742};
    //double highEpolyC12Anu[10] = {1.19153, -0.0273976, 0.00320112,  0.0, 0.0,   0.0, 0.0, 0.0, 0.0, 0.0};




    ////config == 4
    ////updated weights from Rik
    //double lowEpolyC12w4[10] = {-1.09786, 353.897, -14438.2, 249803.0, -961139.0, -3.30078e7, 5.9223e8, -4.42502e9, 1.62258e10, -2.39127e10}; 
    //double midEpolyC12w4[10] = {1.40574, 0.580308, 33.3332, -323.69, 1330.14, -3060.15, 4217.66, -3458.10, 1555.36, -295.373};

    //double highEpolyC12w4[10] = {1.16048, -0.00252581, 0.000616833, 0.0, 0.0,   0.0, 0.0, 0.0, 0.0, 0.0};

    //
    //double lowEpolyO16w4[10] = { -1.60048, 435.953, -19971.4, 452149.0, -5.40787e6, 2.83386e7, 5.7011e7, -1.55825e9, 7.6214e9, -1.28315e10};
    //double midEpolyO16w4[10] = { 1.44861, -0.40647, 42.8448, -374.366, 1495.64, -3402.08, 4662.46, -3809.16, 1708.30, -323.51};
    //double highEpolyO16w4[10] = {1.15438, 0.00410489, -0.000841064, 0.0, 0.0,   0.0, 0.0, 0.0, 0.0, 0.0};

    //double lowEpolyAr40w4[10] = {-15.9188, 2230.47, -116236.0, 3.35816e6, -5.99261e7, 6.88940e8, -5.12413e9, 2.38571e10, -6.32693e10, 7.29876e10}; 
    //double midEpolyAr40w4[10] =  {1.37013, 1.51099, 23.8883, -274.205, 1180.48, -2784.74, 3904.85, -3244.2, 1474.71, -282.505};
    //double highEpolyAr40w4[10] = {1.15711, 0.00116488, -0.000304764, 0.0, 0.0,   0.0, 0.0, 0.0, 0.0, 0.0};

    //
    //double lowEpolyFe56w4[10] = {-17.5586, 2416.34, -125628.0, 3.63612e6, -6.5208e7, 7.55278e8, -5.67114e9, 2.66995e10, -7.16941e10, 8.38294e10};
    //double midEpolyFe56w4[10] = {1.36436, 1.61799, 22.7834, -266.814, 1149.24, -2702.28, 3770.63, -3114.18, 1406.30, -267.505};
    //double highEpolyFe56w4[10] = {1.1583, -0.000677155, 0.000283527, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0};

    //double lowEpolyPb208w4[10] = {-22.7764, 3056.80, -158958.0, 4.6053e6, -8.2631e7, 9.56844e8, -7.17712e9, 3.37296e10, -9.03555e10, 1.05349e11};
    //double midEpolyPb208w4[10] = {1.38084, 1.34013, 24.5788, -272.457, 1158.25, -2708.59, 3770.54, -3112.01, 1405.81, -267.677};
    //double highEpolyPb208w4[10] = {1.15639, 0.00157778, -0.000370309, 0.0, 0.0,   0.0, 0.0, 0.0, 0.0, 0.0};

 
    // QE has a 25 MeV offset for this GENIE version
    // but the weights are really high, so move this forward by 10 MeV.
    //if(inputKEgev < 0.025)return 1.0;
    
    // For protons from Delta, same or 25 MeV offset?

    double *poly;
    double *low, *mid, *high;

    //Helium 4
    if( inputA >3 && inputA<=10 )
    {
      if(config==2){
        // weight elasticFSI to otherFSI
        low = lowEpolyHe4w2; mid = midEpolyHe4w2; high = highEpolyHe4w2;
      } else {
        // weight elasticFSI to noFSI
        low = lowEpolyHe4; mid = midEpolyHe4; high = highEpolyHe4;
        if( ! m_neutrinoMode ) low = lowEpolyHe4Anu; mid = midEpolyHe4Anu; high = highEpolyHe4Anu;
      }
    } else if(inputA > 10 && inputA <= 15){
      // For carbon, and maybe nitrogen
      if(config==2){
        // weight elasticFSI to otherFSI
        low = lowEpolyC12w2; mid = midEpolyC12w2; high = highEpolyC12w2;
      } else {
        // weight elasticFSI to noFSI
        low = lowEpolyC12; mid = midEpolyC12; high = highEpolyC12;
        if( ! m_neutrinoMode ) low = lowEpolyC12Anu; mid = midEpolyC12Anu; high = highEpolyC12Anu;
      }
    } else if(inputA > 15 && inputA <= 23){
      // mostly for oxygen
      if(config==2){
        low = lowEpolyO16w2; mid = midEpolyO16w2; high = highEpolyO16w2;
      } else {
        low = lowEpolyO16; mid = midEpolyO16; high = highEpolyO16;
        if( ! m_neutrinoMode ) low = lowEpolyO16Anu; mid = midEpolyO16Anu; high = highEpolyO16Anu;
        //weight elasticFSI to noFSI
      }
    } else if (inputA  > 23 && inputA < 32 ){
      if(config==2){
        // weight elasticFSI to otherFSI
        low = lowEpolySi28w2; mid = midEpolySi28w2; high = highEpolySi28w2;
      } else {
        //weight elasticFSI to noFSI
        low = lowEpolySi28; mid = midEpolySi28; high = highEpolySi28;
        if( ! m_neutrinoMode ) low = lowEpolySi28Anu; mid = midEpolySi28Anu; high = highEpolySi28Anu;
      }
    } else if(inputA > 32 && inputA <= 48){
      if(config==2){
        // weight elasticFSI to otherFSI
        low = lowEpolyAr40w2; mid = midEpolyAr40w2; high = highEpolyAr40w2;
      } else {
        //weight elasticFSI to noFSI
        low = lowEpolyAr40; mid = midEpolyAr40; high = highEpolyAr40;
        if( ! m_neutrinoMode ) low = lowEpolyAr40Anu; mid = midEpolyAr40Anu; high = highEpolyAr40Anu;
      }
    } else if(inputA > 48 && inputA <= 70){
      // Use iron for nickel and such
      if(config==2){
        // weight elasticFSI to otherFSI
        low = lowEpolyFe56w2; mid = midEpolyFe56w2; high = highEpolyFe56w2;
      } else {
        //weight elasticFSI to noFSI
        low = lowEpolyFe56; mid = midEpolyFe56; high = highEpolyFe56;
        if( ! m_neutrinoMode ) low = lowEpolyFe56Anu; mid = midEpolyFe56Anu; high = highEpolyFe56Anu;
      }
    } else if(inputA > 195 && inputA <= 240){
      if(config==2){
        // weight elasticFSI to otherFSI
        low = lowEpolyPb208w2; mid = midEpolyPb208w2; high = highEpolyPb208w2;
      } else {
        //weight elasticFSI to noFSI
        low = lowEpolyPb208; mid = midEpolyPb208; high = highEpolyPb208;
        if( ! m_neutrinoMode ) low = lowEpolyPb208Anu; mid = midEpolyPb208Anu; high = highEpolyPb208Anu;
      }
    } else {
      // some nuclei are not tested.
      return 1.0;
    }
    poly = polySwitcher( inputKEgev, low, mid, high );

    // Finally calculate the weight. 
    double powke = 1.0;
    double tempweight = 0.0;
    for(int i=0; i<10; i++){
      tempweight += powke * poly[i];
      powke *= inputKEgev;
    }

    // debug print statement
    //std::cout << "ElasticFSIWeight inputKEgev " << inputKEgev << " inputA " << inputA << " config " << config << " tempweight " << tempweight << std::endl;


    if (config == 4)
    {
      if (mode == 1) return 1;
      if (mode == 3) return 0;
      if (mode >= 4) return tempweight;
    }
    if (config == 2)
    {
      if (mode >= 4) return 1;
      if (mode == 3) return 0;
      if (mode == 1) return tempweight;
    }
  }
  return 1.;
}

double* weight_fsi_cai::polySwitcher(double inputKEgev, double *low, double *mid, double *high )
{
  if(inputKEgev > m_highcutoff)inputKEgev = m_highcutoff;
  double *poly;
  poly = high;
  if(inputKEgev < m_breakpointlow)poly = low;
  else if(inputKEgev < m_breakpointhigh)poly = mid;
  return poly;
}

