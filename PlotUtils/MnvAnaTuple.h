//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr 12 22:19:11 2012 by ROOT version 5.30/00
// from TTree NukeCC/Tuple created by an AnaTuple managed by AnaTupleManager
// found on file: /minerva/data/users/tice/mc_production/nogrid/minerva/ana/v10r2/00/01/00/00/SIM_minerva_00010000_0001_Ana_Tuple_v3_v10r2.root
//////////////////////////////////////////////////////////

#ifndef MnvAnaTuple_h
#define MnvAnaTuple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream>
#include <string>
#include <vector>
#include <map>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::map;

namespace PlotUtils{
  class MnvAnaTuple{
    public :
      TTree          *fChain;   //!pointer to the analyzed TTree or TChain
      Int_t           fCurrent; //!current Tree number in a TChain

      static const int MAX_N_SLICES = 3; ///< Max number of slices contributing to a PhysicsEvent
      static const int MAX_N_GENIE_WGT_SHIFTS = 10; ///< Max number of nSigma shifts for GENIEWeights (7 is normal)
      static const int MAX_N_FSPART = 25; ///< Max number of final state particles written by GENIE
      static const int MAX_N_ERPART = 50; ///< Max number of particles in GENIE's event record
      static const int MAX_N_WGT_UNIVERSES = 1000; ///< Max number of random universes used for random weights (all use 1000)
      static const int MAX_N_PRONGS = 15; ///< Max number of primary Prongs attached to a PhysicsEvent

      // Declaration of leaf types
      Double_t        eventID;
      Int_t           physEvtNum;
      Int_t           n_hyps;
      Int_t           processType;
      Int_t           primaryPart;
      Int_t           n_slices;
      Int_t           slice_numbers[MAX_N_SLICES];   //[n_slices]
      Int_t           shared_slice;
      Double_t        vtx[4];
      Double_t        vtxErr[4];
      Double_t        E[4];
      Bool_t          found_truth;
      Bool_t          phys_front_activity;
      Bool_t          prim_vtx_has_misassigned_track_direction;
      Bool_t          prim_vtx_has_broken_track;
      Int_t           n_tracks;
      Int_t           n_tracks_non_prim;
      Int_t           n_tracks_prim;
      Int_t           n_tracks_prim_forked;
      Int_t           n_tracks_prim_kinked;
      Int_t           n_vertices_startpoint;
      Int_t           phys_energy_in_road_downstream_nplanes;
      Int_t           phys_energy_in_road_upstream_nplanes;
      Int_t           phys_n_dead_discr_pair;
      Int_t           phys_n_dead_discr_pair_in_prim_track_region;
      Int_t           phys_n_dead_discr_pair_two_mod_downstream_prim_track;
      Int_t           phys_n_dead_discr_pair_two_mod_upstream_prim_vtx;
      Int_t           phys_n_dead_discr_pair_upstream_prim_track_proj;
      Int_t           phys_vertex_is_fiducial;
      Int_t           broken_track_most_us_plane;
      Double_t        prim_vtx_smallest_opening_angle;
      Double_t        energy_from_mc;
      Double_t        energy_from_mc_fraction;
      Double_t        energy_from_mc_fraction_of_highest;
      Double_t        phys_energy_dispersed;
      Double_t        phys_energy_in_road_downstream;
      Double_t        phys_energy_in_road_upstream;
      Double_t        phys_energy_unattached;
      Double_t        primary_track_minerva_energy;
      Double_t        primary_track_minerva_phi;
      Double_t        primary_track_minerva_theta;
      Double_t        primary_track_minerva_end_position[3];
      Double_t        primary_track_minerva_start_position[3];
      Double_t        muon_phi;
      Double_t        muon_theta;
      Double_t        muon_thetaX;
      Double_t        muon_thetaY;
      Double_t        numi_horn_curr;
      Double_t        numi_pot;
      Double_t        numi_x;
      Double_t        numi_x_width;
      Double_t        numi_y;
      Double_t        numi_y_width;
      Bool_t          truth_has_physics_event;
      Int_t           genie_wgt_n_shifts;
      Double_t        truth_genie_wgt_AGKYxF1pi[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_CCQEPauliSupViaFK[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_EtaNCEL[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_FrAbs_N[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_FrAbs_pi[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_FrCEx_N[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_FrCEx_pi[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_FrElas_N[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_FrElas_pi[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_FrInel_N[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_FrInel_pi[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_FrPiProd_N[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_FrPiProd_pi[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_MFP_N[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_MFP_pi[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_MaCCQE[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_MaCCQEshape[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_MaNCEL[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_MaRES[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_MvRES[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_NormCCQE[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_NormCCRES[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_NormDISCC[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_RDecBR1gamma[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_Rvn1pi[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_Rvn2pi[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_Rvp1pi[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_Rvp2pi[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_Theta_Delta2Npi[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_VecFFCCQEshape[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Double_t        truth_genie_wgt_shifts[MAX_N_GENIE_WGT_SHIFTS];   //[genie_wgt_n_shifts]
      Int_t           ev_run;
      Int_t           ev_subrun;
      Int_t           ev_detector;
      Int_t           ev_triggerType;
      Int_t           ev_gate;
      Int_t           ev_global_gate;
      Int_t           ev_gps_time_sec;
      Int_t           ev_gps_time_usec;
      Int_t           mc_run;
      Int_t           mc_subrun;
      Int_t           mc_nInteractions;
      Int_t           mc_MIState;
      Double_t        mc_pot;
      Int_t           mc_beamConfig;
      Int_t           mc_processType;
      Int_t           mc_nthEvtInSpill;
      Int_t           mc_nthEvtInFile;
      Int_t           mc_intType;
      Int_t           mc_current;
      Int_t           mc_charm;
      Double_t        mc_weight;
      Double_t        mc_XSec;
      Double_t        mc_diffXSec;
      Int_t           mc_incoming;
      Double_t        mc_fluxDriverProb;
      Int_t           mc_targetNucleus;
      Int_t           mc_targetZ;
      Int_t           mc_targetA;
      Int_t           mc_targetNucleon;
      Int_t           mc_struckQuark;
      Int_t           mc_seaQuark;
      Int_t           mc_resID;
      Int_t           mc_primaryLepton;
      Double_t        mc_incomingE;
      Double_t        mc_Bjorkenx;
      Double_t        mc_Bjorkeny;
      Double_t        mc_Q2;
      Double_t        mc_nuT;
      Double_t        mc_w;
      Double_t        mc_vtx[4];
      Double_t        mc_incomingPartVec[4];
      Double_t        mc_initNucVec[4];
      Double_t        mc_primFSLepton[4];
      Int_t           mc_nFSPart;
      Double_t        mc_FSPartPx[MAX_N_FSPART];   //[mc_nFSPart]
      Double_t        mc_FSPartPy[MAX_N_FSPART];   //[mc_nFSPart]
      Double_t        mc_FSPartPz[MAX_N_FSPART];   //[mc_nFSPart]
      Double_t        mc_FSPartE[MAX_N_FSPART];   //[mc_nFSPart]
      Int_t           mc_FSPartPDG[MAX_N_FSPART];   //[mc_nFSPart]
      Int_t           mc_er_nPart;
      Int_t           mc_er_ID[MAX_N_ERPART];   //[mc_er_nPart]
      Int_t           mc_er_status[MAX_N_ERPART];   //[mc_er_nPart]
      Double_t        mc_er_posInNucX[MAX_N_ERPART];   //[mc_er_nPart]
      Double_t        mc_er_posInNucY[MAX_N_ERPART];   //[mc_er_nPart]
      Double_t        mc_er_posInNucZ[MAX_N_ERPART];   //[mc_er_nPart]
      Double_t        mc_er_Px[MAX_N_ERPART];   //[mc_er_nPart]
      Double_t        mc_er_Py[MAX_N_ERPART];   //[mc_er_nPart]
      Double_t        mc_er_Pz[MAX_N_ERPART];   //[mc_er_nPart]
      Double_t        mc_er_E[MAX_N_ERPART];   //[mc_er_nPart]
      Int_t           mc_er_FD[MAX_N_ERPART];   //[mc_er_nPart]
      Int_t           mc_er_LD[MAX_N_ERPART];   //[mc_er_nPart]
      Int_t           mc_er_mother[MAX_N_ERPART];   //[mc_er_nPart]
      Int_t           mc_fr_nuParentID;
      Int_t           mc_fr_decMode;
      Double_t        mc_fr_primProtonVtx[3];
      Double_t        mc_fr_primProtonP[4];
      Double_t        mc_fr_nuParentDecVtx[3];
      Double_t        mc_fr_nuParentProdVtx[3];
      Double_t        mc_fr_nuParentProdP[4];
      Double_t        mc_cvweight_total;
      Double_t        wgt;
      Double_t        mc_cvweight_totalFlux;
      Double_t        mc_cvweight_totalXsec;
      Int_t           mc_wgt_GENIE_sz;
      Double_t        mc_wgt_GENIE[MAX_N_WGT_UNIVERSES];   //[mc_wgt_GENIE_sz]
      Int_t           mc_wgt_Flux_Tertiary_sz;
      Double_t        mc_wgt_Flux_Tertiary[MAX_N_WGT_UNIVERSES];   //[mc_wgt_Flux_Tertiary_sz]
      Int_t           mc_wgt_Flux_BeamFocus_sz;
      Double_t        mc_wgt_Flux_BeamFocus[MAX_N_WGT_UNIVERSES];   //[mc_wgt_Flux_BeamFocus_sz]
      Int_t           mc_wgt_Flux_NA49_sz;
      Double_t        mc_wgt_Flux_NA49[MAX_N_WGT_UNIVERSES];   //[mc_wgt_Flux_NA49_sz]
      Int_t           n_prongs;
      Int_t           prong_nParticles[MAX_N_PRONGS];   //[n_prongs]
      Double_t        prong_part_score[MAX_N_PRONGS];   //[n_prongs]
      Double_t        prong_part_mass[MAX_N_PRONGS];   //[n_prongs]
      Int_t           prong_part_charge[MAX_N_PRONGS];   //[n_prongs]
      Int_t           prong_part_pid[MAX_N_PRONGS];   //[n_prongs]
      Double_t        prong_part_E[MAX_N_PRONGS][4];   //[n_prongs]
      Double_t        prong_part_pos[MAX_N_PRONGS][4];   //[n_prongs]

      // List of branches
      TBranch        *b_eventID;   //!
      TBranch        *b_physEvtNum;   //!
      TBranch        *b_n_hyps;   //!
      TBranch        *b_processType;   //!
      TBranch        *b_primaryPart;   //!
      TBranch        *b_n_slices;   //!
      TBranch        *b_slice_numbers;   //!
      TBranch        *b_shared_slice;   //!
      TBranch        *b_vtx;   //!
      TBranch        *b_vtxErr;   //!
      TBranch        *b_E;   //!
      TBranch        *b_found_truth;   //!
      TBranch        *b_phys_front_activity;   //!
      TBranch        *b_prim_vtx_has_misassigned_track_direction;   //!
      TBranch        *b_prim_vtx_has_broken_track;   //!
      TBranch        *b_n_tracks;   //!
      TBranch        *b_n_tracks_non_prim;   //!
      TBranch        *b_n_tracks_prim;   //!
      TBranch        *b_n_tracks_prim_forked;   //!
      TBranch        *b_n_tracks_prim_kinked;   //!
      TBranch        *b_n_vertices_startpoint;   //!
      TBranch        *b_phys_energy_in_road_downstream_nplanes;   //!
      TBranch        *b_phys_energy_in_road_upstream_nplanes;   //!
      TBranch        *b_phys_n_dead_discr_pair;   //!
      TBranch        *b_phys_n_dead_discr_pair_in_prim_track_region;   //!
      TBranch        *b_phys_n_dead_discr_pair_two_mod_downstream_prim_track;   //!
      TBranch        *b_phys_n_dead_discr_pair_two_mod_upstream_prim_vtx;   //!
      TBranch        *b_phys_n_dead_discr_pair_upstream_prim_track_proj;   //!
      TBranch        *b_prim_vtx_smallest_opening_angle; //!
      TBranch        *b_broken_track_most_us_plane; //!
      TBranch        *b_phys_vertex_is_fiducial;   //!
      TBranch        *b_energy_from_mc;   //!
      TBranch        *b_energy_from_mc_fraction;   //!
      TBranch        *b_energy_from_mc_fraction_of_highest;   //!
      TBranch        *b_muon_phi;   //!
      TBranch        *b_muon_theta;   //!
      TBranch        *b_muon_thetaX;   //!
      TBranch        *b_muon_thetaY;   //!
      TBranch        *b_numi_horn_curr;   //!
      TBranch        *b_numi_pot;   //!
      TBranch        *b_numi_x;   //!
      TBranch        *b_numi_x_width;   //!
      TBranch        *b_numi_y;   //!
      TBranch        *b_numi_y_width;   //!
      TBranch        *b_phys_energy_dispersed;   //!
      TBranch        *b_phys_energy_in_road_downstream;   //!
      TBranch        *b_phys_energy_in_road_upstream;   //!
      TBranch        *b_phys_energy_unattached;   //!
      TBranch        *b_primary_track_minerva_energy;   //!
      TBranch        *b_primary_track_minerva_phi;   //!
      TBranch        *b_primary_track_minerva_theta;   //!
      TBranch        *b_primary_track_minerva_end_position;   //!
      TBranch        *b_primary_track_minerva_start_position;   //!
      TBranch        *b_truth_has_physics_event;   //!
      TBranch        *b_genie_wgt_n_shifts;   //!
      TBranch        *b_truth_genie_wgt_AGKYxF1pi;   //!
      TBranch        *b_truth_genie_wgt_CCQEPauliSupViaFK;   //!
      TBranch        *b_truth_genie_wgt_EtaNCEL;   //!
      TBranch        *b_truth_genie_wgt_FrAbs_N;   //!
      TBranch        *b_truth_genie_wgt_FrAbs_pi;   //!
      TBranch        *b_truth_genie_wgt_FrCEx_N;   //!
      TBranch        *b_truth_genie_wgt_FrCEx_pi;   //!
      TBranch        *b_truth_genie_wgt_FrElas_N;   //!
      TBranch        *b_truth_genie_wgt_FrElas_pi;   //!
      TBranch        *b_truth_genie_wgt_FrInel_N;   //!
      TBranch        *b_truth_genie_wgt_FrInel_pi;   //!
      TBranch        *b_truth_genie_wgt_FrPiProd_N;   //!
      TBranch        *b_truth_genie_wgt_FrPiProd_pi;   //!
      TBranch        *b_truth_genie_wgt_MFP_N;   //!
      TBranch        *b_truth_genie_wgt_MFP_pi;   //!
      TBranch        *b_truth_genie_wgt_MaCCQE;   //!
      TBranch        *b_truth_genie_wgt_MaCCQEshape;   //!
      TBranch        *b_truth_genie_wgt_MaNCEL;   //!
      TBranch        *b_truth_genie_wgt_MaRES;   //!
      TBranch        *b_truth_genie_wgt_MvRES;   //!
      TBranch        *b_truth_genie_wgt_NormCCQE;   //!
      TBranch        *b_truth_genie_wgt_NormCCRES;   //!
      TBranch        *b_truth_genie_wgt_NormDISCC;   //!
      TBranch        *b_truth_genie_wgt_RDecBR1gamma;   //!
      TBranch        *b_truth_genie_wgt_Rvn1pi;   //!
      TBranch        *b_truth_genie_wgt_Rvn2pi;   //!
      TBranch        *b_truth_genie_wgt_Rvp1pi;   //!
      TBranch        *b_truth_genie_wgt_Rvp2pi;   //!
      TBranch        *b_truth_genie_wgt_Theta_Delta2Npi;   //!
      TBranch        *b_truth_genie_wgt_VecFFCCQEshape;   //!
      TBranch        *b_truth_genie_wgt_shifts;   //!
      TBranch        *b_ev_run;   //!
      TBranch        *b_ev_subrun;   //!
      TBranch        *b_ev_detector;   //!
      TBranch        *b_ev_triggerType;   //!
      TBranch        *b_ev_gate;   //!
      TBranch        *b_ev_global_gate;   //!
      TBranch        *b_ev_gps_time_sec;   //!
      TBranch        *b_ev_gps_time_usec;   //!
      TBranch        *b_mc_run;   //!
      TBranch        *b_mc_subrun;   //!
      TBranch        *b_mc_nInteractions;   //!
      TBranch        *b_mc_MIState;   //!
      TBranch        *b_mc_pot;   //!
      TBranch        *b_mc_beamConfig;   //!
      TBranch        *b_mc_processType;   //!
      TBranch        *b_mc_nthEvtInSpill;   //!
      TBranch        *b_mc_nthEvtInFile;   //!
      TBranch        *b_mc_intType;   //!
      TBranch        *b_mc_current;   //!
      TBranch        *b_mc_charm;   //!
      TBranch        *b_mc_weight;   //!
      TBranch        *b_mc_XSec;   //!
      TBranch        *b_mc_diffXSec;   //!
      TBranch        *b_mc_incoming;   //!
      TBranch        *b_mc_fluxDriverProb;   //!
      TBranch        *b_mc_targetNucleus;   //!
      TBranch        *b_mc_targetZ;   //!
      TBranch        *b_mc_targetA;   //!
      TBranch        *b_mc_targetNucleon;   //!
      TBranch        *b_mc_struckQuark;   //!
      TBranch        *b_mc_seaQuark;   //!
      TBranch        *b_mc_resID;   //!
      TBranch        *b_mc_primaryLepton;   //!
      TBranch        *b_mc_incomingE;   //!
      TBranch        *b_mc_Bjorkenx;   //!
      TBranch        *b_mc_Bjorkeny;   //!
      TBranch        *b_mc_Q2;   //!
      TBranch        *b_mc_nuT;   //!
      TBranch        *b_mc_w;   //!
      TBranch        *b_mc_vtx;   //!
      TBranch        *b_mc_incomingPartVec;   //!
      TBranch        *b_mc_initNucVec;   //!
      TBranch        *b_mc_primFSLepton;   //!
      TBranch        *b_mc_nFSPart;   //!
      TBranch        *b_mc_FSPartPx;   //!
      TBranch        *b_mc_FSPartPy;   //!
      TBranch        *b_mc_FSPartPz;   //!
      TBranch        *b_mc_FSPartE;   //!
      TBranch        *b_mc_FSPartPDG;   //!
      TBranch        *b_mc_er_nPart;   //!
      TBranch        *b_mc_er_ID;   //!
      TBranch        *b_mc_er_status;   //!
      TBranch        *b_mc_er_posInNucX;   //!
      TBranch        *b_mc_er_posInNucY;   //!
      TBranch        *b_mc_er_posInNucZ;   //!
      TBranch        *b_mc_er_Px;   //!
      TBranch        *b_mc_er_Py;   //!
      TBranch        *b_mc_er_Pz;   //!
      TBranch        *b_mc_er_E;   //!
      TBranch        *b_mc_er_FD;   //!
      TBranch        *b_mc_er_LD;   //!
      TBranch        *b_mc_er_mother;   //!
      TBranch        *b_mc_fr_nuParentID;   //!
      TBranch        *b_mc_fr_decMode;   //!
      TBranch        *b_mc_fr_primProtonVtx;   //!
      TBranch        *b_mc_fr_primProtonP;   //!
      TBranch        *b_mc_fr_nuParentDecVtx;   //!
      TBranch        *b_mc_fr_nuParentProdVtx;   //!
      TBranch        *b_mc_fr_nuParentProdP;   //!
      TBranch        *b_mc_cvweight_total;   //!
      TBranch        *b_wgt;   //!
      TBranch        *b_mc_cvweight_totalFlux;   //!
      TBranch        *b_mc_cvweight_totalXsec;   //!
      TBranch        *b_mc_wgt_GENIE_sz;   //!
      TBranch        *b_mc_wgt_GENIE;   //!
      TBranch        *b_mc_wgt_Flux_Tertiary_sz;   //!
      TBranch        *b_mc_wgt_Flux_Tertiary;   //!
      TBranch        *b_mc_wgt_Flux_BeamFocus_sz;   //!
      TBranch        *b_mc_wgt_Flux_BeamFocus;   //!
      TBranch        *b_mc_wgt_Flux_NA49_sz;   //!
      TBranch        *b_mc_wgt_Flux_NA49;   //!
      TBranch        *b_n_prongs;   //!
      TBranch        *b_prong_nParticles;   //!
      TBranch        *b_prong_part_score;   //!
      TBranch        *b_prong_part_mass;   //!
      TBranch        *b_prong_part_charge;   //!
      TBranch        *b_prong_part_pid;   //!
      TBranch        *b_prong_part_E;   //!
      TBranch        *b_prong_part_pos;   //!

      MnvAnaTuple() {};
      MnvAnaTuple(TTree *tree);
      MnvAnaTuple(TFile *file, const char *treeName);
      MnvAnaTuple(const char *fileName, const char *treeName);

      virtual ~MnvAnaTuple();

      virtual Int_t    GetEntry(Long64_t entry);
      virtual Long64_t LoadTree(Long64_t entry);
      virtual void     Init(TTree *tree);
      virtual Bool_t   Notify();
      virtual void     Show(Long64_t entry = -1);
  };

  MnvAnaTuple::MnvAnaTuple(TTree *tree)
  {
    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.
    if (tree == 0) 
      cout << "ERROR [MnvAnaTuple]: Cannot construct an MnvAnaTuple from a NULL TTree*" << endl;
    else
      Init(tree);
  }

  MnvAnaTuple::MnvAnaTuple(TFile *file, const char *treeName)
  {
    if(!file)
      cout << "ERROR [MnvAnaTuple::MnvAnaTuple]: Cannot look for a TTree in a NULL TFile*" << endl;
    else {
      TTree *tree = (TTree*)file->Get(treeName);
      if( tree )
        Init(tree);
      else
        cout << "ERROR [MnvAnaTuple::MnvAnaTuple]: Cannot find a tree named '" << treeName << "' in TFile named '" << file->GetName() << "'" << endl;
    }
  }

  MnvAnaTuple::MnvAnaTuple(const char *fileName, const char *treeName)
  {
    TFile *file = new TFile(fileName);
    if(!file)
      cout << "ERROR [MnvAnaTuple::MnvAnaTuple]: Cannot open TFile with filename: " << fileName << endl;
    else {
      TTree *tree = (TTree*)file->Get(treeName);
      if( tree )
        Init(tree);
      else
        cout << "ERROR [MnvAnaTuple::MnvAnaTuple]: Cannot find a tree named '" << treeName << "' in file: " << fileName << endl;
    }
  }


  MnvAnaTuple::~MnvAnaTuple()
  {
    if (!fChain) return;
    delete fChain->GetCurrentFile();
  }

  Int_t MnvAnaTuple::GetEntry(Long64_t entry)
  {
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
  }
  Long64_t MnvAnaTuple::LoadTree(Long64_t entry)
  {
    // Set the environment to read one entry
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
    }
    return centry;
  }

  void MnvAnaTuple::Init(TTree *tree)
  {
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("eventID", &eventID, &b_eventID);
    fChain->SetBranchAddress("physEvtNum", &physEvtNum, &b_physEvtNum);
    fChain->SetBranchAddress("n_hyps", &n_hyps, &b_n_hyps);
    fChain->SetBranchAddress("processType", &processType, &b_processType);
    fChain->SetBranchAddress("primaryPart", &primaryPart, &b_primaryPart);
    fChain->SetBranchAddress("n_slices", &n_slices, &b_n_slices);
    fChain->SetBranchAddress("slice_numbers", slice_numbers, &b_slice_numbers);
    fChain->SetBranchAddress("shared_slice", &shared_slice, &b_shared_slice);
    fChain->SetBranchAddress("vtx", vtx, &b_vtx);
    fChain->SetBranchAddress("vtxErr", vtxErr, &b_vtxErr);
    fChain->SetBranchAddress("E", E, &b_E);
    fChain->SetBranchAddress("found_truth", &found_truth, &b_found_truth);
    fChain->SetBranchAddress("phys_front_activity", &phys_front_activity, &b_phys_front_activity);
    fChain->SetBranchAddress("prim_vtx_has_misassigned_track_direction", &prim_vtx_has_misassigned_track_direction, &b_prim_vtx_has_misassigned_track_direction);
    fChain->SetBranchAddress("prim_vtx_has_broken_track", &prim_vtx_has_broken_track, &b_prim_vtx_has_broken_track);
    fChain->SetBranchAddress("n_tracks", &n_tracks, &b_n_tracks);
    fChain->SetBranchAddress("n_tracks_non_prim", &n_tracks_non_prim, &b_n_tracks_non_prim);
    fChain->SetBranchAddress("n_tracks_prim", &n_tracks_prim, &b_n_tracks_prim);
    fChain->SetBranchAddress("n_tracks_prim_forked", &n_tracks_prim_forked, &b_n_tracks_prim_forked);
    fChain->SetBranchAddress("n_tracks_prim_kinked", &n_tracks_prim_kinked, &b_n_tracks_prim_kinked);
    fChain->SetBranchAddress("n_vertices_startpoint", &n_vertices_startpoint, &b_n_vertices_startpoint);
    fChain->SetBranchAddress("phys_energy_in_road_downstream_nplanes", &phys_energy_in_road_downstream_nplanes, &b_phys_energy_in_road_downstream_nplanes);
    fChain->SetBranchAddress("phys_energy_in_road_upstream_nplanes", &phys_energy_in_road_upstream_nplanes, &b_phys_energy_in_road_upstream_nplanes);
    fChain->SetBranchAddress("phys_n_dead_discr_pair", &phys_n_dead_discr_pair, &b_phys_n_dead_discr_pair); 
    fChain->SetBranchAddress("phys_n_dead_discr_pair_in_prim_track_region", &phys_n_dead_discr_pair_in_prim_track_region, &b_phys_n_dead_discr_pair_in_prim_track_region);
    fChain->SetBranchAddress("phys_n_dead_discr_pair_two_mod_downstream_prim_track", &phys_n_dead_discr_pair_two_mod_downstream_prim_track, &b_phys_n_dead_discr_pair_two_mod_downstream_prim_track);
    fChain->SetBranchAddress("phys_n_dead_discr_pair_two_mod_upstream_prim_vtx", &phys_n_dead_discr_pair_two_mod_upstream_prim_vtx, &b_phys_n_dead_discr_pair_two_mod_upstream_prim_vtx);
    fChain->SetBranchAddress("phys_n_dead_discr_pair_upstream_prim_track_proj", &phys_n_dead_discr_pair_upstream_prim_track_proj, &b_phys_n_dead_discr_pair_upstream_prim_track_proj);
    fChain->SetBranchAddress("prim_vtx_smallest_opening_angle", &prim_vtx_smallest_opening_angle, &b_prim_vtx_smallest_opening_angle);
    fChain->SetBranchAddress("broken_track_most_us_plane", &broken_track_most_us_plane, &b_broken_track_most_us_plane);
    fChain->SetBranchAddress("phys_vertex_is_fiducial", &phys_vertex_is_fiducial, &b_phys_vertex_is_fiducial);
    fChain->SetBranchAddress("energy_from_mc", &energy_from_mc, &b_energy_from_mc);
    fChain->SetBranchAddress("energy_from_mc_fraction", &energy_from_mc_fraction, &b_energy_from_mc_fraction);
    fChain->SetBranchAddress("energy_from_mc_fraction_of_highest", &energy_from_mc_fraction_of_highest, &b_energy_from_mc_fraction_of_highest);
    fChain->SetBranchAddress("muon_phi", &muon_phi, &b_muon_phi);
    fChain->SetBranchAddress("muon_theta", &muon_theta, &b_muon_theta);
    fChain->SetBranchAddress("muon_thetaX", &muon_thetaX, &b_muon_thetaX);
    fChain->SetBranchAddress("muon_thetaY", &muon_thetaY, &b_muon_thetaY);
    fChain->SetBranchAddress("numi_horn_curr", &numi_horn_curr, &b_numi_horn_curr);
    fChain->SetBranchAddress("numi_pot", &numi_pot, &b_numi_pot);
    fChain->SetBranchAddress("numi_x", &numi_x, &b_numi_x);
    fChain->SetBranchAddress("numi_x_width", &numi_x_width, &b_numi_x_width);
    fChain->SetBranchAddress("numi_y", &numi_y, &b_numi_y);
    fChain->SetBranchAddress("numi_y_width", &numi_y_width, &b_numi_y_width);
    fChain->SetBranchAddress("phys_energy_dispersed", &phys_energy_dispersed, &b_phys_energy_dispersed);
    fChain->SetBranchAddress("phys_energy_in_road_downstream", &phys_energy_in_road_downstream, &b_phys_energy_in_road_downstream);
    fChain->SetBranchAddress("phys_energy_in_road_upstream", &phys_energy_in_road_upstream, &b_phys_energy_in_road_upstream);
    fChain->SetBranchAddress("phys_energy_unattached", &phys_energy_unattached, &b_phys_energy_unattached);
    fChain->SetBranchAddress("primary_track_minerva_energy", &primary_track_minerva_energy, &b_primary_track_minerva_energy);
    fChain->SetBranchAddress("primary_track_minerva_phi", &primary_track_minerva_phi, &b_primary_track_minerva_phi);
    fChain->SetBranchAddress("primary_track_minerva_theta", &primary_track_minerva_theta, &b_primary_track_minerva_theta);
    fChain->SetBranchAddress("primary_track_minerva_end_position", primary_track_minerva_end_position, &b_primary_track_minerva_end_position);
    fChain->SetBranchAddress("primary_track_minerva_start_position", primary_track_minerva_start_position, &b_primary_track_minerva_start_position);
    fChain->SetBranchAddress("truth_has_physics_event", &truth_has_physics_event, &b_truth_has_physics_event);
    fChain->SetBranchAddress("genie_wgt_n_shifts", &genie_wgt_n_shifts, &b_genie_wgt_n_shifts);
    fChain->SetBranchAddress("truth_genie_wgt_AGKYxF1pi", &truth_genie_wgt_AGKYxF1pi, &b_truth_genie_wgt_AGKYxF1pi);
    fChain->SetBranchAddress("truth_genie_wgt_CCQEPauliSupViaFK", &truth_genie_wgt_CCQEPauliSupViaFK, &b_truth_genie_wgt_CCQEPauliSupViaFK);
    fChain->SetBranchAddress("truth_genie_wgt_EtaNCEL", &truth_genie_wgt_EtaNCEL, &b_truth_genie_wgt_EtaNCEL);
    fChain->SetBranchAddress("truth_genie_wgt_FrAbs_N", &truth_genie_wgt_FrAbs_N, &b_truth_genie_wgt_FrAbs_N);
    fChain->SetBranchAddress("truth_genie_wgt_FrAbs_pi", &truth_genie_wgt_FrAbs_pi, &b_truth_genie_wgt_FrAbs_pi);
    fChain->SetBranchAddress("truth_genie_wgt_FrCEx_N", &truth_genie_wgt_FrCEx_N, &b_truth_genie_wgt_FrCEx_N);
    fChain->SetBranchAddress("truth_genie_wgt_FrCEx_pi", &truth_genie_wgt_FrCEx_pi, &b_truth_genie_wgt_FrCEx_pi);
    fChain->SetBranchAddress("truth_genie_wgt_FrElas_N", &truth_genie_wgt_FrElas_N, &b_truth_genie_wgt_FrElas_N);
    fChain->SetBranchAddress("truth_genie_wgt_FrElas_pi", &truth_genie_wgt_FrElas_pi, &b_truth_genie_wgt_FrElas_pi);
    fChain->SetBranchAddress("truth_genie_wgt_FrInel_N", &truth_genie_wgt_FrInel_N, &b_truth_genie_wgt_FrInel_N);
    fChain->SetBranchAddress("truth_genie_wgt_FrInel_pi", &truth_genie_wgt_FrInel_pi, &b_truth_genie_wgt_FrInel_pi);
    fChain->SetBranchAddress("truth_genie_wgt_FrPiProd_N", &truth_genie_wgt_FrPiProd_N, &b_truth_genie_wgt_FrPiProd_N);
    fChain->SetBranchAddress("truth_genie_wgt_FrPiProd_pi", &truth_genie_wgt_FrPiProd_pi, &b_truth_genie_wgt_FrPiProd_pi);
    fChain->SetBranchAddress("truth_genie_wgt_MFP_N", &truth_genie_wgt_MFP_N, &b_truth_genie_wgt_MFP_N);
    fChain->SetBranchAddress("truth_genie_wgt_MFP_pi", &truth_genie_wgt_MFP_pi, &b_truth_genie_wgt_MFP_pi);
    fChain->SetBranchAddress("truth_genie_wgt_MaCCQE", &truth_genie_wgt_MaCCQE, &b_truth_genie_wgt_MaCCQE);
    fChain->SetBranchAddress("truth_genie_wgt_MaCCQEshape", &truth_genie_wgt_MaCCQEshape, &b_truth_genie_wgt_MaCCQEshape);
    fChain->SetBranchAddress("truth_genie_wgt_MaNCEL", &truth_genie_wgt_MaNCEL, &b_truth_genie_wgt_MaNCEL);
    fChain->SetBranchAddress("truth_genie_wgt_MaRES", &truth_genie_wgt_MaRES, &b_truth_genie_wgt_MaRES);
    fChain->SetBranchAddress("truth_genie_wgt_MvRES", &truth_genie_wgt_MvRES, &b_truth_genie_wgt_MvRES);
    fChain->SetBranchAddress("truth_genie_wgt_NormCCQE", &truth_genie_wgt_NormCCQE, &b_truth_genie_wgt_NormCCQE);
    fChain->SetBranchAddress("truth_genie_wgt_NormCCRES", &truth_genie_wgt_NormCCRES, &b_truth_genie_wgt_NormCCRES);
    fChain->SetBranchAddress("truth_genie_wgt_NormDISCC", &truth_genie_wgt_NormDISCC, &b_truth_genie_wgt_NormDISCC);
    fChain->SetBranchAddress("truth_genie_wgt_RDecBR1gamma", &truth_genie_wgt_RDecBR1gamma, &b_truth_genie_wgt_RDecBR1gamma);
    fChain->SetBranchAddress("truth_genie_wgt_Rvn1pi", &truth_genie_wgt_Rvn1pi, &b_truth_genie_wgt_Rvn1pi);
    fChain->SetBranchAddress("truth_genie_wgt_Rvn2pi", &truth_genie_wgt_Rvn2pi, &b_truth_genie_wgt_Rvn2pi);
    fChain->SetBranchAddress("truth_genie_wgt_Rvp1pi", &truth_genie_wgt_Rvp1pi, &b_truth_genie_wgt_Rvp1pi);
    fChain->SetBranchAddress("truth_genie_wgt_Rvp2pi", &truth_genie_wgt_Rvp2pi, &b_truth_genie_wgt_Rvp2pi);
    fChain->SetBranchAddress("truth_genie_wgt_Theta_Delta2Npi", &truth_genie_wgt_Theta_Delta2Npi, &b_truth_genie_wgt_Theta_Delta2Npi);
    fChain->SetBranchAddress("truth_genie_wgt_VecFFCCQEshape", &truth_genie_wgt_VecFFCCQEshape, &b_truth_genie_wgt_VecFFCCQEshape);
    fChain->SetBranchAddress("truth_genie_wgt_shifts", &truth_genie_wgt_shifts, &b_truth_genie_wgt_shifts);
    fChain->SetBranchAddress("ev_run", &ev_run, &b_ev_run);
    fChain->SetBranchAddress("ev_subrun", &ev_subrun, &b_ev_subrun);
    fChain->SetBranchAddress("ev_detector", &ev_detector, &b_ev_detector);
    fChain->SetBranchAddress("ev_triggerType", &ev_triggerType, &b_ev_triggerType);
    fChain->SetBranchAddress("ev_gate", &ev_gate, &b_ev_gate);
    fChain->SetBranchAddress("ev_global_gate", &ev_global_gate, &b_ev_global_gate);
    fChain->SetBranchAddress("ev_gps_time_sec", &ev_gps_time_sec, &b_ev_gps_time_sec);
    fChain->SetBranchAddress("ev_gps_time_usec", &ev_gps_time_usec, &b_ev_gps_time_usec);
    fChain->SetBranchAddress("mc_run", &mc_run, &b_mc_run);
    fChain->SetBranchAddress("mc_subrun", &mc_subrun, &b_mc_subrun);
    fChain->SetBranchAddress("mc_nInteractions", &mc_nInteractions, &b_mc_nInteractions);
    fChain->SetBranchAddress("mc_MIState", &mc_MIState, &b_mc_MIState);
    fChain->SetBranchAddress("mc_pot", &mc_pot, &b_mc_pot);
    fChain->SetBranchAddress("mc_beamConfig", &mc_beamConfig, &b_mc_beamConfig);
    fChain->SetBranchAddress("mc_processType", &mc_processType, &b_mc_processType);
    fChain->SetBranchAddress("mc_nthEvtInSpill", &mc_nthEvtInSpill, &b_mc_nthEvtInSpill);
    fChain->SetBranchAddress("mc_nthEvtInFile", &mc_nthEvtInFile, &b_mc_nthEvtInFile);
    fChain->SetBranchAddress("mc_intType", &mc_intType, &b_mc_intType);
    fChain->SetBranchAddress("mc_current", &mc_current, &b_mc_current);
    fChain->SetBranchAddress("mc_charm", &mc_charm, &b_mc_charm);
    fChain->SetBranchAddress("mc_weight", &mc_weight, &b_mc_weight);
    fChain->SetBranchAddress("mc_XSec", &mc_XSec, &b_mc_XSec);
    fChain->SetBranchAddress("mc_diffXSec", &mc_diffXSec, &b_mc_diffXSec);
    fChain->SetBranchAddress("mc_incoming", &mc_incoming, &b_mc_incoming);
    fChain->SetBranchAddress("mc_fluxDriverProb", &mc_fluxDriverProb, &b_mc_fluxDriverProb);
    fChain->SetBranchAddress("mc_targetNucleus", &mc_targetNucleus, &b_mc_targetNucleus);
    fChain->SetBranchAddress("mc_targetZ", &mc_targetZ, &b_mc_targetZ);
    fChain->SetBranchAddress("mc_targetA", &mc_targetA, &b_mc_targetA);
    fChain->SetBranchAddress("mc_targetNucleon", &mc_targetNucleon, &b_mc_targetNucleon);
    fChain->SetBranchAddress("mc_struckQuark", &mc_struckQuark, &b_mc_struckQuark);
    fChain->SetBranchAddress("mc_seaQuark", &mc_seaQuark, &b_mc_seaQuark);
    fChain->SetBranchAddress("mc_resID", &mc_resID, &b_mc_resID);
    fChain->SetBranchAddress("mc_primaryLepton", &mc_primaryLepton, &b_mc_primaryLepton);
    fChain->SetBranchAddress("mc_incomingE", &mc_incomingE, &b_mc_incomingE);
    fChain->SetBranchAddress("mc_Bjorkenx", &mc_Bjorkenx, &b_mc_Bjorkenx);
    fChain->SetBranchAddress("mc_Bjorkeny", &mc_Bjorkeny, &b_mc_Bjorkeny);
    fChain->SetBranchAddress("mc_Q2", &mc_Q2, &b_mc_Q2);
    fChain->SetBranchAddress("mc_nuT", &mc_nuT, &b_mc_nuT);
    fChain->SetBranchAddress("mc_w", &mc_w, &b_mc_w);
    fChain->SetBranchAddress("mc_vtx", mc_vtx, &b_mc_vtx);
    fChain->SetBranchAddress("mc_incomingPartVec", mc_incomingPartVec, &b_mc_incomingPartVec);
    fChain->SetBranchAddress("mc_initNucVec", mc_initNucVec, &b_mc_initNucVec);
    fChain->SetBranchAddress("mc_primFSLepton", mc_primFSLepton, &b_mc_primFSLepton);
    fChain->SetBranchAddress("mc_nFSPart", &mc_nFSPart, &b_mc_nFSPart);
    fChain->SetBranchAddress("mc_FSPartPx", mc_FSPartPx, &b_mc_FSPartPx);
    fChain->SetBranchAddress("mc_FSPartPy", mc_FSPartPy, &b_mc_FSPartPy);
    fChain->SetBranchAddress("mc_FSPartPz", mc_FSPartPz, &b_mc_FSPartPz);
    fChain->SetBranchAddress("mc_FSPartE", mc_FSPartE, &b_mc_FSPartE);
    fChain->SetBranchAddress("mc_FSPartPDG", mc_FSPartPDG, &b_mc_FSPartPDG);
    fChain->SetBranchAddress("mc_er_nPart", &mc_er_nPart, &b_mc_er_nPart);
    fChain->SetBranchAddress("mc_er_ID", mc_er_ID, &b_mc_er_ID);
    fChain->SetBranchAddress("mc_er_status", mc_er_status, &b_mc_er_status);
    fChain->SetBranchAddress("mc_er_posInNucX", mc_er_posInNucX, &b_mc_er_posInNucX);
    fChain->SetBranchAddress("mc_er_posInNucY", mc_er_posInNucY, &b_mc_er_posInNucY);
    fChain->SetBranchAddress("mc_er_posInNucZ", mc_er_posInNucZ, &b_mc_er_posInNucZ);
    fChain->SetBranchAddress("mc_er_Px", mc_er_Px, &b_mc_er_Px);
    fChain->SetBranchAddress("mc_er_Py", mc_er_Py, &b_mc_er_Py);
    fChain->SetBranchAddress("mc_er_Pz", mc_er_Pz, &b_mc_er_Pz);
    fChain->SetBranchAddress("mc_er_E", mc_er_E, &b_mc_er_E);
    fChain->SetBranchAddress("mc_er_FD", mc_er_FD, &b_mc_er_FD);
    fChain->SetBranchAddress("mc_er_LD", mc_er_LD, &b_mc_er_LD);
    fChain->SetBranchAddress("mc_er_mother", mc_er_mother, &b_mc_er_mother);
    fChain->SetBranchAddress("mc_fr_nuParentID", &mc_fr_nuParentID, &b_mc_fr_nuParentID);
    fChain->SetBranchAddress("mc_fr_decMode", &mc_fr_decMode, &b_mc_fr_decMode);
    fChain->SetBranchAddress("mc_fr_primProtonVtx", mc_fr_primProtonVtx, &b_mc_fr_primProtonVtx);
    fChain->SetBranchAddress("mc_fr_primProtonP", mc_fr_primProtonP, &b_mc_fr_primProtonP);
    fChain->SetBranchAddress("mc_fr_nuParentDecVtx", mc_fr_nuParentDecVtx, &b_mc_fr_nuParentDecVtx);
    fChain->SetBranchAddress("mc_fr_nuParentProdVtx", mc_fr_nuParentProdVtx, &b_mc_fr_nuParentProdVtx);
    fChain->SetBranchAddress("mc_fr_nuParentProdP", mc_fr_nuParentProdP, &b_mc_fr_nuParentProdP);
    fChain->SetBranchAddress("mc_cvweight_total", &mc_cvweight_total, &b_mc_cvweight_total);
    fChain->SetBranchAddress("wgt", &wgt, &b_wgt);
    fChain->SetBranchAddress("mc_cvweight_totalFlux", &mc_cvweight_totalFlux, &b_mc_cvweight_totalFlux);
    fChain->SetBranchAddress("mc_cvweight_totalXsec", &mc_cvweight_totalXsec, &b_mc_cvweight_totalXsec);
    fChain->SetBranchAddress("mc_wgt_GENIE_sz", &mc_wgt_GENIE_sz, &b_mc_wgt_GENIE_sz);
    fChain->SetBranchAddress("mc_wgt_GENIE", mc_wgt_GENIE, &b_mc_wgt_GENIE);
    fChain->SetBranchAddress("mc_wgt_Flux_Tertiary_sz", &mc_wgt_Flux_Tertiary_sz, &b_mc_wgt_Flux_Tertiary_sz);
    fChain->SetBranchAddress("mc_wgt_Flux_Tertiary", mc_wgt_Flux_Tertiary, &b_mc_wgt_Flux_Tertiary);
    fChain->SetBranchAddress("mc_wgt_Flux_BeamFocus_sz", &mc_wgt_Flux_BeamFocus_sz, &b_mc_wgt_Flux_BeamFocus_sz);
    fChain->SetBranchAddress("mc_wgt_Flux_BeamFocus", mc_wgt_Flux_BeamFocus, &b_mc_wgt_Flux_BeamFocus);
    fChain->SetBranchAddress("mc_wgt_Flux_NA49_sz", &mc_wgt_Flux_NA49_sz, &b_mc_wgt_Flux_NA49_sz);
    fChain->SetBranchAddress("mc_wgt_Flux_NA49", mc_wgt_Flux_NA49, &b_mc_wgt_Flux_NA49);
    fChain->SetBranchAddress("n_prongs", &n_prongs, &b_n_prongs);
    fChain->SetBranchAddress("prong_nParticles", prong_nParticles, &b_prong_nParticles);
    fChain->SetBranchAddress("prong_part_score", prong_part_score, &b_prong_part_score);
    fChain->SetBranchAddress("prong_part_mass", prong_part_mass, &b_prong_part_mass);
    fChain->SetBranchAddress("prong_part_charge", prong_part_charge, &b_prong_part_charge);
    fChain->SetBranchAddress("prong_part_pid", prong_part_pid, &b_prong_part_pid);
    fChain->SetBranchAddress("prong_part_E", prong_part_E, &b_prong_part_E);
    fChain->SetBranchAddress("prong_part_pos", prong_part_pos, &b_prong_part_pos);
    Notify();
  }

  Bool_t MnvAnaTuple::Notify()
  {
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
  }

  void MnvAnaTuple::Show(Long64_t entry)
  {
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
  }

}//end PlotUtils namespace
#endif
