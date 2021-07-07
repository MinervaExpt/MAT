from collections import OrderedDict

error_bands = OrderedDict()
error_bands["Flux"]                   = ["Flux"]
error_bands["Muon"]                   = ["Muon_Energy_MINERvA",
                                         "Muon_Energy_MINOS",
                                         "Muon_Energy_Resolution",
                                         "BeamAngleX",
                                         "BeamAngleY"]
error_bands["RPA"]                    = ["RPA_HighQ2",
                                         "RPA_LowQ2"]
error_bands["Low_recoil_2p2h_Tune"]   = ["Low_Recoil_2p2h_Tune"]
error_bands["Response"]               = ["Proton_Response",
                                         "Pion_Response",
                                         "EM_Response",
                                         "Other_Response"]
error_bands["MINOS"]                  = ["MINOS_Reconstruction_Efficiency"]
error_bands["GENIE"]                  = ["GENIE_FrAbs_N",
                                         "GENIE_FrAbs_N",
                                         "GENIE_FrAbs_pi",
                                         "GENIE_FrCEx_N",
                                         "GENIE_FrCEx_pi",
                                         "GENIE_FrElas_N",
                                         "GENIE_FrElas_pi",
                                         "GENIE_FrInel_N",
                                         "GENIE_FrInel_pi",
                                         "GENIE_FrPiProd_N",
                                         "GENIE_FrPiProd_pi",
                                         "GENIE_MFP_N",
                                         "GENIE_MFP_pi",
                                         "GENIE_AGKYxF1pi",
                                         "GENIE_AhtBY",
                                         "GENIE_BhtBY",
                                         "GENIE_CCQEPauliSupViaKF",
                                         "GENIE_CV1uBY",
                                         "GENIE_CV2uBY",
                                         "GENIE_EtaNCEL",
                                         "GENIE_MaCCQE",
                                         "GENIE_MaCCQEshape",
                                         "GENIE_MaNCEL",
                                         "GENIE_MaRES",
                                         "GENIE_MvRES",
                                         "GENIE_NormCCQE",
                                         "GENIE_NormCCRES",
                                         "GENIE_NormDISCC",
                                         "GENIE_NormNCRES",
                                         "GENIE_RDecBR1gamma",
                                         "GENIE_Rvn1pi",
                                         "GENIE_Rvn2pi",
                                         "GENIE_Rvn3pi",
                                         "GENIE_Rvp1pi",
                                         "GENIE_Rvp2pi",
                                         "GENIE_Theta_Delta2Npi",
                                         "GENIE_VecFFCCQEshape"]

