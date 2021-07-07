#include "MnvNormalization.h"

#include <iostream>
#include <cstdlib>
#include <cmath>

#include <TString.h>
#include <TError.h>


using namespace PlotUtils;
using namespace MnvNorm;

//================================

MnvNormalizer::MnvNormalizer(const std::string& processing_name, const std::string& playlist_name )
    : __the_correction_struct(0)
{
    ReadCorrectionTable();
    m_processing = processing_name;
    LoadPlaylist( playlist_name );
}

//================================
// Destructor
//================================
MnvNormalizer::~MnvNormalizer() {};

void MnvNormalizer::ReadCorrectionTable()
{
    
    correction_struct* correction_data1 = new correction_struct;
    correction_data1->mnv_mu_reco_eff             = resurrection::minerva1::mnv_mu_reco_eff;
    correction_data1->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - resurrection::minerva1::mnv_mu_reco_eff); 
    correction_data1->minos_mu_reco_eff_lowp      = resurrection::minerva1::minos_mu_reco_eff_lowp;
    correction_data1->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - resurrection::minerva1::minos_mu_reco_eff_lowp);
    correction_data1->minos_mu_reco_eff_highp     = resurrection::minerva1::minos_mu_reco_eff_highp;
    correction_data1->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - resurrection::minerva1::minos_mu_reco_eff_highp);

    __correction_table["Resurrection"]["minerva1"] = correction_data1;

    correction_struct* correction_data1e = new correction_struct;
    correction_data1e->mnv_mu_reco_eff             = eroica::minerva1::mnv_mu_reco_eff;
    correction_data1e->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minerva1::mnv_mu_reco_eff);
    correction_data1e->minos_mu_reco_eff_lowp      = eroica::minerva1::minos_mu_reco_eff_lowp;
    correction_data1e->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minerva1::minos_mu_reco_eff_lowp);
    correction_data1e->minos_mu_reco_eff_highp     = eroica::minerva1::minos_mu_reco_eff_highp;
    correction_data1e->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minerva1::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minerva1"] = correction_data1e;

    correction_struct* correction_data1low = new correction_struct;
    correction_data1low->mnv_mu_reco_eff             = resurrection::minerva1low::mnv_mu_reco_eff;
    correction_data1low->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - resurrection::minerva1low::mnv_mu_reco_eff);
    correction_data1low->minos_mu_reco_eff_lowp      = resurrection::minerva1low::minos_mu_reco_eff_lowp;
    correction_data1low->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - resurrection::minerva1low::minos_mu_reco_eff_lowp);
    correction_data1low->minos_mu_reco_eff_highp     = resurrection::minerva1low::minos_mu_reco_eff_highp;
    correction_data1low->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - resurrection::minerva1low::minos_mu_reco_eff_highp);

    __correction_table["Resurrection"]["minerva1low"] = correction_data1low;

    correction_struct* correction_data1lowe = new correction_struct;
    correction_data1lowe->mnv_mu_reco_eff             = eroica::minerva1low::mnv_mu_reco_eff;
    correction_data1lowe->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minerva1low::mnv_mu_reco_eff);
    correction_data1lowe->minos_mu_reco_eff_lowp      = eroica::minerva1low::minos_mu_reco_eff_lowp;
    correction_data1lowe->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minerva1low::minos_mu_reco_eff_lowp);
    correction_data1lowe->minos_mu_reco_eff_highp     = eroica::minerva1low::minos_mu_reco_eff_highp;
    correction_data1lowe->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minerva1low::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minerva1low"] = correction_data1lowe;
    
    correction_struct* correction_data2 = new correction_struct;
    correction_data2->mnv_mu_reco_eff             = resurrection::minerva2::mnv_mu_reco_eff;
    correction_data2->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - resurrection::minerva2::mnv_mu_reco_eff);
    correction_data2->minos_mu_reco_eff_lowp      = resurrection::minerva2::minos_mu_reco_eff_lowp;
    correction_data2->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - resurrection::minerva2::minos_mu_reco_eff_lowp);
    correction_data2->minos_mu_reco_eff_highp     = resurrection::minerva2::minos_mu_reco_eff_highp;
    correction_data2->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - resurrection::minerva2::minos_mu_reco_eff_highp);

    __correction_table["Resurrection"]["minerva2"] = correction_data2;

    correction_struct* correction_data2e = new correction_struct;
    correction_data2e->mnv_mu_reco_eff             = eroica::minerva2::mnv_mu_reco_eff;
    correction_data2e->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minerva2::mnv_mu_reco_eff);
    correction_data2e->minos_mu_reco_eff_lowp      = eroica::minerva2::minos_mu_reco_eff_lowp;
    correction_data2e->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minerva2::minos_mu_reco_eff_lowp);
    correction_data2e->minos_mu_reco_eff_highp     = eroica::minerva2::minos_mu_reco_eff_highp;
    correction_data2e->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minerva2::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minerva2"] = correction_data2e;

    correction_struct* correction_data3 = new correction_struct;
    correction_data3->mnv_mu_reco_eff             = resurrection::minerva3::mnv_mu_reco_eff;
    correction_data3->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - resurrection::minerva3::mnv_mu_reco_eff);
    correction_data3->minos_mu_reco_eff_lowp      = resurrection::minerva3::minos_mu_reco_eff_lowp;
    correction_data3->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - resurrection::minerva3::minos_mu_reco_eff_lowp);
    correction_data3->minos_mu_reco_eff_highp     = resurrection::minerva3::minos_mu_reco_eff_highp;
    correction_data3->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - resurrection::minerva3::minos_mu_reco_eff_highp);

    __correction_table["Resurrection"]["minerva3"] = correction_data3;

    correction_struct* correction_data3e = new correction_struct;
    correction_data3e->mnv_mu_reco_eff             = eroica::minerva3::mnv_mu_reco_eff;
    correction_data3e->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minerva3::mnv_mu_reco_eff);
    correction_data3e->minos_mu_reco_eff_lowp      = eroica::minerva3::minos_mu_reco_eff_lowp;
    correction_data3e->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minerva3::minos_mu_reco_eff_lowp);
    correction_data3e->minos_mu_reco_eff_highp     = eroica::minerva3::minos_mu_reco_eff_highp;
    correction_data3e->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minerva3::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minerva3"] = correction_data3e;

    correction_struct* correction_data4 = new correction_struct;
    correction_data4->mnv_mu_reco_eff             = resurrection::minerva4::mnv_mu_reco_eff;
    correction_data4->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - resurrection::minerva4::mnv_mu_reco_eff);
    correction_data4->minos_mu_reco_eff_lowp      = resurrection::minerva4::minos_mu_reco_eff_lowp;
    correction_data4->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - resurrection::minerva4::minos_mu_reco_eff_lowp);
    correction_data4->minos_mu_reco_eff_highp     = resurrection::minerva4::minos_mu_reco_eff_highp;
    correction_data4->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - resurrection::minerva4::minos_mu_reco_eff_highp);

    __correction_table["Resurrection"]["minerva4"] = correction_data4;

    correction_struct* correction_data4e = new correction_struct;
    correction_data4e->mnv_mu_reco_eff             = eroica::minerva4::mnv_mu_reco_eff;
    correction_data4e->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minerva4::mnv_mu_reco_eff);
    correction_data4e->minos_mu_reco_eff_lowp      = eroica::minerva4::minos_mu_reco_eff_lowp;
    correction_data4e->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minerva4::minos_mu_reco_eff_lowp);
    correction_data4e->minos_mu_reco_eff_highp     = eroica::minerva4::minos_mu_reco_eff_highp;
    correction_data4e->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minerva4::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minerva4"] = correction_data4e;

    correction_struct* correction_data5 = new correction_struct;
    correction_data5->mnv_mu_reco_eff             = resurrection::minerva5::mnv_mu_reco_eff;
    correction_data5->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - resurrection::minerva5::mnv_mu_reco_eff);
    correction_data5->minos_mu_reco_eff_lowp      = resurrection::minerva5::minos_mu_reco_eff_lowp;
    correction_data5->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - resurrection::minerva5::minos_mu_reco_eff_lowp);
    correction_data5->minos_mu_reco_eff_highp     = resurrection::minerva5::minos_mu_reco_eff_highp;
    correction_data5->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - resurrection::minerva5::minos_mu_reco_eff_highp);

    __correction_table["Resurrection"]["minerva5"] = correction_data5;
    
    correction_struct* correction_data5e = new correction_struct;
    correction_data5e->mnv_mu_reco_eff             = eroica::minerva5::mnv_mu_reco_eff;
    correction_data5e->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minerva5::mnv_mu_reco_eff);
    correction_data5e->minos_mu_reco_eff_lowp      = eroica::minerva5::minos_mu_reco_eff_lowp;
    correction_data5e->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minerva5::minos_mu_reco_eff_lowp);
    correction_data5e->minos_mu_reco_eff_highp     = eroica::minerva5::minos_mu_reco_eff_highp;
    correction_data5e->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minerva5::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minerva5"] = correction_data5e;
    
    correction_struct* correction_data6 = new correction_struct;
    correction_data6->mnv_mu_reco_eff             = resurrection::minerva6::mnv_mu_reco_eff;
    correction_data6->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - resurrection::minerva6::mnv_mu_reco_eff);
    correction_data6->minos_mu_reco_eff_lowp      = resurrection::minerva6::minos_mu_reco_eff_lowp;
    correction_data6->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - resurrection::minerva6::minos_mu_reco_eff_lowp);
    correction_data6->minos_mu_reco_eff_highp     = resurrection::minerva6::minos_mu_reco_eff_highp;
    correction_data6->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - resurrection::minerva6::minos_mu_reco_eff_highp);

    __correction_table["Resurrection"]["minerva6"] = correction_data6;

    correction_struct* correction_data6e = new correction_struct;
    correction_data6e->mnv_mu_reco_eff             = eroica::minerva6::mnv_mu_reco_eff;
    correction_data6e->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minerva6::mnv_mu_reco_eff); 
    correction_data6e->minos_mu_reco_eff_lowp      = eroica::minerva6::minos_mu_reco_eff_lowp;
    correction_data6e->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minerva6::minos_mu_reco_eff_lowp);
    correction_data6e->minos_mu_reco_eff_highp     = eroica::minerva6::minos_mu_reco_eff_highp;
    correction_data6e->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minerva6::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minerva6"] = correction_data6e;

    correction_struct* correction_data7 = new correction_struct;
    correction_data7->mnv_mu_reco_eff             = resurrection::minerva7::mnv_mu_reco_eff;
    correction_data7->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - resurrection::minerva7::mnv_mu_reco_eff);
    correction_data7->minos_mu_reco_eff_lowp      = resurrection::minerva7::minos_mu_reco_eff_lowp;
    correction_data7->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - resurrection::minerva7::minos_mu_reco_eff_lowp);
    correction_data7->minos_mu_reco_eff_highp     = resurrection::minerva7::minos_mu_reco_eff_highp;
    correction_data7->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - resurrection::minerva7::minos_mu_reco_eff_highp);

    __correction_table["Resurrection"]["minerva7"] = correction_data7;

    correction_struct* correction_data7e = new correction_struct;
    correction_data7e->mnv_mu_reco_eff             = eroica::minerva7::mnv_mu_reco_eff;
    correction_data7e->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minerva7::mnv_mu_reco_eff);
    correction_data7e->minos_mu_reco_eff_lowp      = eroica::minerva7::minos_mu_reco_eff_lowp;
    correction_data7e->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minerva7::minos_mu_reco_eff_lowp);
    correction_data7e->minos_mu_reco_eff_highp     = eroica::minerva7::minos_mu_reco_eff_highp;
    correction_data7e->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minerva7::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minerva7"] = correction_data7e;

    correction_struct* correction_data8 = new correction_struct;
    correction_data8->mnv_mu_reco_eff             = resurrection::minerva8::mnv_mu_reco_eff;
    correction_data8->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - resurrection::minerva8::mnv_mu_reco_eff);
    correction_data8->minos_mu_reco_eff_lowp      = resurrection::minerva8::minos_mu_reco_eff_lowp;
    correction_data8->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - resurrection::minerva8::minos_mu_reco_eff_lowp);
    correction_data8->minos_mu_reco_eff_highp     = resurrection::minerva8::minos_mu_reco_eff_highp;

    __correction_table["Resurrection"]["minerva8"] = correction_data8;

    correction_struct* correction_data8e = new correction_struct;
    correction_data8e->mnv_mu_reco_eff             = eroica::minerva8::mnv_mu_reco_eff;
    correction_data8e->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minerva8::mnv_mu_reco_eff);
    correction_data8e->minos_mu_reco_eff_lowp      = eroica::minerva8::minos_mu_reco_eff_lowp;
    correction_data8e->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minerva8::minos_mu_reco_eff_lowp);
    correction_data8e->minos_mu_reco_eff_highp     = eroica::minerva8::minos_mu_reco_eff_highp;
    correction_data8e->minos_mu_reco_eff_highp_err =0.5 * fabs(1.0 - eroica::minerva8::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minerva8e"] = correction_data8e;

    correction_struct* correction_data9 = new correction_struct;
    correction_data9->mnv_mu_reco_eff             = resurrection::minerva9::mnv_mu_reco_eff;
    correction_data9->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - resurrection::minerva9::mnv_mu_reco_eff);
    correction_data9->minos_mu_reco_eff_lowp      = resurrection::minerva9::minos_mu_reco_eff_lowp;
    correction_data9->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 -  resurrection::minerva9::minos_mu_reco_eff_lowp);
    correction_data9->minos_mu_reco_eff_highp     = resurrection::minerva9::minos_mu_reco_eff_highp;
    correction_data9->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - resurrection::minerva9::minos_mu_reco_eff_highp);

    __correction_table["Resurrection"]["minerva9"] = correction_data9;

    correction_struct* correction_data9e = new correction_struct;
    correction_data9e->mnv_mu_reco_eff             = eroica::minerva9::mnv_mu_reco_eff;
    correction_data9e->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minerva9::mnv_mu_reco_eff);
    correction_data9e->minos_mu_reco_eff_lowp      = eroica::minerva9::minos_mu_reco_eff_lowp;
    correction_data9e->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minerva9::minos_mu_reco_eff_lowp);
    correction_data9e->minos_mu_reco_eff_highp     = eroica::minerva9::minos_mu_reco_eff_highp;
    correction_data9e->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minerva9::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minerva9"] = correction_data9e;

    correction_struct* correction_data10 = new correction_struct;
    correction_data10->mnv_mu_reco_eff             = resurrection::minerva10::mnv_mu_reco_eff;
    correction_data10->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - resurrection::minerva10::mnv_mu_reco_eff);
    correction_data10->minos_mu_reco_eff_lowp      = resurrection::minerva10::minos_mu_reco_eff_lowp;
    correction_data10->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - resurrection::minerva10::minos_mu_reco_eff_lowp);
    correction_data10->minos_mu_reco_eff_highp     = resurrection::minerva10::minos_mu_reco_eff_highp;
    correction_data10->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 -  resurrection::minerva10::minos_mu_reco_eff_highp);

    __correction_table["Resurrection"]["minerva10"] = correction_data10;


    correction_struct* correction_data10e = new correction_struct;
    correction_data10e->mnv_mu_reco_eff             = eroica::minerva10::mnv_mu_reco_eff;
    correction_data10e->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minerva10::mnv_mu_reco_eff);
    correction_data10e->minos_mu_reco_eff_lowp      = eroica::minerva10::minos_mu_reco_eff_lowp;
    correction_data10e->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minerva10::minos_mu_reco_eff_lowp);
    correction_data10e->minos_mu_reco_eff_highp     = eroica::minerva10::minos_mu_reco_eff_highp;
    correction_data10e->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minerva10::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minerva10"] = correction_data10e;

    correction_struct* correction_data13A = new correction_struct;
    correction_data13A->mnv_mu_reco_eff             = resurrection::minerva13A::mnv_mu_reco_eff;
    correction_data13A->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - resurrection::minerva13A::mnv_mu_reco_eff);
    correction_data13A->minos_mu_reco_eff_lowp      = resurrection::minerva13A::minos_mu_reco_eff_lowp;
    correction_data13A->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - resurrection::minerva13A::minos_mu_reco_eff_lowp);
    correction_data13A->minos_mu_reco_eff_highp     = resurrection::minerva13A::minos_mu_reco_eff_highp;
    correction_data13A->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 -  resurrection::minerva13A::minos_mu_reco_eff_highp);
    
    __correction_table["Resurrection"]["minerva13A"] = correction_data13A;
    
    
    correction_struct* correction_data13Ae = new correction_struct;
    correction_data13Ae->mnv_mu_reco_eff             = eroica::minerva13A::mnv_mu_reco_eff;
    correction_data13Ae->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minerva13A::mnv_mu_reco_eff);
    correction_data13Ae->minos_mu_reco_eff_lowp      = eroica::minerva13A::minos_mu_reco_eff_lowp;
    correction_data13Ae->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minerva13A::minos_mu_reco_eff_lowp);
    correction_data13Ae->minos_mu_reco_eff_highp     = eroica::minerva13A::minos_mu_reco_eff_highp;
    correction_data13Ae->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minerva13A::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minerva13A"] = correction_data13Ae;

    correction_struct* correction_data13B = new correction_struct;
    correction_data13B->mnv_mu_reco_eff             = resurrection::minerva13B::mnv_mu_reco_eff;
    correction_data13B->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - resurrection::minerva13B::mnv_mu_reco_eff);
    correction_data13B->minos_mu_reco_eff_lowp      = resurrection::minerva13B::minos_mu_reco_eff_lowp;
    correction_data13B->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - resurrection::minerva13B::minos_mu_reco_eff_lowp);
    correction_data13B->minos_mu_reco_eff_highp     = resurrection::minerva13B::minos_mu_reco_eff_highp;
    correction_data13B->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 -  resurrection::minerva13B::minos_mu_reco_eff_highp);

    __correction_table["Resurrection"]["minerva13B"] = correction_data13B;


    correction_struct* correction_data13Be = new correction_struct;
    correction_data13Be->mnv_mu_reco_eff             = eroica::minerva13B::mnv_mu_reco_eff;
    correction_data13Be->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minerva13B::mnv_mu_reco_eff);
    correction_data13Be->minos_mu_reco_eff_lowp      = eroica::minerva13B::minos_mu_reco_eff_lowp;
    correction_data13Be->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minerva13B::minos_mu_reco_eff_lowp);
    correction_data13Be->minos_mu_reco_eff_highp     = eroica::minerva13B::minos_mu_reco_eff_highp;
    correction_data13Be->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minerva13B::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minerva13B"] = correction_data13Be;   


   correction_struct* correction_data13C = new correction_struct;
    correction_data13C->mnv_mu_reco_eff             = resurrection::minerva13C::mnv_mu_reco_eff;
    correction_data13C->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - resurrection::minerva13C::mnv_mu_reco_eff);
    correction_data13C->minos_mu_reco_eff_lowp      = resurrection::minerva13C::minos_mu_reco_eff_lowp;
    correction_data13C->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - resurrection::minerva13C::minos_mu_reco_eff_lowp);
    correction_data13C->minos_mu_reco_eff_highp     = resurrection::minerva13C::minos_mu_reco_eff_highp;
    correction_data13C->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 -  resurrection::minerva13C::minos_mu_reco_eff_highp);
    
    __correction_table["Resurrection"]["minerva13C"] = correction_data13C;
    
    
    correction_struct* correction_data13Ce = new correction_struct;
    correction_data13Ce->mnv_mu_reco_eff             = eroica::minerva13C::mnv_mu_reco_eff;
    correction_data13Ce->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minerva13C::mnv_mu_reco_eff);
    correction_data13Ce->minos_mu_reco_eff_lowp      = eroica::minerva13C::minos_mu_reco_eff_lowp;
    correction_data13Ce->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minerva13C::minos_mu_reco_eff_lowp);
    correction_data13Ce->minos_mu_reco_eff_highp     = eroica::minerva13C::minos_mu_reco_eff_highp;
    correction_data13Ce->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minerva13C::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minerva13C"] = correction_data13Ce;    


    correction_struct* correction_data2p2he = new correction_struct;
    correction_data2p2he->mnv_mu_reco_eff             = eroica::minerva2p2h::mnv_mu_reco_eff;
    correction_data2p2he->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minerva2p2h::mnv_mu_reco_eff);
    correction_data2p2he->minos_mu_reco_eff_lowp      = eroica::minerva2p2h::minos_mu_reco_eff_lowp;
    correction_data2p2he->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minerva2p2h::minos_mu_reco_eff_lowp);
    correction_data2p2he->minos_mu_reco_eff_highp     = eroica::minerva2p2h::minos_mu_reco_eff_highp;
    correction_data2p2he->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minerva2p2h::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minerva2p2h"] = correction_data2p2he; 

    correction_struct* correction_datanonreweightablese = new correction_struct;
    correction_datanonreweightablese->mnv_mu_reco_eff             = eroica::minervanonreweightables::mnv_mu_reco_eff;
    correction_datanonreweightablese->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minervanonreweightables::mnv_mu_reco_eff);
    correction_datanonreweightablese->minos_mu_reco_eff_lowp      = eroica::minervanonreweightables::minos_mu_reco_eff_lowp;
    correction_datanonreweightablese->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minervanonreweightables::minos_mu_reco_eff_lowp);
    correction_datanonreweightablese->minos_mu_reco_eff_highp     = eroica::minervanonreweightables::minos_mu_reco_eff_highp;
    correction_datanonreweightablese->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minervanonreweightables::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minervanonreweightables"] = correction_datanonreweightablese; 


   correction_struct* correction_data13D = new correction_struct;
    correction_data13D->mnv_mu_reco_eff             = resurrection::minerva13D::mnv_mu_reco_eff;
    correction_data13D->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - resurrection::minerva13D::mnv_mu_reco_eff);
    correction_data13D->minos_mu_reco_eff_lowp      = resurrection::minerva13D::minos_mu_reco_eff_lowp;
    correction_data13D->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - resurrection::minerva13D::minos_mu_reco_eff_lowp);
    correction_data13D->minos_mu_reco_eff_highp     = resurrection::minerva13D::minos_mu_reco_eff_highp;
    correction_data13D->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 -  resurrection::minerva13D::minos_mu_reco_eff_highp);
    
    __correction_table["Resurrection"]["minerva13D"] = correction_data13D;
    
    
    correction_struct* correction_data13De = new correction_struct;
    correction_data13De->mnv_mu_reco_eff             = eroica::minerva13D::mnv_mu_reco_eff;
    correction_data13De->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minerva13D::mnv_mu_reco_eff);
    correction_data13De->minos_mu_reco_eff_lowp      = eroica::minerva13D::minos_mu_reco_eff_lowp;
    correction_data13De->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minerva13D::minos_mu_reco_eff_lowp);
    correction_data13De->minos_mu_reco_eff_highp     = eroica::minerva13D::minos_mu_reco_eff_highp;
    correction_data13De->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minerva13D::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minerva13D"] = correction_data13De;

    correction_struct* correction_data13E = new correction_struct;
    correction_data13E->mnv_mu_reco_eff             = resurrection::minerva13E::mnv_mu_reco_eff;
    correction_data13E->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - resurrection::minerva13E::mnv_mu_reco_eff);
    correction_data13E->minos_mu_reco_eff_lowp      = resurrection::minerva13E::minos_mu_reco_eff_lowp;
    correction_data13E->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - resurrection::minerva13E::minos_mu_reco_eff_lowp);
    correction_data13E->minos_mu_reco_eff_highp     = resurrection::minerva13E::minos_mu_reco_eff_highp;
    correction_data13E->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 -  resurrection::minerva13E::minos_mu_reco_eff_highp);
    
    __correction_table["Resurrection"]["minerva13E"] = correction_data13E;
    
    
    correction_struct* correction_data13Ee = new correction_struct;
    correction_data13Ee->mnv_mu_reco_eff             = eroica::minerva13E::mnv_mu_reco_eff;
    correction_data13Ee->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minerva13E::mnv_mu_reco_eff);
    correction_data13Ee->minos_mu_reco_eff_lowp      = eroica::minerva13E::minos_mu_reco_eff_lowp;
    correction_data13Ee->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minerva13E::minos_mu_reco_eff_lowp);
    correction_data13Ee->minos_mu_reco_eff_highp     = eroica::minerva13E::minos_mu_reco_eff_highp;
    correction_data13Ee->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minerva13E::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minerva13E"] = correction_data13Ee;

   correction_struct* correction_data13r = new correction_struct;
    correction_data13r->mnv_mu_reco_eff             = resurrection::minerva13::mnv_mu_reco_eff;
    correction_data13r->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - resurrection::minerva13::mnv_mu_reco_eff);
    correction_data13r->minos_mu_reco_eff_lowp      = resurrection::minerva13::minos_mu_reco_eff_lowp;
    correction_data13r->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - resurrection::minerva13::minos_mu_reco_eff_lowp);
    correction_data13r->minos_mu_reco_eff_highp     = resurrection::minerva13::minos_mu_reco_eff_highp;
    correction_data13r->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 -  resurrection::minerva13::minos_mu_reco_eff_highp);

    __correction_table["Resurrection"]["minerva13"] = correction_data13r;

   correction_struct* correction_data13e = new correction_struct;
    correction_data13e->mnv_mu_reco_eff             = eroica::minerva13::mnv_mu_reco_eff;
    correction_data13e->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minerva13::mnv_mu_reco_eff);
    correction_data13e->minos_mu_reco_eff_lowp      = eroica::minerva13::minos_mu_reco_eff_lowp;
    correction_data13e->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minerva13::minos_mu_reco_eff_lowp);
    correction_data13e->minos_mu_reco_eff_highp     = eroica::minerva13::minos_mu_reco_eff_highp;
    correction_data13e->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minerva13::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minerva13"] = correction_data13e;

   correction_struct* correction_datame1A = new correction_struct;
    correction_datame1A->mnv_mu_reco_eff             = eroica::minervame1A::mnv_mu_reco_eff;
    correction_datame1A->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minervame1A::mnv_mu_reco_eff);
    correction_datame1A->minos_mu_reco_eff_lowp      = eroica::minervame1A::minos_mu_reco_eff_lowp;
    correction_datame1A->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minervame1A::minos_mu_reco_eff_lowp);
    correction_datame1A->minos_mu_reco_eff_highp     = eroica::minervame1A::minos_mu_reco_eff_highp;
    correction_datame1A->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minervame1A::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minervame1A"] = correction_datame1A;

   correction_struct* correction_datame1B = new correction_struct;
    correction_datame1B->mnv_mu_reco_eff             = eroica::minervame1B::mnv_mu_reco_eff;
    correction_datame1B->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minervame1B::mnv_mu_reco_eff);
    correction_datame1B->minos_mu_reco_eff_lowp      = eroica::minervame1B::minos_mu_reco_eff_lowp;
    correction_datame1B->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minervame1B::minos_mu_reco_eff_lowp);
    correction_datame1B->minos_mu_reco_eff_highp     = eroica::minervame1B::minos_mu_reco_eff_highp;
    correction_datame1B->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minervame1B::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minervame1B"] = correction_datame1B;

   correction_struct* correction_datame1C = new correction_struct;
    correction_datame1C->mnv_mu_reco_eff             = eroica::minervame1C::mnv_mu_reco_eff;
    correction_datame1C->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minervame1C::mnv_mu_reco_eff);
    correction_datame1C->minos_mu_reco_eff_lowp      = eroica::minervame1C::minos_mu_reco_eff_lowp;
    correction_datame1C->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minervame1C::minos_mu_reco_eff_lowp);
    correction_datame1C->minos_mu_reco_eff_highp     = eroica::minervame1C::minos_mu_reco_eff_highp;
    correction_datame1C->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minervame1C::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minervame1C"] = correction_datame1C;

    correction_struct* correction_datame1D = new correction_struct;
    correction_datame1D->mnv_mu_reco_eff             = eroica::minervame1D::mnv_mu_reco_eff;
    correction_datame1D->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minervame1D::mnv_mu_reco_eff);
    correction_datame1D->minos_mu_reco_eff_lowp      = eroica::minervame1D::minos_mu_reco_eff_lowp;
    correction_datame1D->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minervame1D::minos_mu_reco_eff_lowp);
    correction_datame1D->minos_mu_reco_eff_highp     = eroica::minervame1D::minos_mu_reco_eff_highp;
    correction_datame1D->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minervame1D::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minervame1D"] = correction_datame1D;

    correction_struct* correction_datame1E = new correction_struct;
    correction_datame1E->mnv_mu_reco_eff             = eroica::minervame1E::mnv_mu_reco_eff;
    correction_datame1E->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minervame1E::mnv_mu_reco_eff);
    correction_datame1E->minos_mu_reco_eff_lowp      = eroica::minervame1E::minos_mu_reco_eff_lowp;
    correction_datame1E->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minervame1E::minos_mu_reco_eff_lowp);
    correction_datame1E->minos_mu_reco_eff_highp     = eroica::minervame1E::minos_mu_reco_eff_highp;
    correction_datame1E->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minervame1E::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minervame1E"] = correction_datame1E;

    correction_struct* correction_datame1F = new correction_struct;
    correction_datame1F->mnv_mu_reco_eff             = eroica::minervame1F::mnv_mu_reco_eff;
    correction_datame1F->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minervame1F::mnv_mu_reco_eff);
    correction_datame1F->minos_mu_reco_eff_lowp      = eroica::minervame1F::minos_mu_reco_eff_lowp;
    correction_datame1F->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minervame1F::minos_mu_reco_eff_lowp);
    correction_datame1F->minos_mu_reco_eff_highp     = eroica::minervame1F::minos_mu_reco_eff_highp;
    correction_datame1F->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minervame1F::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minervame1F"] = correction_datame1F;

    correction_struct* correction_datame1G = new correction_struct;
    correction_datame1G->mnv_mu_reco_eff             = eroica::minervame1G::mnv_mu_reco_eff;
    correction_datame1G->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minervame1G::mnv_mu_reco_eff);
    correction_datame1G->minos_mu_reco_eff_lowp      = eroica::minervame1G::minos_mu_reco_eff_lowp;
    correction_datame1G->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minervame1G::minos_mu_reco_eff_lowp);
    correction_datame1G->minos_mu_reco_eff_highp     = eroica::minervame1G::minos_mu_reco_eff_highp;
    correction_datame1G->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minervame1G::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minervame1G"] = correction_datame1G;

    correction_struct* correction_datame1L = new correction_struct;
    correction_datame1L->mnv_mu_reco_eff             = eroica::minervame1L::mnv_mu_reco_eff;
    correction_datame1L->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minervame1L::mnv_mu_reco_eff);
    correction_datame1L->minos_mu_reco_eff_lowp      = eroica::minervame1L::minos_mu_reco_eff_lowp;
    correction_datame1L->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minervame1L::minos_mu_reco_eff_lowp);
    correction_datame1L->minos_mu_reco_eff_highp     = eroica::minervame1L::minos_mu_reco_eff_highp;
    correction_datame1L->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minervame1L::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minervame1L"] = correction_datame1L;

    correction_struct* correction_datame1M = new correction_struct;
    correction_datame1M->mnv_mu_reco_eff             = eroica::minervame1M::mnv_mu_reco_eff;
    correction_datame1M->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minervame1M::mnv_mu_reco_eff);
    correction_datame1M->minos_mu_reco_eff_lowp      = eroica::minervame1M::minos_mu_reco_eff_lowp;
    correction_datame1M->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minervame1M::minos_mu_reco_eff_lowp);
    correction_datame1M->minos_mu_reco_eff_highp     = eroica::minervame1M::minos_mu_reco_eff_highp;
    correction_datame1M->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minervame1M::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minervame1M"] = correction_datame1M;

    correction_struct* correction_datame1A2p2h = new correction_struct;
    correction_datame1A2p2h->mnv_mu_reco_eff             = eroica::minervame1A2p2h::mnv_mu_reco_eff;
    correction_datame1A2p2h->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minervame1A2p2h::mnv_mu_reco_eff);
    correction_datame1A2p2h->minos_mu_reco_eff_lowp      = eroica::minervame1A2p2h::minos_mu_reco_eff_lowp;
    correction_datame1A2p2h->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minervame1A2p2h::minos_mu_reco_eff_lowp);
    correction_datame1A2p2h->minos_mu_reco_eff_highp     = eroica::minervame1A2p2h::minos_mu_reco_eff_highp;
    correction_datame1A2p2h->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minervame1A2p2h::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minervame1A2p2h"] = correction_datame1A2p2h;

    correction_struct* correction_datame1B2p2h = new correction_struct;
    correction_datame1B2p2h->mnv_mu_reco_eff             = eroica::minervame1B2p2h::mnv_mu_reco_eff;
    correction_datame1B2p2h->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minervame1B2p2h::mnv_mu_reco_eff);
    correction_datame1B2p2h->minos_mu_reco_eff_lowp      = eroica::minervame1B2p2h::minos_mu_reco_eff_lowp;
    correction_datame1B2p2h->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minervame1B2p2h::minos_mu_reco_eff_lowp);
    correction_datame1B2p2h->minos_mu_reco_eff_highp     = eroica::minervame1B2p2h::minos_mu_reco_eff_highp;
    correction_datame1B2p2h->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minervame1B2p2h::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minervame1B2p2h"] = correction_datame1B2p2h;

    correction_struct* correction_datame1C2p2h = new correction_struct;
    correction_datame1C2p2h->mnv_mu_reco_eff             = eroica::minervame1C2p2h::mnv_mu_reco_eff;
    correction_datame1C2p2h->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minervame1C2p2h::mnv_mu_reco_eff);
    correction_datame1C2p2h->minos_mu_reco_eff_lowp      = eroica::minervame1C2p2h::minos_mu_reco_eff_lowp;
    correction_datame1C2p2h->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minervame1C2p2h::minos_mu_reco_eff_lowp);
    correction_datame1C2p2h->minos_mu_reco_eff_highp     = eroica::minervame1C2p2h::minos_mu_reco_eff_highp;
    correction_datame1C2p2h->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minervame1C2p2h::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minervame1C2p2h"] = correction_datame1C2p2h;

    correction_struct* correction_datame1D2p2h = new correction_struct;
    correction_datame1D2p2h->mnv_mu_reco_eff             = eroica::minervame1D2p2h::mnv_mu_reco_eff;
    correction_datame1D2p2h->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minervame1D2p2h::mnv_mu_reco_eff);
    correction_datame1D2p2h->minos_mu_reco_eff_lowp      = eroica::minervame1D2p2h::minos_mu_reco_eff_lowp;
    correction_datame1D2p2h->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minervame1D2p2h::minos_mu_reco_eff_lowp);
    correction_datame1D2p2h->minos_mu_reco_eff_highp     = eroica::minervame1D2p2h::minos_mu_reco_eff_highp;
    correction_datame1D2p2h->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minervame1D2p2h::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minervame1D2p2h"] = correction_datame1D2p2h;

    correction_struct* correction_datame1E2p2h = new correction_struct;
    correction_datame1E2p2h->mnv_mu_reco_eff             = eroica::minervame1E2p2h::mnv_mu_reco_eff;
    correction_datame1E2p2h->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minervame1E2p2h::mnv_mu_reco_eff);
    correction_datame1E2p2h->minos_mu_reco_eff_lowp      = eroica::minervame1E2p2h::minos_mu_reco_eff_lowp;
    correction_datame1E2p2h->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minervame1E2p2h::minos_mu_reco_eff_lowp);
    correction_datame1E2p2h->minos_mu_reco_eff_highp     = eroica::minervame1E2p2h::minos_mu_reco_eff_highp;
    correction_datame1E2p2h->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minervame1E2p2h::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minervame1E2p2h"] = correction_datame1E2p2h;

    correction_struct* correction_datame1F2p2h = new correction_struct;
    correction_datame1F2p2h->mnv_mu_reco_eff             = eroica::minervame1F2p2h::mnv_mu_reco_eff;
    correction_datame1F2p2h->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minervame1F2p2h::mnv_mu_reco_eff);
    correction_datame1F2p2h->minos_mu_reco_eff_lowp      = eroica::minervame1F2p2h::minos_mu_reco_eff_lowp;
    correction_datame1F2p2h->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minervame1F2p2h::minos_mu_reco_eff_lowp);
    correction_datame1F2p2h->minos_mu_reco_eff_highp     = eroica::minervame1F2p2h::minos_mu_reco_eff_highp;
    correction_datame1F2p2h->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minervame1F2p2h::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minervame1F2p2h"] = correction_datame1F2p2h;

    correction_struct* correction_datame1G2p2h = new correction_struct;
    correction_datame1G2p2h->mnv_mu_reco_eff             = eroica::minervame1G2p2h::mnv_mu_reco_eff;
    correction_datame1G2p2h->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minervame1G2p2h::mnv_mu_reco_eff);
    correction_datame1G2p2h->minos_mu_reco_eff_lowp      = eroica::minervame1G2p2h::minos_mu_reco_eff_lowp;
    correction_datame1G2p2h->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minervame1G2p2h::minos_mu_reco_eff_lowp);
    correction_datame1G2p2h->minos_mu_reco_eff_highp     = eroica::minervame1G2p2h::minos_mu_reco_eff_highp;
    correction_datame1G2p2h->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minervame1G2p2h::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minervame1G2p2h"] = correction_datame1G2p2h;

    correction_struct* correction_datame1L2p2h = new correction_struct;
    correction_datame1L2p2h->mnv_mu_reco_eff             = eroica::minervame1L2p2h::mnv_mu_reco_eff;
    correction_datame1L2p2h->mnv_mu_reco_eff_err         = 0.5 * fabs(1.0 - eroica::minervame1L2p2h::mnv_mu_reco_eff);
    correction_datame1L2p2h->minos_mu_reco_eff_lowp      = eroica::minervame1L2p2h::minos_mu_reco_eff_lowp;
    correction_datame1L2p2h->minos_mu_reco_eff_lowp_err  = 0.5 * fabs(1.0 - eroica::minervame1L2p2h::minos_mu_reco_eff_lowp);
    correction_datame1L2p2h->minos_mu_reco_eff_highp     = eroica::minervame1L2p2h::minos_mu_reco_eff_highp;
    correction_datame1L2p2h->minos_mu_reco_eff_highp_err = 0.5 * fabs(1.0 - eroica::minervame1L2p2h::minos_mu_reco_eff_highp);

    __correction_table["Eroica"]["minervame1L2p2h"] = correction_datame1L2p2h;


    /////////////////////////////////////////////////
    //Inextinguishable

   correction_struct* correction_datame1Ai = new correction_struct;
    correction_datame1Ai->mnv_mu_reco_eff             = inextinguishable::minervame1A::mnv_mu_reco_eff;
    correction_datame1Ai->mnv_mu_reco_eff_err         = 0.6666667 * fabs(1.0 - inextinguishable::minervame1A::mnv_mu_reco_eff);
    correction_datame1Ai->minos_mu_reco_eff_lowp      = inextinguishable::minervame1A::minos_mu_reco_eff_lowp;
    correction_datame1Ai->minos_mu_reco_eff_lowp_err  = 0.6666667 * fabs(1.0 - inextinguishable::minervame1A::minos_mu_reco_eff_lowp);
    correction_datame1Ai->minos_mu_reco_eff_highp     = inextinguishable::minervame1A::minos_mu_reco_eff_highp;
    correction_datame1Ai->minos_mu_reco_eff_highp_err = 0.6666667 * fabs(1.0 - inextinguishable::minervame1A::minos_mu_reco_eff_highp);

    __correction_table["Inextinguishable"]["minervame1A"] = correction_datame1Ai;

   correction_struct* correction_datame1Bi = new correction_struct;
    correction_datame1Bi->mnv_mu_reco_eff             = inextinguishable::minervame1B::mnv_mu_reco_eff;
    correction_datame1Bi->mnv_mu_reco_eff_err         = 0.6666667 * fabs(1.0 - inextinguishable::minervame1B::mnv_mu_reco_eff);
    correction_datame1Bi->minos_mu_reco_eff_lowp      = inextinguishable::minervame1B::minos_mu_reco_eff_lowp;
    correction_datame1Bi->minos_mu_reco_eff_lowp_err  = 0.6666667 * fabs(1.0 - inextinguishable::minervame1B::minos_mu_reco_eff_lowp);
    correction_datame1Bi->minos_mu_reco_eff_highp     = inextinguishable::minervame1B::minos_mu_reco_eff_highp;
    correction_datame1Bi->minos_mu_reco_eff_highp_err = 0.6666667 * fabs(1.0 - inextinguishable::minervame1B::minos_mu_reco_eff_highp);

    __correction_table["Inextinguishable"]["minervame1B"] = correction_datame1Bi;

   correction_struct* correction_datame1Ci = new correction_struct;
    correction_datame1Ci->mnv_mu_reco_eff             = inextinguishable::minervame1C::mnv_mu_reco_eff;
    correction_datame1Ci->mnv_mu_reco_eff_err         = 0.6666667 * fabs(1.0 - inextinguishable::minervame1C::mnv_mu_reco_eff);
    correction_datame1Ci->minos_mu_reco_eff_lowp      = inextinguishable::minervame1C::minos_mu_reco_eff_lowp;
    correction_datame1Ci->minos_mu_reco_eff_lowp_err  = 0.6666667 * fabs(1.0 - inextinguishable::minervame1C::minos_mu_reco_eff_lowp);
    correction_datame1Ci->minos_mu_reco_eff_highp     = inextinguishable::minervame1C::minos_mu_reco_eff_highp;
    correction_datame1Ci->minos_mu_reco_eff_highp_err = 0.6666667 * fabs(1.0 - inextinguishable::minervame1C::minos_mu_reco_eff_highp);

    __correction_table["Inextinguishable"]["minervame1C"] = correction_datame1Ci;

    correction_struct* correction_datame1Di = new correction_struct;
    correction_datame1Di->mnv_mu_reco_eff             = inextinguishable::minervame1D::mnv_mu_reco_eff;
    correction_datame1Di->mnv_mu_reco_eff_err         = 0.6666667 * fabs(1.0 - inextinguishable::minervame1D::mnv_mu_reco_eff);
    correction_datame1Di->minos_mu_reco_eff_lowp      = inextinguishable::minervame1D::minos_mu_reco_eff_lowp;
    correction_datame1Di->minos_mu_reco_eff_lowp_err  = 0.6666667 * fabs(1.0 - inextinguishable::minervame1D::minos_mu_reco_eff_lowp);
    correction_datame1Di->minos_mu_reco_eff_highp     = inextinguishable::minervame1D::minos_mu_reco_eff_highp;
    correction_datame1Di->minos_mu_reco_eff_highp_err = 0.6666667 * fabs(1.0 - inextinguishable::minervame1D::minos_mu_reco_eff_highp);

    __correction_table["Inextinguishable"]["minervame1D"] = correction_datame1Di;

    correction_struct* correction_datame1Ei = new correction_struct;
    correction_datame1Ei->mnv_mu_reco_eff             = inextinguishable::minervame1E::mnv_mu_reco_eff;
    correction_datame1Ei->mnv_mu_reco_eff_err         = 0.6666667 * fabs(1.0 - inextinguishable::minervame1E::mnv_mu_reco_eff);
    correction_datame1Ei->minos_mu_reco_eff_lowp      = inextinguishable::minervame1E::minos_mu_reco_eff_lowp;
    correction_datame1Ei->minos_mu_reco_eff_lowp_err  = 0.6666667 * fabs(1.0 - inextinguishable::minervame1E::minos_mu_reco_eff_lowp);
    correction_datame1Ei->minos_mu_reco_eff_highp     = inextinguishable::minervame1E::minos_mu_reco_eff_highp;
    correction_datame1Ei->minos_mu_reco_eff_highp_err = 0.6666667 * fabs(1.0 - inextinguishable::minervame1E::minos_mu_reco_eff_highp);

    __correction_table["Inextinguishable"]["minervame1E"] = correction_datame1Ei;

    correction_struct* correction_datame1Fi = new correction_struct;
    correction_datame1Fi->mnv_mu_reco_eff             = inextinguishable::minervame1F::mnv_mu_reco_eff;
    correction_datame1Fi->mnv_mu_reco_eff_err         = 0.6666667 * fabs(1.0 - inextinguishable::minervame1F::mnv_mu_reco_eff);
    correction_datame1Fi->minos_mu_reco_eff_lowp      = inextinguishable::minervame1F::minos_mu_reco_eff_lowp;
    correction_datame1Fi->minos_mu_reco_eff_lowp_err  = 0.6666667 * fabs(1.0 - inextinguishable::minervame1F::minos_mu_reco_eff_lowp);
    correction_datame1Fi->minos_mu_reco_eff_highp     = inextinguishable::minervame1F::minos_mu_reco_eff_highp;
    correction_datame1Fi->minos_mu_reco_eff_highp_err = 0.6666667 * fabs(1.0 - inextinguishable::minervame1F::minos_mu_reco_eff_highp);

    __correction_table["Inextinguishable"]["minervame1F"] = correction_datame1Fi;

    correction_struct* correction_datame1Gi = new correction_struct;
    correction_datame1Gi->mnv_mu_reco_eff             = inextinguishable::minervame1G::mnv_mu_reco_eff;
    correction_datame1Gi->mnv_mu_reco_eff_err         = 0.6666667 * fabs(1.0 - inextinguishable::minervame1G::mnv_mu_reco_eff);
    correction_datame1Gi->minos_mu_reco_eff_lowp      = inextinguishable::minervame1G::minos_mu_reco_eff_lowp;
    correction_datame1Gi->minos_mu_reco_eff_lowp_err  = 0.6666667 * fabs(1.0 - inextinguishable::minervame1G::minos_mu_reco_eff_lowp);
    correction_datame1Gi->minos_mu_reco_eff_highp     = inextinguishable::minervame1G::minos_mu_reco_eff_highp;
    correction_datame1Gi->minos_mu_reco_eff_highp_err = 0.6666667 * fabs(1.0 - inextinguishable::minervame1G::minos_mu_reco_eff_highp);

    __correction_table["Inextinguishable"]["minervame1G"] = correction_datame1Gi;

    correction_struct* correction_datame1Li = new correction_struct;
    correction_datame1Li->mnv_mu_reco_eff             = inextinguishable::minervame1L::mnv_mu_reco_eff;
    correction_datame1Li->mnv_mu_reco_eff_err         = 0.6666667 * fabs(1.0 - inextinguishable::minervame1L::mnv_mu_reco_eff);
    correction_datame1Li->minos_mu_reco_eff_lowp      = inextinguishable::minervame1L::minos_mu_reco_eff_lowp;
    correction_datame1Li->minos_mu_reco_eff_lowp_err  = 0.6666667 * fabs(1.0 - inextinguishable::minervame1L::minos_mu_reco_eff_lowp);
    correction_datame1Li->minos_mu_reco_eff_highp     = inextinguishable::minervame1L::minos_mu_reco_eff_highp;
    correction_datame1Li->minos_mu_reco_eff_highp_err = 0.6666667 * fabs(1.0 - inextinguishable::minervame1L::minos_mu_reco_eff_highp);

    __correction_table["Inextinguishable"]["minervame1L"] = correction_datame1Li;

    correction_struct* correction_datame1Mi = new correction_struct;
    correction_datame1Mi->mnv_mu_reco_eff             = inextinguishable::minervame1M::mnv_mu_reco_eff;
    correction_datame1Mi->mnv_mu_reco_eff_err         = 0.6666667 * fabs(1.0 - inextinguishable::minervame1M::mnv_mu_reco_eff);
    correction_datame1Mi->minos_mu_reco_eff_lowp      = inextinguishable::minervame1M::minos_mu_reco_eff_lowp;
    correction_datame1Mi->minos_mu_reco_eff_lowp_err  = 0.6666667 * fabs(1.0 - inextinguishable::minervame1M::minos_mu_reco_eff_lowp);
    correction_datame1Mi->minos_mu_reco_eff_highp     = inextinguishable::minervame1M::minos_mu_reco_eff_highp;
    correction_datame1Mi->minos_mu_reco_eff_highp_err = 0.6666667 * fabs(1.0 - inextinguishable::minervame1M::minos_mu_reco_eff_highp);

    __correction_table["Inextinguishable"]["minervame1M"] = correction_datame1Mi;

    correction_struct* correction_datame1Ni = new correction_struct;
    correction_datame1Ni->mnv_mu_reco_eff             = inextinguishable::minervame1N::mnv_mu_reco_eff;
    correction_datame1Ni->mnv_mu_reco_eff_err         = 0.6666667 * fabs(1.0 - inextinguishable::minervame1N::mnv_mu_reco_eff);
    correction_datame1Ni->minos_mu_reco_eff_lowp      = inextinguishable::minervame1N::minos_mu_reco_eff_lowp;
    correction_datame1Ni->minos_mu_reco_eff_lowp_err  = 0.6666667 * fabs(1.0 - inextinguishable::minervame1N::minos_mu_reco_eff_lowp);
    correction_datame1Ni->minos_mu_reco_eff_highp     = inextinguishable::minervame1N::minos_mu_reco_eff_highp;
    correction_datame1Ni->minos_mu_reco_eff_highp_err = 0.6666667 * fabs(1.0 - inextinguishable::minervame1N::minos_mu_reco_eff_highp);

    __correction_table["Inextinguishable"]["minervame1N"] = correction_datame1Ni;

    correction_struct* correction_datame1Oi = new correction_struct;
    correction_datame1Oi->mnv_mu_reco_eff             = inextinguishable::minervame1O::mnv_mu_reco_eff;
    correction_datame1Oi->mnv_mu_reco_eff_err         = 0.6666667 * fabs(1.0 - inextinguishable::minervame1O::mnv_mu_reco_eff);
    correction_datame1Oi->minos_mu_reco_eff_lowp      = inextinguishable::minervame1O::minos_mu_reco_eff_lowp;
    correction_datame1Oi->minos_mu_reco_eff_lowp_err  = 0.6666667 * fabs(1.0 - inextinguishable::minervame1O::minos_mu_reco_eff_lowp);
    correction_datame1Oi->minos_mu_reco_eff_highp     = inextinguishable::minervame1O::minos_mu_reco_eff_highp;
    correction_datame1Oi->minos_mu_reco_eff_highp_err = 0.6666667 * fabs(1.0 - inextinguishable::minervame1O::minos_mu_reco_eff_highp);

    __correction_table["Inextinguishable"]["minervame1O"] = correction_datame1Oi;

    correction_struct* correction_datame1Pi = new correction_struct;
    correction_datame1Pi->mnv_mu_reco_eff             = inextinguishable::minervame1P::mnv_mu_reco_eff;
    correction_datame1Pi->mnv_mu_reco_eff_err         = 0.6666667 * fabs(1.0 - inextinguishable::minervame1P::mnv_mu_reco_eff);
    correction_datame1Pi->minos_mu_reco_eff_lowp      = inextinguishable::minervame1P::minos_mu_reco_eff_lowp;
    correction_datame1Pi->minos_mu_reco_eff_lowp_err  = 0.6666667 * fabs(1.0 - inextinguishable::minervame1P::minos_mu_reco_eff_lowp);
    correction_datame1Pi->minos_mu_reco_eff_highp     = inextinguishable::minervame1P::minos_mu_reco_eff_highp;
    correction_datame1Pi->minos_mu_reco_eff_highp_err = 0.6666667 * fabs(1.0 - inextinguishable::minervame1P::minos_mu_reco_eff_highp);

    __correction_table["Inextinguishable"]["minervame1P"] = correction_datame1Pi;

    correction_struct* correction_datame5Ai = new correction_struct;
    correction_datame5Ai->mnv_mu_reco_eff             = inextinguishable::minervame5A::mnv_mu_reco_eff;
    correction_datame5Ai->mnv_mu_reco_eff_err         = 0.6666667 * fabs(1.0 - inextinguishable::minervame5A::mnv_mu_reco_eff);
    correction_datame5Ai->minos_mu_reco_eff_lowp      = inextinguishable::minervame5A::minos_mu_reco_eff_lowp;
    correction_datame5Ai->minos_mu_reco_eff_lowp_err  = 0.6666667 * fabs(1.0 - inextinguishable::minervame5A::minos_mu_reco_eff_lowp);
    correction_datame5Ai->minos_mu_reco_eff_highp     = inextinguishable::minervame5A::minos_mu_reco_eff_highp;
    correction_datame5Ai->minos_mu_reco_eff_highp_err = 0.6666667 * fabs(1.0 - inextinguishable::minervame5A::minos_mu_reco_eff_highp);
  
    __correction_table["Inextinguishable"]["minervame5A"] = correction_datame5Ai;

    correction_struct* correction_datame6Ai = new correction_struct;
    correction_datame6Ai->mnv_mu_reco_eff             = inextinguishable::minervame6A::mnv_mu_reco_eff;
    correction_datame6Ai->mnv_mu_reco_eff_err         = 0.6666667 * fabs(1.0 - inextinguishable::minervame6A::mnv_mu_reco_eff);
    correction_datame6Ai->minos_mu_reco_eff_lowp      = inextinguishable::minervame6A::minos_mu_reco_eff_lowp;
    correction_datame6Ai->minos_mu_reco_eff_lowp_err  = 0.6666667 * fabs(1.0 - inextinguishable::minervame6A::minos_mu_reco_eff_lowp);
    correction_datame6Ai->minos_mu_reco_eff_highp     = inextinguishable::minervame6A::minos_mu_reco_eff_highp;
    correction_datame6Ai->minos_mu_reco_eff_highp_err = 0.6666667 * fabs(1.0 - inextinguishable::minervame6A::minos_mu_reco_eff_highp);
  
    __correction_table["Inextinguishable"]["minervame6A"] = correction_datame6Ai;

    correction_struct* correction_datame6Bi = new correction_struct;
    correction_datame6Bi->mnv_mu_reco_eff             = inextinguishable::minervame6B::mnv_mu_reco_eff;
    correction_datame6Bi->mnv_mu_reco_eff_err         = 0.6666667 * fabs(1.0 - inextinguishable::minervame6B::mnv_mu_reco_eff);
    correction_datame6Bi->minos_mu_reco_eff_lowp      = inextinguishable::minervame6B::minos_mu_reco_eff_lowp;
    correction_datame6Bi->minos_mu_reco_eff_lowp_err  = 0.6666667 * fabs(1.0 - inextinguishable::minervame6B::minos_mu_reco_eff_lowp);
    correction_datame6Bi->minos_mu_reco_eff_highp     = inextinguishable::minervame6B::minos_mu_reco_eff_highp;
    correction_datame6Bi->minos_mu_reco_eff_highp_err = 0.6666667 * fabs(1.0 - inextinguishable::minervame6B::minos_mu_reco_eff_highp);
  
    __correction_table["Inextinguishable"]["minervame6B"] = correction_datame6Bi;
}

//================================================
// LoadPlaylist
//================================================
void MnvNormalizer::LoadPlaylist( const std::string& playlist_name )
{

  std::string playlist_lookup = playlist_name;
  // Convert "playlistX" to "minervaX"
  const std::string playlistword = "playlist";
  if ( playlist_lookup.compare( 0, playlistword.length(), playlistword ) == 0 )
  {
    playlist_lookup.replace( 0, playlistword.length(), "minerva" );
  }
  // Retrieve corrections for specified processing and playlist
  __the_correction_struct = __correction_table[m_processing][playlist_lookup];
  if ( !__the_correction_struct )
  {
    std::cerr << Form( "MnvNormalizer: Corrections not found for processing \"%s\", playlist \"%s\"", m_processing.c_str(), playlist_name.c_str() ) << std::endl;
    exit(1);
  }

  //std::cout << Form( "MnvNormalizer: Loaded corrections for processing \"%s\", playlist \"%s\"", m_processing.c_str(), playlist_lookup.c_str() ) << std::endl;

  // Do these member variables serve a purpose?
  //set current config member variables
  playlist_ = playlist_name;
  analysis_ = "";

  // Mass model scale + error
  mass_scale = MnvNorm::mass_scale;
  mass_scale_err = MnvNorm::mass_scale_err;

  //Dead time / pile-up in MINERvA + error due to tdead cut.
  dead_time_tdead = MnvNorm::dead_time_tdead;
  dead_time_tdead_err = MnvNorm::dead_time_tdead_err;

}

//================================================
// LoadAnalysis
//================================================
void MnvNormalizer::LoadAnalysis( const std::string& analysis )
{
  //set current config member variable for analysis only
  analysis_ = analysis;

  /*
     don't use this.  just an example
     dead_time_tdead = MnvNorm::ccqe::dead_time_tdead;
     dead_time_tdead_err = MnvNorm::ccqe::dead_time_tdead_err;
   */

  std::cout << "MnvNormalizer::LoadAnalysis ERROR: Cannot load unknown analysis: " << analysis << std::endl;
  throw 23;
}

//================================================
// Load
//================================================
void MnvNormalizer::Load( const std::string& playlist, const std::string& analysis )
{
  //load the playlist defaults
  LoadPlaylist( playlist );

  //load the analysis defaults
  LoadAnalysis( analysis );

  //set current config member variables
  playlist_ = playlist;
  analysis_ = analysis;

  //now set anything is specific to both playlist and analysis

  /* not actually used.  just an example.
     if( "minerva5" == playlist_ && "ccqe" == analysis_ )
     {
     dead_time_tdead = MnvNorm::minerva5::ccqe::dead_time_tdead;
     dead_time_tdead_err = MnvNorm::minerva5::ccqe::dead_time_tdead_err;
     }
   */
}



/*! Get the total correction factor
  @param[in] minosP momentum of muon at face of MINOS in MeV/c
  @return total correction factor
 */
double MnvNormalizer::Correction( double minosP ) const
{
  //multiply all corrections.  the 1 is just to make lines easier to comment out if necessary
  double corr = 1.
//  * mass_scale // mass model should be separate from reconstruction efficiency
    * __the_correction_struct->mnv_mu_reco_eff
    ;

  if( minosP < MnvNorm::minosP_minos_trk_eff_threshP )
  {
    corr *= __the_correction_struct->minos_mu_reco_eff_lowp;
  }else{
    corr *= __the_correction_struct->minos_mu_reco_eff_highp;
  }

  return corr;
}

double MnvNormalizer::GetCorrection(double minosP) const
{
  return Correction(minosP);
}


/*! Get the error on the total correction factor
  @param[in] minosP momentum of muon at face of MINOS in MeV/c
  @return absolute error on the correction factor
 */
double MnvNormalizer::CorrectionErr( double minosP ) const
{
  //sum all errors in quadrature.  the 0 is just to make lines easier to comment out if necessary
  double err = sqrt( 0.
//    + mass_scale_err*mass_scale_err // mass model should be separate from reconstruction efficiency
      + __the_correction_struct->mnv_mu_reco_eff_err*__the_correction_struct->mnv_mu_reco_eff_err
      + dead_time_tdead_err*dead_time_tdead_err
      );

  if( minosP < MnvNorm::minosP_minos_trk_eff_threshP )
  {
    err = sqrt( err*err + __the_correction_struct->minos_mu_reco_eff_lowp_err*__the_correction_struct->minos_mu_reco_eff_lowp_err );
  }else{
    err = sqrt( err*err + __the_correction_struct->minos_mu_reco_eff_highp_err*__the_correction_struct->minos_mu_reco_eff_highp_err );
  }

  return err;
}

double MnvNormalizer::GetCorrectionErr(double minosP) const
{
  return CorrectionErr(minosP);
}

/*! Get the mass model scale
@return mass model scale
*/
double MnvNormalizer::MassModelScale() const
{
    return MnvNorm::mass_scale;
}

double MnvNormalizer::GetMassModelScale() const
{
    return MassModelScale();
}

/*! Get the error on the mass model scale
@return error on the mass model scale
*/
double MnvNormalizer::MassModelErr() const
{
    return MnvNorm::mass_scale_err;
}

double MnvNormalizer::GetMassModelErr() const
{
    return MassModelErr();
}
