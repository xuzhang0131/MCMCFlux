/*
 * Constants.h
 *
 *  Created on: Oct 15, 2018
 *      Author: zhangxu
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <math.h>
#include <time.h>
#include <Eigen/Dense>
#include <vector>

using namespace std;

extern const int NumberOfIsotopomers (512);

extern const double CellVolume (0.5E-6);      //  u Liter
extern const double MediumVolume (7.8e3);       // u Liter        7.8 mL
extern const int TissueVolume (26);    // 0.026 mL = 26 u Liter
extern const double CytoplasmVolume = TissueVolume * 0.64;     // 0.64 * 26 u Liter
extern const double MitochonVolume = TissueVolume * 0.16;     // 0.16 * 26 u Liter

//extern const double MoleculeNumberInOneNanoMole = 6.02214129e14;

extern const int CO2 = 1;
extern const int NH3 = 1;
extern const int cCO2 = 1;
extern const int cNH3 = 1;
extern const int mCO2 = 1;
extern const int mNH3 = 1;
extern const double AvogadroConstant = 6.02214129E+23;

const double Kf_GLC_Transport = 0.00073 * 0.3;                 // / min
const double Kr_GLC_Transport = 0.288 * 0.3;                 // / min
const double Kf_GLC_SerPool1 = 0.0011 * 0.4;                    // / min
const double Kf_GLC_PYR = 0.1833;                    // / min
const double Kf_GLC_PRPP = 0.000083 * 1.5;                   // / min
const double Kf_SerPool3_Synthesis = 6.4e-6 * 1.2;            // M /min
const double Kf_SerPool3_Degradation = 1.77e-2;                     // /min
const double Kf_GlyPool3_Synthesis = 6.4e-6 * 3;                 // M /min
const double Kf_GlyPool3_Degradation = 0.53e-3;                // /min
const double Kf_SerPool1_GlyPool1 = 176 * 1.5;         // /M /min
const double Kr_SerPool1_GlyPool1 = 3.55 * 0.75;         // /M /min
const double Kf_Ser_Transport = 0.000139 * 1.5;            // /min
const double Kr_Ser_Transport = 0.00667 * 1.5;                // /min
const double Kf_Gly_Transport = 0.0000603 * 1.5;            // /min
const double Kr_Gly_Transport = 0.000694 * 1.5;            // /min
const double Kf_SerPool2_SerPoolMitochon = 0.0013 * 3;                // /min
const double Kr_SerPool2_SerPoolMitochon = 0.00325 * 3;            // /min
const double Kf_GlyPool2_GlyPoolMitochon = 0.00195 * 3;            // /min
const double Kr_GlyPool2_GlyPoolMitochon = 0.004875 * 3;            // /min
const double Kf_SerPoolMitochon_GlyPoolMitochon = 277 * 3;        // /M /min
const double Kr_SerPoolMitochon_GlyPoolMitochon = 5.5 * 3;        // /M /min
const double Kf_GlyPool1_CO2 = 17.7 * 3;                        // /M /min
const double Kr_GlyPool1_CO2 = 0.213 * 3;                        // /min
const double Kf_MethyleneTHF_FormylTHF_Cytoplasm = 0.32;            // /min
const double Kr_MethyleneTHF_FormylTHF_Cytoplasm = 0.8;            // /min
const double Kf_FormylTHF_Formate_Cytoplasm = 0.32;                        // /min
const double Kr_FormylTHF_Formate_Cytoplasm = 1.0e6;             // /M /min
const double Kf_GlyPoolMitochon_CO2 = 434;                        // /M /min
const double Kr_GlyPoolMitochon_CO2 = 0.42;                     // /min
const double Kf_MethyleneTHF_FormylTHF_Mitochon = 2.67;                    // /min
const double Kr_MethyleneTHF_FormylTHF_Mitochon = 0.267;                    // /min
const double Kf_FormylTHF_Formate_Mitochon = 2.67;                            // /min
const double Kr_FormylTHF_Formate_Mitochon = 4.44e4;            // /M /min
const double Kf_PRPP1_GAR = 44 * 2;             // /M /min
const double Kf_GAR_FGAR = 88e3 * 0.5;             // /M /min
const double Kf_FGAR_AMP = 88e3 * 0.5;             // /M /min
const double Kf_IMP_AMP = 0.025 * 0.5;        // /min
const double Kf_AMP_Degradation = 0.00125 * 0.5;         // /min
const double Kr_MethyleneTHF_Transport = 0.013 * 3;             // /min
const double Kf_MethyleneTHF_Transport = 0.026 * 0.3;         // /min
const double Kr_FormylTHF_Transport = 2.34 * 2;             // /min
const double Kf_FormylTHF_Transport = 0.468;                 // /min
const double Kr_THF_Transport = 0.78;                         // /min
const double Kf_THF_Transport = 0.156 * 0.5;         // /min
const double Kr_Formate_Transport = 4.68;                         // /min
const double Kf_Formate_Transport = 0.936 * 0.5;         // /min
const double Kf_PRPP2_GAR = 22 * 2;              // /M /min
const double Kf_PRPP3_GAR = 11 * 2;              // /M /min
const double Kf_SerPool2_GlyPool2 = 176 * 1.5;            // /M /min
const double Kr_SerPool2_GlyPool2 = 3.55 * 0.75;        // /M /min
const double Kf_GlyPool2_CO2 = 17.7 * 3;                // /M /min
const double Kr_GlyPool2_CO2 = 0.213 * 3;             // /min
const double& Kf_SerPool3_GlyPool3 = Kf_SerPool2_GlyPool2;
const double& Kr_SerPool3_GlyPool3 = Kr_SerPool2_GlyPool2;
const double& Kf_GlyPool3_CO2 = Kf_GlyPool2_CO2;
const double& Kr_GlyPool3_CO2 = Kr_GlyPool2_CO2;
const double& Kf_Ser_Transport3 = Kf_Ser_Transport;
const double& Kr_Ser_Transport3 = Kr_Ser_Transport;
const double& Kf_Gly_Transport3 = Kf_Gly_Transport;
const double& Kr_Gly_Transport3 = Kr_Gly_Transport;

//typedef boost::array< double , NumberOfIsotopomers > state_type;

extern const double OptimizeParameters[57] = {0.000461657750477753  ,  0.427503086639990  ,  0.000138933974778344 ,   0.912957765988904 ,   0.000302264234311588  ,
    0.000264216988034605 ,   0.528433976069209  ,  0.000178629396143585  ,  0.0105686795213842 ,   389.997313257504 ,
    3.65988124056536 ,   0.000422747180855367 ,   0.0496192767065513   , 7.44289150598270e-05 ,   0.0178629396143585  ,
    0.0109460124414652 ,   0.0140352064043982 ,   0.00966781987350446  ,  0.0490562874450916  ,  830.998723620212,
    17.0631986313531 ,   54.1720082557871,    1.25536770067575 ,   1.67337425755250 ,   2.03451173487994,    0.399124270787164 ,
    1000000.00000000 ,   434.000150742915 ,   1.64330339153493  ,  2.22492836752176 ,   0.600557183706194 ,   2.29538774044506 ,
    44400.0000000432  ,  89.9889179513641  ,  44000.0000016264 ,   43999.9999979280  ,  0.0266554754467594 ,   0.00151731008995339,
    0.189003218774087  ,  0.00751136610783774  ,  11.6715973413360  ,  0.495498097191622 ,   2.51038222698186  ,
    0.0677700081258078  ,  4.57755449879482,    0.579744450285708 ,   42.5463885164964 ,   21.0649543138687 ,   85.0037419701686  ,
    120.001190410154 ,   65.0980133729533  ,  1.78629396143585 ,   CellVolume  ,  MediumVolume  ,  TissueVolume ,   CytoplasmVolume ,   MitochonVolume};

extern const double OptimizeParameters2[57] = {0.000465165222015355, 0.242665045280639, 0.000036, 0.518326766571717, 0.0002, 0.00009, 0.2,
    0.000075, 0.075, 130.0,   0.6, 0.00018 , 0.075, 0.00005, 0.12 , 0.001,0.0008 , 0.00974237918290052,   0.0278492039049790 , 780.0 ,
    0.6, 54.2492566486518, 0.03 , 0.95 ,   2.00461377104628, 0.398825587055585, 999999.999999978,  434.005733100667,  1.05698640506409,
    2.24239831194499 , 0.551424353965021, 2.31291778438686, 44400.0000092611,  80.0, 44000.0001051987,  44000.0000615634, 0.0268575860192817,
    0.000861408841917836,   0.107341967375275, 0.00756915225856982, 11.6618608584811,  0.499324532242705, 2.01627621832653 ,  0.0682897697038028,
    4.57835999241956, 0.579450741578679, 42.0, 21.0, 120.0, 1.2,  65.0,  0.72, CellVolume,  MediumVolume,  TissueVolume,  CytoplasmVolume,
    MitochonVolume};

#endif /* CONSTANTS_H_ */

