#ifndef vicfuns_H__
#define vicfuns_H__
typedef double m_elem;
#include <brtenv.h>		/* Basic real-time environment */
#include <io1103.h>		/* io1103 macros*/
#include <int1103.h>	/* Interrupt handling */
#include <math.h>		/* Math macros */
#include <dsstd.h>		/* Standard macros */
#include <slvdsp1103.h> /* Slave macro */
#include <stdio.h>
#include <stdlib.h>
#include <io1103.h>		/* io1103 macros*/

//typedef double Float64;
//typedef int Int16;
#define pi 3.14159265359
#define _2pi3 2.0943951023932
#define _4pi3 4.18879020478639
 Float64 kgain1=2.653, kgain2=2.656, kgain3=2.651, kgain4=2.644, kgain5=2.647, kgain6=2.644, kgain7=1.0006, kgain8=1.0005, kgain9=1.0005; //variac

Float64 fv=2.0;
Float64 fi=200.0;
Float64 fio=200.0;
Float64 fis=130.0;
 
////////Potencias
Float64 Po[3]={0.0,0.0,0.0};
Float64 Ps=0.0;
Float64 Pm=0.0;
Float64 Pz[3]={0.0,0.0,0.0};
Float64 Pxy[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Podq[2] = {0.0,0.0};

Float64 Pod[5]={0.0,0.0,0.0,0.0,0.0};
Float64 Poq[5]={0.0,0.0,0.0,0.0,0.0};

Float64 PoqF[5]={0.0,0.0,0.0,0.0,0.0};
Float64 PodF[5]={0.0,0.0,0.0,0.0,0.0};
Float64 PodqF[2] = {0.0,0.0};
Float64 pod_ref = 0.0;
Float64 poq_ref = 0.0;
Float64 err_iz[2] = {0.0,0.0};
Float64 err_iz_1[2] = {0.0,0.0};
Float64 Out_iz[2] = {0.0,0.0};
Float64 Out_iz_1[2] = {0.0,0.0};
Float64 max_iz = 3.0*197.0;
Float64 min_iz = -3.0*197.0;
Float64 Izdq_iny[2] = {0.0,0.0};
Float64 kp_iz = -0.4;
Float64 ki_iz = -5.0;

Float64 Pzd[5]={0.0,0.0,0.0,0.0,0.0};
Float64 Pzq[5]={0.0,0.0,0.0,0.0,0.0};
Float64 Pzdq[2] = {0.0,0.0};
Float64 PzqF[5]={0.0,0.0,0.0,0.0,0.0};
Float64 PzdF[5]={0.0,0.0,0.0,0.0,0.0};
Float64 PzdqF[2] = {0.0,0.0};
Float64 pzd_ref = 0.0;
Float64 pzq_ref = 0.0;
Float64 err_vm[2] = {0.0,0.0};
Float64 err_vm_1[2] = {0.0,0.0};
Float64 Out_vm[2] = {0.0,0.0};
Float64 Out_vm_1[2] = {0.0,0.0};
Float64 max_vm = 60.0;
Float64 min_vm = -60.0;
Float64 Vmdq_iny[2] = {0.0,0.0};
Float64 kp_vm = -0.6;
Float64 ki_vm = -4.0;

//[Bpo,Apo]=ellip(4,.1,20,(8*Ts)*2); 
UInt32 inyect_IZ=0;
UInt32 inyect_VM=0;

// Float64 Bpo[5]={0.09989361743041,-0.399573038300251,0.599358841742866, -0.399573038300251, 0.099893617430412};
// Float64 Apo[5]={1.0,-3.997863102043092,5.993593487062134,-3.993597663977122, 0.997867278961305};
//Float64 Bpo[4] = {0.0000755358809314749,-0.0000755358516796430,-0.0000755358516796429, 0.0000755358809314749};
//Float64 Apo[4]={1.0,-2.999411707976560, 2.998823655069235,-0.999411947034171};
Float64 Bpo[4] = {0.000753368607151969,-0.000753339432624164,-0.000753339432624161,  0.000753368607151969};
Float64 Apo[4]={1.0,-2.994111137176799,2.988246201525302,-0.994135005999447};


Float64 Bpo2[7]={0.001502352101877,-0.006008464009620,0.007509871797898,0.0, -0.007509871797898, 0.006008464009620,-0.001502352101877};
Float64 Apo2[7]={1.0,-5.987499887989729,14.938310150730249,-19.878235336037520, 14.879844175158484,-5.940723504480928,0.988304402632659};
Float64 Iod[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Ioq[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
Float64 IodF[4] = {0.0,0.0,0.0,0.0};
Float64 IoqF[4] = {0.0,0.0,0.0,0.0};
Float64 IodqT[2] = {0.0,0.0};
Float64 Bio[4] = {0.039212174010877,-0.033899238060784 ,-0.033899238060784, 0.039212174010877};
Float64 Aio[4]={1.0,-2.624958536285712,2.337998370104467,-0.702413961918569};
Float64 BioF[6] = {0.028761583201554,0.143095989407480 ,0.328308857607255, 0.328308857607255, 0.143010335270111,0.028680807122634};

// Float64 Bpo2[5]={0.099329596170799,-0.397141600876095,0.595624014981739,-0.397141600876095, 0.099329596170799};
// Float64 Apo2[5]={1.0,-5.987499887989729,14.938310150730249,-19.878235336037520, 14.879844175158484,-5.940723504480928,0.988304402632659};

Float64 Bpz[5]={0.099853913144369,-0.399412947377363,0.599118068477383,-0.399412947377363, 0.099853913144369};
Float64 Apz[5]={1.0,-3.997060786669827,5.991190265693001,-3.991198160929517, 0.997068681917871};
int Po_order = 3;
int Po_order_BP = 6;
int Pz_order = 3;
Float64 Poa[7] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Pob[7] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Poc[7] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
Float64 PoaF[7] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
Float64 PobF[7] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
Float64 PocF[7] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Poabc[3] = {0.0,0.0,0.0};

Float64 ang=0.0;
Float64 ang_1=0.0;
Float64 Wred=0.0;
Float64 Umax_pll = 9900000.0;
Float64 Umin_pll = -9900000.0;
Float64 Vq_1 = 0.0;
Float64 kp_pll=2.0;
Float64 ki_pll=0.01;
Float64 Wred_1 = 0.0;
Float64 Vdq[2]={0.0,0.0};
Float64 Vq = 0.0;
Float64 Ts=25e-6;
Float64 Ti = 4.0;
Float64 Tv = 9.0;
Float64 IodFir;
Float64 IoqFir;

UInt32 L1,L2,L3,U3,U2,U1;
UInt32 signpa=0;
UInt32 signpb=0;
UInt32 signpc=0;
UInt32 protect=0;
UInt32 protect2 = 0;
UInt32 signna=0;
UInt32 signnb=0;
UInt32 signnc=0;
UInt32 update_control=0;
UInt32 levelAp=4;
UInt32 levelBp=4;
UInt32 levelCp=4;

UInt32 levelAn=4;
UInt32 levelBn=4;
UInt32 levelCn=4;
Float64 B1,B2,B3,B4,B5,B6,B7,A2,A3,A4,A5,A6,A7;
Float64 B1i,B2i,B3i,B4i,B5i,B6i,B7i,A2i,A3i,A4i,A5i,A6i,A7i;
Float64 B1z,B2z,B3z,B4z,B5z,B6z,B7z,A2z,A3z,A4z,A5z,A6z,A7z;

Float64 F[9]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Fii[9]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

Float64 Vcfir[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Vcfir_1[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Vcfir_2[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Vcfir_3[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Vcfir_4[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Vcfir_5[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Vcfir_6[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Vcfir_7[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Vcfir_8[6]={0.0,0.0,0.0,0.0,0.0,0.0};


Float64 VcxyIn[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 VcxyIn_1[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 VcxyIn_2[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 VcxyIn_3[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 VcxyIn_4[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 VcxyIn_5[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 VcxyIn_6[6]={0.0,0.0,0.0,0.0,0.0,0.0};

Float64 VcxyO[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 VcxyO_1[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 VcxyO_2[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 VcxyO_3[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 VcxyO_4[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 VcxyO_5[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 VcxyO_6[6]={0.0,0.0,0.0,0.0,0.0,0.0};

Float64 IxyF[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 IxyO[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 IxyF_1[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 IxyF_2[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 IxyF_3[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 IxyF_4[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 IxyF_5[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 IxyF_6[6]={0.0,0.0,0.0,0.0,0.0,0.0};



Float64 IxyO_1[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 IxyO_2[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 IxyO_3[6]={0.0,0.0,0.0,0.0,0.0,0.0};

Float64 IZz[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 IZz_1[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 IZz_2[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 IZz_3[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 IZz_4[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 IZz_5[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 IZz_6[6]={0.0,0.0,0.0,0.0,0.0,0.0};

Float64 IzO[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 IzO_1[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 IzO_2[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 IzO_3[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 IzO_4[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 IzO_5[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 IzO_6[6]={0.0,0.0,0.0,0.0,0.0,0.0};

UInt32 Control_maximus=0;
Float64 vcpa,vcpb,vcpc,vcna,vcnb,vcnc;
////////////////////////
UInt32 fpga;
UInt32 readi;
int read_state = 0;
UInt32 error_counter=0;
UInt32 warn_cero=0;
///////////////////////



//Parametros del filtro con frecuencia de corte en 60 Hz con muestreo de 192e-6, modificar de acuerdo al muestreo
double b[4]={0.0052,-0.0051,-0.0051,0.0052};
double a[4]={1.0,-2.9258,2.8582,-0.9322};

Float64 Ixy[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Vxy[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Exy[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Mxy[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 modxy[6]={0.0,0.0,0.0,0.0,0.0,0.0};

Float64 vcap;
Float64 vcbp;
Float64 vccp;
Float64 vcan;
Float64 vcbn;
Float64 vccn;
Float64 is_sign=1.0;
Float64 Im=0;
Float64 Io[3]={0.0,0.0,0.0};
Float64 Iz[3]={0.0,0.0,0.0};
Float64 Is=0.0;
Float64 Is_est = 0.3;
Float64 Vo[3]={0.0,0.0,0.0};
Float64 Vz[3]={0.0,0.0,0.0};
Float64 Vs=0.0;
Float64 Vm=0.0;
Float64 Eo[3]={0.0,0.0,0.0};
Float64 Ez[3]={0.0,0.0,0.0};
Float64 Es=0.0;
Float64 Em=0.0;

Float64 Iodq_ref[2]={0.0,0.0};
Float64 Izdq_ref[2]={0.0,0.0};

Float64 Io_ref[3]={0.0,0.0,0.0};
Float64 Iz_refT[3]={0.0,0.0,0.0};
Float64 Io_refT[3]={0.0,0.0,0.0};

Float64 Iodq[2]={0.0,0.0};
Float64 Izdq[2]={0.0,0.0};
Float64 Vmdq[2]={0.0,0.0};
Float64 d_ioref[3]={0.0,0.0,0.0};
Float64 d_izref[3]={0.0,0.0,0.0};
Float64 d_isref=0.0;
Float64 d_vmref=0.0;
Float64 Iz_ref[3]={0.0,0.0,0.0};
Float64 Vm_ripple=0.0;
Float64 Vmref=0.0;
Float64 Exy_ref ;
Float64 fo=50.0; //frecuencia de corriente de salida
Float64 wn_v ;
Float64 wn_i,wn_io,wn_is ;
Float64 xi=0.707;
Float64 h;
Float64 hv;
Float64 m1;
Float64 m2;
Float64 Vo_mod[2]={30.0,0.0};
Float64 vcpa_off,vcpb_off,vcpc_off,vcna_off,vcnb_off,vcnc_off;
 int counter = 0;
 UInt32 enable_mod = 0;
UInt32 voltage[6]={0,0,0,0,0,0};
Float64 voltage_t[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 kpv ;
Float64 kiv;
Float64 kpvo,kpvs,kpvz,kpvm;
Float64 kivo,kivs,kivz,kivm;

Float64 a1,a1o,a1s,a2s;
Float64 a2,a2o;

Float64 kp_z;
Float64 ki_z;

Float64 ki_o;
Float64 kp_o;

Float64 ki_s;
Float64 kp_s;

/////////////////////////////////////////////////
Float64 po_ref_1[3]={0.0,0.0,0.0};
Float64 po_ref[3]={0.0,0.0,0.0};
Float64 err_ov[3]={0.0,0.0,0.0};
Float64 err_ov_1[3]={0.0,0.0,0.0};
Float64 xov[3]={0.0,0.0,0.0};
Float64 xov_1[3]={0.0,0.0,0.0};

Float64 pz_ref_1[3]={0.0,0.0,0.0};
Float64 pz_ref[3]={0.0,0.0,0.0};
Float64 err_zv[3]={0.0,0.0,0.0};
Float64 err_zv_1[3]={0.0,0.0,0.0};
Float64 xzv[3]={0.0,0.0,0.0};
Float64 xzv_1[3]={0.0,0.0,0.0};

Float64 ps_ref_1=0.0;
Float64 ps_ref=0.0;
Float64 err_sv=0.0;
Float64 err_sv_1=0.0;
Float64 xsv=0.0;
Float64 xsv_1=0.0;

Float64 pm_ref_1=0.0;
Float64 pm_ref=0.0;
Float64 err_mv=0.0;
Float64 err_mv_1=0.0;
Float64 xmv=0.0;
Float64 xmv_1=0.0;

////////////////////////////////////////////////
Float64 Voref[3]={0.0,0.0};
Float64 Voref_dq[2]={0.0,0.0};
Float64 Voref_dq_1[2]={0.0,0.0};
Float64 Vzref_dq[2]={0.0,0.0};
Float64 Vzref_dq_1[2]={0.0,0.0};

Float64 err_oi_dq_1[2]={0.0,0.0};
Float64 err_oi_dq[2]={0.0,0.0};
Float64 err_zi_dq_1[2]={0.0,0.0};
Float64 err_zi_dq[2]={0.0,0.0};
Float64 d_io_dq[2]={0.0,0.0};
Float64 d_iz_dq[2]={0.0,0.0};

Float64 Voref_1[3]={0.0,0.0,0.0};
Float64 err_oi[3]={0.0,0.0,0.0};
Float64 err_oi_1[3]={0.0,0.0,0.0};
Float64 xoi[2]={0.0,0.0};
Float64 xoi_1[2]={0.0,0.0};

Float64 Vzref[3]={0.0,0.0,0.0};
Float64 Vzref_1[3]={0.0,0.0,0.0};
Float64 err_zi[3]={0.0,0.0,0.0};
Float64 err_zi_1[3]={0.0,0.0,0.0};
Float64 xzi[3]={0.0,0.0,0.0};
Float64 xzi_1[3]={0.0,0.0,0.0};

Float64 Vsref=0.0;
Float64 Vsref_1=0.0;
Float64 err_si=0.0;
Float64 err_si_1=0.0;
Float64 xsi_1 = 0.0;
Float64 xsi = 0.0;

Float64 L = 0.015; 
Float64 Lac =0.015;
Float64 Ldc = 0.030;
Float64 Vdc = 162.5;
Float64 Vx;
Float64 r = 0.5;
Float64 rac = 0.2;
Float64 rdc = 0.8;
Float64 nc=2.0; //Celdas por rama
Float64 ro;
Float64 Lo;
Float64 rs;
Float64 Ls;
Float64 C = 0.001 ;
Float64 Rc = 20000.0;
Float64 Vc_ref = 75.0;
UInt32 Control_corrientes=0;
UInt32 Control_Voltajes=0;
UInt32 ready_volt = 0;
UInt32 ready_current = 0;
//Estados de proteccion
UInt32 over_voltage = 0;
UInt32 over_current = 0;
///////////////////////Saturaciones////////////////////
UInt32 tt[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 diff[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Prom[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Prom_1[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Prom_2[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Prom_3[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 comp[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Uprom=0.0;
Float64 k = 0.3;
Float64 DV [6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Rvolt_1[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Rvolt_2[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Rvolt_3[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Rvolt_4[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Rvolt_5[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Rvolt_6[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Float64 Bp[4]={0.0,0.0,0.0,0.0};
Float64 Ap[4]={0.0,0.0,0.0,0.0};

Float64 RvoltF[6]={0.0,0.0,0.0,0.0,0.0,0.0};


Float64 Io_max=8.0;
Float64 Iosat=10.0;
Float64 Umax_i;
Float64 Umax_v;
Float64 Umin_i;
Float64 Umin_v;
Float64 sinc;
Float64 sinc1;
Float64 sinc2,sinc3;
UInt32 count1=0;
UInt32 count2=0;
Float64 Vabc[3];
Float64 Valfa,Vbeta;



void FIR(Float64 *B,Float64 *In,Float64 *Out,int Order);
void ellip(Float64 *B,Float64 *A,Float64 *In,Float64 *Out,int Order);
void PI_back_dq(Float64 kp,Float64 ki,Float64 h,Float64 *in_1,Float64 *out_1, Float64 *in,Float64 *out,Float64 umax,Float64 umin);
void Sinc_ang(Float64 *Vabc,Float64 *Valp,Float64 *Vbeta,Float64 *ang);
void PI_tustin(Float64 kp,Float64 ki,Float64 h,Float64 *in_1,Float64 *out_1, Float64 *in,Float64 *out,Float64 umax,Float64 umin);
void PI_tustin_3ph(Float64 kp,Float64 ki,Float64 h,Float64 *in_1,Float64 *out_1, Float64 *in,Float64 *out,Float64 umax,Float64 umin);
void constantes();
void PI_back(Float64 kp,Float64 ki,Float64 h,Float64 *in_1,Float64 *out_1, Float64 *in,Float64 *out,Float64 umax,Float64 umin);
//matrix C = matrix A x matrix B , A(a_rows x a_cols), B(a_cols x b_cols)

void s_bal(Float64 *o,Float64 *z);
void PI_forward_AW(Float64 kp,Float64 ki,Float64 h,Float64 *x_1,Float64 *x,Float64 *out_1, Float64 *in,Float64 *out,Float64 umax,Float64 umin);
void PI_tustin_AW_3ph(Float64 kp,Float64 ki,Float64 h,Float64 *x_1,Float64 *x,Float64 *out_1, Float64 *in,Float64 *out,Float64 umax,Float64 umin);
void PI_forward_AW_3ph(Float64 kp,Float64 ki,Float64 h,Float64 *x_1,Float64 *x,Float64 *out_1, Float64 *in,Float64 *out,Float64 umax,Float64 umin);
void PI_tustin_AW(Float64 kp,Float64 ki,Float64 h,Float64 *x_1,Float64 *x,Float64 *out_1, Float64 *in,Float64 *out,Float64 umax,Float64 umin);
void abc2dq_pos(Float64 *abc, Float64 *dq,Float64 sinc);
void abc2dq_neg(Float64 *abc, Float64 *dq,Float64 sinc);
void dq2abc_pos(Float64 *dq,Float64 *abc,Float64 sinc);
void dq2abc_neg(Float64 *dq,Float64 *abc,Float64 sinc);
void dq2Vm(double *dq,double *m,double sinc);
void s_Dec2xy(Float64 *xy,Float64 *o,Float64 s,Float64 *z, Float64 m);
void s_xy2Dec(Float64 *xy,Float64 *o,Float64 *s,Float64 *z, Float64 *m);

#endif
