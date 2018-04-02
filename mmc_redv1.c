
#include "vicfuns.c"
#include <math.h>

 Int16 index1 = -1;
 Int16 task_id = 0; // comunication channel
 UInt16 chanels1[4]={1,2,3,4}; //Lectura CHANNELS 1,2,3,4 del conversor 1
 UInt16 chanels2[4]={5,6,7,8}; //Lectura CHANNELS 1,2,3,4 del conversor 1

 Float64 chmux1[4];
 Float64 chmux2[4];

 Float64 ch001, ch002, ch003, ch004, ch005, ch006, ch007, ch008, ch009, ch010;
Float64 IoR[3] = {0.0,0.0,0.0};
Float64 Gi = 1.0;
 UInt32 Offset=1;
 UInt32 bitmapB=0;
 UInt32 bitmap=0;
 Float64 exec_time_B=0;
 Float64 exec_time_A=0;
 Float64 k_offset1, k_offset2, k_offset3, k_offset4, k_offset5, k_offset6, k_offset7, k_offset8, k_offset9;//variac
 Float64 Offsetmed, Gain, x1, x2, x3, x4, x5, x6, x7, x8, x9; 
 UInt32 voltage_count=0;
 UInt32 current_count = 0;
Float64 Vz_mod[2]={0.0,0.0};
UInt32 testa;
 
 void controlMMC(void)
{

		  host_service(1,0);                          /* Data Acquisition service */
		  RTLIB_TIC_START();	                       /* Check interrupt overrun */
  /* ADC READING */
  // ds1103_adc_start(DS1103_ADC1 | DS1103_ADC2 | DS1103_ADC3 | DS1103_ADC5 | DS1103_ADC6 | DS1103_ADC7| DS1103_ADC8);
	ds1103_adc_start(DS1103_ADC1 | DS1103_ADC2 | DS1103_ADC5 | DS1103_ADC6 | DS1103_ADC7 | DS1103_ADC8);  
	ds1103_adc_read_mux(chanels1,4,chmux1);
	ds1103_adc_read_mux(chanels2,4,chmux2);

	if(Vdc < 1.0) Vdc=1.0;
	Vx=Vdc*0.5;
	if(Vc_ref< Vdc/4.0*1.2) 
		Vc_ref = Vdc/4.0*1.2;
	Exy_ref= Vc_ref*Vc_ref*C*0.5;
	Umax_i=Vc_ref*nc;
	Umax_v=Vx*Io_max;
	Umin_i=-Umax_i;
	Umin_v=-Umax_v;
	int i;
	bitmapB=ds1103_bit_io_read();

	readi =(bitmapB&0xff000000)>>24;
	
	
if(update_control){
	wn_v = 2.0*pi*fv;
wn_i = 2.0*pi*fi;//85
wn_io = 2.0*pi*fio;//85
wn_is = 2.0*pi*fis;//85


hv = Ts*Tv;		
 h =  Ts*Ti;	
m1=-2.0*exp(-xi*wn_v*hv)*cos(sqrt(1.0 -xi*xi)*wn_v*hv);
m2=exp(-2.0*xi*wn_v*hv);
////////////////////Con controlador utilizando aprox de tustin////////////////////
kpv =(2.0*exp(-2.0*hv/(C*Rc))+1.0+m1-m2)/(Rc*C*(1.0-exp(-2.0*hv/(C*Rc))));
kiv = 2.0*(m1+m2+1.0)/(Rc*C*hv*(1.0-exp(-2.0*hv/(C*Rc))));
kpvo = kpv;
kpvs = kpv;
kpvm = kpv;
kpvz = kpv;

kivo = kiv;
kivs = kiv;
kivz = kiv;
kivm = kiv;
//kiv=0.0;
a1=-2.0*exp(-xi*wn_i*h)*cos(sqrt(1.0-xi*xi)*wn_i*h);
a2=exp(-2.0*xi*wn_i*h);

a1o=-2.0*exp(-xi*wn_io*h)*cos(sqrt(1.0-xi*xi)*wn_io*h);
a2o=exp(-2.0*xi*wn_io*h);
kp_z=r*(a1-a2+1.0+2.0*exp(-h*r/L))/(2.0*(1.0-exp(-h*r/L)));
ki_z=r/h*(1.0+a1+a2)/(1.0-exp(-h*r/L));
//ki_z = 0.0;

ki_o=ro/h*(1.0+a1o+a2o)/(1.0-exp(-h*ro/Lo));
kp_o=ro*(a1o-a2o+1.0+2.0*exp(-h*ro/Lo))/(2.0*(1.0-exp(-h*ro/Lo)));
//ki_o = 0.0;
a1s=-2.0*exp(-xi*wn_is*h)*cos(sqrt(1.0-xi*xi)*wn_is*h);
a2s=exp(-2.0*xi*wn_is*h);

ki_s=rs/h*(1.0+a1s+a2s)/(1.0-exp(-h*rs/Ls));
kp_s=rs*(a1s-a2s+1.0+2.0*exp(-h*rs/Ls))/(2.0*(1.0-exp(-h*rs/Ls)));
	
	}
	//Uprom=(Prom[0]+Prom[1]+Prom[2]+Prom[3]+Prom[4]+Prom[5])/6.0;

	
	if(voltage_count>=Tv){
		voltage_count=0;
		ready_volt=1;
	}
	else{
		ready_volt=0;
	}
	voltage_count++;

	if(current_count>=Ti){
		current_count=0;
		ready_current=1;
	}
	else{
		ready_current=0;
		
	}
	current_count++;
switch(read_state){
	case 0: // espera de inicio de lectura
	{
		counter = 0;
	
		if(readi==0)
			read_state = 1;
		else read_state =0;
		break;
	}
	case 1:
	{ 
		
	if(counter<9){
		if(counter >= 7 & readi == 0xff){
			
			for(i=0;i<6;i++){
				
				if(voltage[i]==255 || voltage[i]==0){
					voltage_t[i]=Rvolt_2[i];
					testa=i;
					}
				else
					voltage_t[i]=(Float64)voltage[i]*0.784313;
				//if(voltage_t[i]<=4)voltage_t[i]=Prom[i];
				diff[i]= fabs(voltage_t[i]-Prom[i]);
				// if(diff[i]>9.0){
				// tt[i]=20;
				// voltage_t[i]=Prom[i];
				// }
				// else tt[i]=0;
				//comp[i]=Uprom-Prom[i];
				
				Prom[i] = voltage_t[i]*Bp[0]+Rvolt_1[i]*Bp[1]+Rvolt_2[i]*Bp[2]+Rvolt_3[i]*Bp[3]-(Prom_1[i]*Ap[0]+Prom_2[i]*Ap[1]+Prom_3[i]*Ap[2]);
				//Rvolt_6[i] = Rvolt_5[i];
				//Rvolt_5[i] = Rvolt_4[i];
				//Rvolt_4[i] = Rvolt_3[i];
				Rvolt_3[i] = Rvolt_2[i];
				Rvolt_2[i] = Rvolt_1[i];
				Rvolt_1[i]=voltage_t[i];
				Prom_3[i]=Prom_2[i];
				Prom_2[i]=Prom_1[i];
				Prom_1[i]=Prom[i];
					
				
				if(Prom[i]<2.0)Prom[i]=2.0;
				DV[i] = Prom[i]*k;
				if(diff[i]>DV[i]){
					RvoltF[i] = Prom[i];
					//voltage_t[i]=Prom[i];
					}
				else 
					RvoltF[i] = voltage_t[i];
				//RvoltF[i] = voltage_t[i];	
					if(RvoltF[i]>190.0)
							over_voltage = 1;
						else over_voltage = 0;
				//IIR orden 3 	
				// VcxyIn[i]=voltage_t[i];
				// VcxyO[i]=VcxyIn[i]*B1 +VcxyIn_1[i]*B2+VcxyIn_2[i]*B3- (VcxyO_1[i]*A2 + VcxyO_2[i]*A3);
				// VcxyIn_2[i]=VcxyIn_1[i];
				// VcxyIn_1[i]=VcxyIn[i];
				// VcxyO_2[i]=VcxyO_1[i];
				// VcxyO_1[i]=VcxyO[i];
				
				//orden 3 IIR [b,a]=ellip(3,.1,30,(200*180e-6)*2);
				// VcxyIn[i]=RvoltF[i];
				// VcxyO[i]=VcxyIn[i]*B1 +VcxyIn_1[i]*B2+VcxyIn_2[i]*B3+VcxyIn_3[i]*B4- (VcxyO_1[i]*A2 + VcxyO_2[i]*A3+VcxyO_3[i]*A4);
				
				// VcxyIn_3[i]=VcxyIn_2[i];
				// VcxyIn_2[i]=VcxyIn_1[i];
				// VcxyIn_1[i]=VcxyIn[i];
			
				
				// VcxyO_3[i]=VcxyO_2[i];
				// VcxyO_2[i]=VcxyO_1[i];
				// VcxyO_1[i]=VcxyO[i];
				
				
				
				
				// VcxyIn[i]=RvoltF[i];
				// VcxyO[i]=VcxyIn[i]*B1 +VcxyIn_1[i]*B2+VcxyIn_2[i]*B3+VcxyIn_3[i]*B4+VcxyIn_4[i]*B5+VcxyIn_5[i]*B6+VcxyIn_6[i]*B7- (VcxyO_1[i]*A2 + VcxyO_2[i]*A3+VcxyO_3[i]*A4 + VcxyO_4[i]*A5+VcxyO_5[i]*A6+VcxyO_6[i]*A7);
				// VcxyIn_6[i]=VcxyIn_5[i];
				// VcxyIn_5[i]=VcxyIn_4[i];				
				// VcxyIn_4[i]=VcxyIn_3[i];
				// VcxyIn_3[i]=VcxyIn_2[i];
				// VcxyIn_2[i]=VcxyIn_1[i];
				// VcxyIn_1[i]=VcxyIn[i];
			
				// VcxyO_6[i]=VcxyO_5[i];			
				// VcxyO_5[i]=VcxyO_4[i];
				// VcxyO_4[i]=VcxyO_3[i];
				// VcxyO_3[i]=VcxyO_2[i];
				// VcxyO_2[i]=VcxyO_1[i];
				// VcxyO_1[i]=VcxyO[i];
				
				//Filtro FIR orden 5
				Vcfir[i] = RvoltF[i];
				VcxyO[i]=Vcfir[i]*F[0]+Vcfir_1[i]*F[1]+Vcfir_2[i]*F[2]+Vcfir_3[i]*F[3]+Vcfir_4[i]*F[4]+Vcfir_5[i]*F[5];
				
				Vcfir_5[i]=Vcfir_4[i];
				Vcfir_4[i]=Vcfir_3[i];
				Vcfir_3[i]=Vcfir_2[i];
				Vcfir_2[i]=Vcfir_1[i];
				Vcfir_1[i]=Vcfir[i];
				
				
				
				//Filtro FIR orden 8
				// Vcfir[i] = RvoltF[i];
				// Vcfir[i] = RvoltF[i];
				// VcxyO[i]=Vcfir[i]*F[0]+Vcfir_1[i]*F[1]+Vcfir_2[i]*F[2]+Vcfir_3[i]*F[3]+Vcfir_4[i]*F[4]+Vcfir_5[i]*F[5]+Vcfir_6[i]*F[6]+Vcfir_7[i]*F[7]+Vcfir_8[i]*F[8];
				// Vcfir_8[i]=Vcfir_7[i];
				// Vcfir_7[i]=Vcfir_6[i];
				// Vcfir_6[i]=Vcfir_5[i];
				// Vcfir_5[i]=Vcfir_4[i];
				// Vcfir_4[i]=Vcfir_3[i];
				// Vcfir_3[i]=Vcfir_2[i];
				// Vcfir_2[i]=Vcfir_1[i];
				// Vcfir_1[i]=Vcfir[i];
				
				}
				read_state = 0; //lectura correcta
				break;
			}
		else if(counter==8&readi!=0xff){
			//error de lectura
			read_state = 0;
			
			error_counter++;
			break;
			}
		
		else {
			
			read_state = 1;
			voltage[counter]=readi;
			counter++;
			break;
			}
		}
		else read_state = 0;
		break;
	}
	
	default:{read_state=0;
		counter = 0;
		for(i=0;i<6;i++)
				voltage_t[i]=0;
				break;
	
	}

}

 

  /* channel assignation */
  ch001=(ds1103_adc_read_ch(17)*10.0);
  ch002=(ds1103_adc_read_ch(18)*10.0);
  ch003=(ds1103_adc_read_ch(19)*10.0);
  ch004=(ds1103_adc_read_ch(20)*10.0);
  ch005=(chmux1[0]*10.0);
  ch006=(chmux1[1]*10.0);
  ch007=(chmux2[0]*10.0)*200.0;//Va
  ch008=(chmux2[1]*10.0)*200.0;//Vb
  ch009=(chmux2[2]*10.0)*200.0;//Vb

  
  if (Offset) {
		x1 = (1-Ts*0.2)*x1 + Ts*0.2*ch001;
		x2 = (1-Ts*0.2)*x2 + Ts*0.2*ch002;
		x3 = (1-Ts*0.2)*x3 + Ts*0.2*ch003;
		x4 = (1-Ts*0.2)*x4 + Ts*0.2*ch004;
		x5 = (1-Ts*0.2)*x5 + Ts*0.2*ch005;
		x6 = (1-Ts*0.2)*x6 + Ts*0.2*ch006;
		x7 = (1-Ts*0.2)*x7 + Ts*0.2*ch007;
		x8 = (1-Ts*0.2)*x8 + Ts*0.2*ch008;
		x9 = (1-Ts*0.2)*x9 + Ts*0.2*ch009;
	}
	else{
		k_offset1=x1;
		k_offset2=x2;
		k_offset3=x3;
		k_offset4=x4;
		k_offset5=x5;
		k_offset6=x6;
		k_offset7=x7;
		k_offset8=x8;
		k_offset9=x9;
	}
	
	
	Vabc[0]=(ch007-k_offset7)*kgain7;
	Vabc[1]=(ch008-k_offset8)*kgain8;
	Vabc[2]=(ch009-k_offset9)*kgain9;
	//Sinc_ang(Vabc,&Valfa,&Vbeta,&sinc);
	
	//////////////////////////////////////////////////////////////////AQUI CORRIENTES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
  
    
//Filtro de corrientes orden 3 pasabanda 100 Hz
  
if(ready_volt){
  vcpa=VcxyO[0]-vcpa_off;
  vcpb=VcxyO[1]-vcpb_off;
  vcpc=VcxyO[2]-vcpc_off;
  vcna=VcxyO[3]-vcna_off;
  vcnb=VcxyO[4]-vcnb_off;
  vcnc=VcxyO[5]-vcnc_off;

////////////////////////////////////////////////////////////
////													////
////				Control de Tension					////
////								                    ////
////////////////////////////////////////////////////////////

  Exy[0]=vcpa*vcpa*C*0.5;
  Exy[1]=vcpb*vcpb*C*0.5;
  Exy[2]=vcpc*vcpc*C*0.5;
  Exy[3]=vcna*vcna*C*0.5;
  Exy[4]=vcnb*vcnb*C*0.5;
  Exy[5]=vcnc*vcnc*C*0.5;

  

	
  s_xy2Dec(Exy,Eo,&Es,Ez,&Em);

  err_ov[0]=-Eo[0];
  err_ov[1]=-Eo[1];
  err_ov[2]=-Eo[2];
 
  err_zv[0]=-Ez[0];
  err_zv[1]=-Ez[1];
  err_zv[2]=-Ez[2];
  
  err_sv = -Es;
  err_mv =Exy_ref-Em;
////////////////////////////////////////////////////////////////////////////
  if(Control_maximus&&!protect2){
  //Control de Eo
  //(Float64 kp,Float64 ki,Float64 h,Float64 *in_1,Float64 *out_1, Float64 *in,Float64 *out,Float64 umax,Float64 umin)
  //PI_tustin_3ph(kpv,kiv,hv,err_ov_1,po_ref_1,err_ov,po_ref,Umax_v,Umin_v);
  //(Float64 kp,Float64 ki,Float64 h,Float64 *x_1,Float64 *x,Float64 *out_1, Float64 *in,Float64 *out,Float64 umax,Float64 umin)
  // PI_tustin_AW_3ph(kpv,kiv,hv,xov_1,xov,po_ref_1,err_ov,po_ref,Umax_v,Umin_v);
  //Control de Es
   //PI_tustin_AW(kpv,kiv,hv,&xsv_1,&xsv,&ps_ref_1,&err_sv,&ps_ref,Umax_v,Umin_v);
  //Control de Ez
  // PI_tustin_AW_3ph(kpv,kiv,hv,xzv_1,xzv,pz_ref_1,err_zv,pz_ref,Umax_v,Umin_v);
  //Control de Em
  // PI_tustin_AW(kpv,kiv,hv,&xmv_1,&xmv,&pm_ref_1,&err_mv,&pm_ref,Umax_v,Umin_v);
 
  PI_tustin_3ph(kpvo,kivo,hv,err_ov_1,po_ref_1,err_ov,po_ref,Umax_v,Umin_v);
  
  PI_tustin(kpvs,kivs,hv,&err_sv_1,&ps_ref_1,&err_sv,&ps_ref,Umax_v,Umin_v);
  
  PI_tustin_3ph(kpvz,kivz,hv,err_zv_1,pz_ref_1,err_zv,pz_ref,Umax_v,Umin_v);

  PI_tustin(kpvm,kivm,hv,&err_mv_1,&pm_ref_1,&err_mv,&pm_ref,Umax_v,Umin_v);

  }
  else{
  err_ov_1[0]=0.0;
  err_ov_1[1]=0.0;
  err_ov_1[2]=0.0;
 
  err_zv_1[0]=0.0;
  err_zv_1[1]=0.0;
  err_zv_1[2]=0.0;
  
  err_sv_1=0.0;
  xsv_1 = 0.0;
  
  err_mv_1=0.0;
    
  po_ref[0]=0.0;
  po_ref[1]=0.0;
  po_ref[2]=0.0;
 
  pz_ref[0]=0.0;
  pz_ref[1]=0.0;
  pz_ref[2]=0.0;
  
  ps_ref=0.0;
  pm_ref=0.0;
  
  po_ref_1[0]=0.0;
  po_ref_1[1]=0.0;
  po_ref_1[2]=0.0;
 
  pz_ref_1[0]=0.0;
  pz_ref_1[1]=0.0;
  pz_ref_1[2]=0.0;
  
  ps_ref_1=0.0;
  pm_ref_1=0.0;
  }
}

if(ready_current){
 if(!over_current){
  Ixy[0]=(ch001-k_offset1)*kgain1; //iap
  Ixy[1]=(ch002-k_offset2)*kgain2; //ibp
  Ixy[2]=(ch003-k_offset3)*kgain3; //icp
  Ixy[3]=-1.0*(ch004-k_offset4)*kgain4; //ian
  Ixy[4]=-1.0*(ch005-k_offset5)*kgain5; //ibn
  Ixy[5]=-1.0*(ch006-k_offset6)*kgain6; //icn
  }
  else{
	Ixy[0]=Ixy[0];
	Ixy[1]=Ixy[1];
	Ixy[2]=Ixy[2];
	Ixy[3]=Ixy[3];
	Ixy[4]=Ixy[4];
	Ixy[5]=Ixy[5];	
  }
  
  //Filtro de corrientes orden 3 pasabanda 100 Hz
  for(i=0;i<6;i++){
  IxyO[i]=Ixy[i];
		
  		// IxyO[i]=Ixy[i]*Fii[0]+IxyF_1[i]*Fii[1]+IxyF_2[i]*Fii[2]+IxyF_3[i]*Fii[3]+IxyF_4[i]*Fii[4]+IxyF_5[i]*Fii[5];//+IxyF_6[i]*Fii[6];
		
		//IxyF_6[i]=IxyF_5[i];
		// IxyF_5[i]=IxyF_4[i];
		// IxyF_4[i]=IxyF_3[i];
		// IxyF_3[i]=IxyF_2[i];
		// IxyF_2[i]=IxyF_1[i];
		// IxyF_1[i]=Ixy[i];		
		
		// IxyF[i]=Ixy[i];
		// IxyO[i]=IxyF[i]*B1i +IxyF_1[i]*B2i+IxyF_2[i]*B3i+IxyF_3[i]*B4i- (IxyO_1[i]*A2i + IxyO_2[i]*A3i+IxyO_3[i]*A4i);
		
		// IxyF_3[i]=IxyF_2[i];
		// IxyF_2[i]=IxyF_1[i];
		// IxyF_1[i]=IxyF[i];
	
		// IxyO_3[i]=IxyO_2[i];
		// IxyO_2[i]=IxyO_1[i];
		// IxyO_1[i]=IxyO[i]; 
}
  
  if(fabs(IxyO[0])>Iosat)
	over_current = 1;
  if(fabs(IxyO[1])>Iosat)
	over_current = 1;
  if(fabs(IxyO[2])>Iosat)
	over_current = 1;
  if(fabs(IxyO[3])>Iosat)
	over_current = 1;
  if(fabs(IxyO[4])>Iosat)
	over_current = 1;
  if(fabs(IxyO[5])>Iosat)
	over_current = 1;  
	
  protect = over_current||over_voltage;
	if(protect)
		protect2=1;
  //0 para corriente positiva
  if(IxyO[0]>=0)
	signpa = 0;
  else	
	signpa = 1;

  if(IxyO[1]>=0)
	signpb = 0;
  else	
	signpb = 1;
	
  if(IxyO[2]>=0)
	signpc = 0;
  else	
	signpc = 1;
	//Invierto signo para control rama negativa
  if(IxyO[3]>=0)
	signna = 1;
  else	
	signna = 0;
	
  if(IxyO[4]>=0)
	signnb = 1;
  else	
	signnb = 0;
	
  if(IxyO[5]>=0)
	signnc = 1;
  else	
	signnc = 0;
	
  // signpa=(Ixy[0]>=0) ? 0 : 1;
  // signpb=(Ixy[1]>=0) ? 0 : 1;
  // signpc=(Ixy[2]>=0) ? 0 : 1;
  // signna=(Ixy[3]>=0) ? 0 : 1;
  // signnb=(Ixy[4]>=0) ? 0 : 1;
  // signnc=(Ixy[5]>=0) ? 0 : 1; 
  for(i=0;i<6;i++)
	Pxy[i]=Mxy[i]*Ixy[i];
	
	
	s_xy2Dec(Pxy,Po,&Ps,Pz,&Pm);
	
	// Poa[0] = Po[0];
	// Pob[0] = Po[1];
	// Poc[0] = Po[2];
	
	// ellip(Bpo2,Apo2,Poa,PoaF,Po_order_BP);
	// ellip(Bpo2,Apo2,Pob,PobF,Po_order_BP);
	// ellip(Bpo2,Apo2,Poc,PocF,Po_order_BP);
	
	Poabc[0] = Po[0];
	Poabc[1] = Po[1];
	Poabc[2] = Po[2];

	abc2dq_neg(Poabc,Podq,sinc2);//revisar
	abc2dq_pos(Pz,Pzdq,sinc1);//revisar

	Pod[0]=Podq[0];
	Poq[0]=Podq[1];

	
	ellip(Bpo,Apo,Pod,PodF,Po_order);
	ellip(Bpo,Apo,Poq,PoqF,Po_order);
	
	PodqF[0] = PodF[0];
	PodqF[1] = PoqF[0];
	
	Pzd[0]=Pzdq[0];
	Pzq[0]=Pzdq[1];
	
	ellip(Bpo,Apo,Pzd,PzdF,Pz_order);
	ellip(Bpo,Apo,Pzq,PzqF,Pz_order);
	
	PzdqF[0] = PzdF[0];
	PzdqF[1] = PzqF[0];	
	
	
	
  s_xy2Dec(IxyO,Io,&Is,Iz,&Im);
  IoR[0]=Io[0]*Gi;
  IoR[1]=Io[1]*Gi;
  IoR[1]=Io[2]*Gi;
  
	abc2dq_pos(Vabc,Vdq,ang);
	Vq=Vdq[1];
	PI_tustin(kp_pll,ki_pll,h,&Vq_1,&Wred_1,&Vq,&Wred,Umax_pll,Umin_pll);
    Wred = Wred + 2.0*pi*50.0; 
	ang = Wred*h + ang_1;
	ang=fmod(ang,2.0*pi);
  if(ang >= 2.0*pi) ang = ang - 2.0*pi;
  else ang = ang;
  ang_1 = ang;
  sinc1=ang;
  
  if(ang<=pi)
	sinc2=(Float64)2.0*ang; 
  else
	sinc2=(Float64)2.0*ang-2.0*pi;
	
  if(ang<=_2pi3)
	sinc3=(Float64)3.0*ang;
	else if(ang>_2pi3 & ang<= _4pi3)
	sinc3=(Float64)3.0*ang-2.0*pi;
	else 
	sinc3=(Float64)3.0*ang-4.0*pi;
	
  d_izref[0]=po_ref[0]/Vx;
  d_izref[1]=po_ref[1]/Vx;
  d_izref[2]=po_ref[2]/Vx;
	
  d_ioref[0]=pz_ref[0]/Vx;
  d_ioref[1]=pz_ref[1]/Vx;
  d_ioref[2]=pz_ref[2]/Vx;
  
  
if(inyect_IZ)
 {
	err_iz[0] =pzd_ref-PzdqF[0];// pod_ref-PodqF[0];
	err_iz[1] =pzq_ref-PzdqF[1];// poq_ref-PodqF[1];

	PI_tustin_dq(kp_iz,ki_iz,h,err_iz_1,Out_iz_1,err_iz,Out_iz,max_iz,min_iz);
	
	Izdq_ref[0] = Out_iz[0]/Vdc;
	Izdq_ref[1] = Out_iz[1]/Vdc;

}
else
{
	err_iz_1[0] = 0.0;
	err_iz_1[1] = 0.0;
	Out_iz_1[0] = 0.0;
	Out_iz_1[1] = 0.0;	
	Out_iz[0] = 0.0;	
	Out_iz[1] = 0.0;	
	

}	
  
if(inyect_VM)
 {
	err_vm[0] =pod_ref-PodqF[0];// pzd_ref-PzdqF[0];
	err_vm[1] =poq_ref-PodqF[1];// pzq_ref-PzdqF[1];

	PI_tustin_dq(kp_vm,ki_vm,h,err_vm_1,Out_vm_1,err_vm,Out_vm,max_vm,min_vm);
	
	Vmdq_iny[0] = Out_vm[0];
	Vmdq_iny[1] = Out_vm[1];

}
else
{
	err_vm_1[0] = 0.0;
	err_vm_1[1] = 0.0;
	Out_vm_1[0] = 0.0;
	Out_vm_1[1] = 0.0;	
	Out_vm[0] = 0.0;	
	Out_vm[1] = 0.0;	
	//Vmdq_iny[0] =0.0;
	//Vmdq_iny[1] = 0.0;
}
  
  d_isref = pm_ref/Vx;
  if(d_isref>0)is_sign = 1.0;
	else is_sign=-1.0;
	
  // if(fabs(d_isref)<0.01)Is_est=0.01;
  // else Is_est = fabs(d_isref);
  //revisar
  //Is_est=1.0;
  d_vmref = -1.0*ps_ref/Is_est*is_sign; //Dividir por idc...calcular
    
	
  abc2dq_pos(d_ioref,d_io_dq,sinc1);
  abc2dq_pos(Io,IodqT,sinc1);
  // dq2abc_pos(Iodq_ref,Io_ref,sinc1);
  // Iod[0]= IodqT[0];
  // Ioq[0]=IodqT[1];
 // ellip(Bio,Aio,Iod,IodF,3);
  // ellip(Bio,Aio,Ioq,IoqF,3);
// FIR(BioF,Iod,&IodFir,5);
  // FIR(BioF,Ioq,&IoqFir,5);
  
  Iodq[0] = IodqT[0];
  Iodq[1] = IodqT[1];
  Io_refT[0]=Iodq_ref[0]+d_io_dq[0];
  Io_refT[1]=Iodq_ref[1]+d_io_dq[1];
  
  // Io_refT[0]=Io_ref[0]+d_ioref[0];
  // Io_refT[1]=Io_ref[1]+d_ioref[1];
  // Io_refT[2]=Io_ref[2]+d_ioref[2];
  
  // err_oi[0]=Io_refT[0]-Io[0];
  // err_oi[1]=Io_refT[1]-Io[1];
  // err_oi[2]=Io_refT[2]-Io[2];
  
  err_oi_dq[0]=Io_refT[0]-Iodq[0];
  err_oi_dq[1]=Io_refT[1]-Iodq[1];
  
  
  abc2dq_neg(d_izref,d_iz_dq,sinc2);
  abc2dq_neg(Iz,Izdq,sinc2);
  // dq2abc_neg(Izdq_ref,Iz_ref,sinc2);
  Iz_refT[0]=Izdq_ref[0]+d_iz_dq[0];
  Iz_refT[1]=Izdq_ref[1]+d_iz_dq[1];
  
  err_zi_dq[0]=Iz_refT[0]-Izdq[0];
  err_zi_dq[1]=Iz_refT[1]-Izdq[1];
  

  // Iz_refT[0]=Iz_ref[0]+d_izref[0];
  // Iz_refT[1]=Iz_ref[1]+d_izref[1];
  // Iz_refT[2]=Iz_ref[2]+d_izref[2];
  
  // err_zi[0]=Iz_refT[0]-Iz[0];/////*********
  // err_zi[1]=Iz_refT[1]-Iz[1];
  // err_zi[2]=Iz_refT[2]-Iz[2];
  
  err_si=d_isref-Is;
  
  
/////////////////////////////////////////////////////
////											 ////
////			Control de corrientes			 ////
////											 ////
/////////////////////////////////////////////////////



	  if(Control_maximus&&!protect2){
	 // PI_tustin_AW(kpv,kiv,hv,&xsv_1,&xsv,&ps_ref_1,&err_sv,&ps_ref,Umax_v,Umin_v);
	 //void PI_tustin_dq_AW(Float64 kp,Float64 ki,Float64 h,Float64 *x_1,Float64 *x,Float64 *out_1, Float64 *in,Float64 *out,Float64 umax,Float64 umin)

	  //Control de Io  (Float64 kp,Float64 ki,Float64 h,Float64 *in_1,Float64 *out_1, Float64 *in,Float64 *out,Float64 umax,Float64 umin)
	  //PI_tustin_3ph(kp_o,ki_o,h,err_oi_1,Voref_1,err_oi,Voref,Umax_i,Umin_i);
	  
	  PI_tustin_dq(kp_o,ki_o,h,err_oi_dq_1,Voref_dq_1,err_oi_dq,Voref_dq,Umax_i,Umin_i);
	  //(Float64 kp,Float64 ki,Float64 h,Float64 *in_1,Float64 *out_1, Float64 *in,Float64 *out,Float64 umax,Float64 umin)
	  //PI_back_dq(kp_o,ki_o,h,err_oi_dq_1,Voref_dq_1,err_oi_dq,Voref_dq,Umax_i,Umin_i);
	  
	  //PI_tustin_dq_AW(kp_o,ki_o,h,xoi_1,xoi,Voref_dq_1,err_oi_dq,Voref_dq,Umax_i,Umin_i);
	  //Control de Is
	  PI_tustin(kp_s,ki_s,h,&err_si_1,&Vsref_1,&err_si,&Vsref,Umax_i,Umin_i);
	  //Control de Iz
	  //PI_tustin_3ph(kp_z,ki_z,h,err_zi_1,Vzref_1,err_zi,Vzref,Umax_i,Umin_i);
	  PI_tustin_dq(kp_z,ki_z,h,err_zi_dq_1,Vzref_dq_1,err_zi_dq,Vzref_dq,Umax_i,Umin_i);
	  dq2Vm(Vmdq_iny,&Vm_ripple,sinc3);
	  Vmref=d_vmref+Vm_ripple;
	  //dq2abc_pos(Vo_mod,Voref,sinc1);
	  //dq2abc_neg(Vz_mod,Vzref,sinc2);
	  }
	  else{
	  
	  err_si_1=0.0; 
	  Vmref=0.0;//+Vm_ripple;
	  Vm_ripple=0.0;
	  err_oi_dq_1[0]=0.0;
	  err_oi_dq_1[1]=0.0;
	  Voref_dq[0]=0.0;
	  Voref_dq[1]=0.0;
	  Voref_dq_1[0]=0.0;
	  Voref_dq_1[1]=0.0;
	  // Voref[0]=0.0;
	  // Voref[1]=0.0;
	  // Voref[2]=0.0;
	  // Voref_1[0]=0.0;
	  // Voref_1[1]=0.0;
	  // Voref_1[2]=0.0;
	  // err_oi_1[0]=0.0;
  	  // err_oi_1[1]=0.0;
	  // err_oi_1[2]=0.0;
	  
	  err_zi_dq_1[0]=0.0;
	  err_zi_dq_1[1]=0.0;
	  Vzref_dq[0]=0.0;
	  Vzref_dq[1]=0.0;
	  Vzref_dq_1[0]=0.0;
	  Vzref_dq_1[1]=0.0;

	  // err_zi_1[0]=0.0;
  	  // err_zi_1[1]=0.0;
	  // err_zi_1[2]=0.0;
	  // Vzref[0]=0.0;
	  // Vzref[1]=0.0;
	  // Vzref[2]=0.0;
	  // Vzref_1[0]=0.0;
	  // Vzref_1[1]=0.0;
	  // Vzref_1[2]=0.0;
	  
	  // xoi_1[0] = 0.0;
	  // xoi_1[1] = 0.0;
	  // xoi_1[2] = 0.0;

	  // xzi_1[0] = 0.0;
	  // xzi_1[1] = 0.0;
	  // xzi_1[2] = 0.0;
	  
	  // xsv_1 = 0.0;
	  // xsv = 0.0;

	  Vsref=0.0;
	  
	  Vsref_1=0.0;
	 
	/////////////////////////////////////////////////////
		}
		dq2abc_pos(Voref_dq,Voref,sinc1);
		dq2abc_neg(Vzref_dq,Vzref,sinc2);

	//Vsref= 60.0;
	//Vmref=40.0;	
	
		/////////////////////////////////////////////////////
		s_Dec2xy(Vxy,Voref,Vsref,Vzref,Vmref);
		
	  for(i=0;i<3;i++){
	   Mxy[i]=Vx-Vabc[i]-Vxy[i];
	   if(Mxy[i]>Umax_i)
			Mxy[i]=Umax_i;
	   else if(Mxy[i]<Umin_i)
			Mxy[i]=Umin_i;
	   else
			Mxy[i]=Vx-Vabc[i]-Vxy[i];
			
	   //Mxy[i]=Vx-Vabc[i]-Vxy[i];
	   
	   }
	  for(i=3;i<6;i++){
	  Mxy[i]=-Vx-Vabc[i-3]-Vxy[i]; 
	   if(Mxy[i]>Umax_i)
			Mxy[i]=Umax_i;
	   else if(Mxy[i]<Umin_i)
			Mxy[i]=Umin_i;
	   else
			Mxy[i]=-Vx-Vabc[i-3]-Vxy[i]; 
	   //Mxy[i]=-Vx-Vabc[i-3]-Vxy[i]; 
	   }
	}

modxy[0]=Mxy[0]/Vc_ref/nc;
modxy[1]=Mxy[1]/Vc_ref/nc;
modxy[2]=Mxy[2]/Vc_ref/nc;
modxy[3]=Mxy[3]/Vc_ref/nc;
modxy[4]=Mxy[4]/Vc_ref/nc;
modxy[5]=Mxy[5]/Vc_ref/nc;

if(!Control_maximus||protect2)  
{
  levelAp = 4;//revisar 4 genera 2vdc todas las celdas abiertas condensadores se cargan a Vdc
  levelBp = 4;
  levelCp = 4;
  levelAn = 4;
  levelBn = 4;
  levelCn = 4;
}
else {
  levelAp = floor((Mxy[0]/vcpa)+2.0+0.5);//revisar
  levelBp = floor((Mxy[1]/vcpb)+2.0+0.5);
  levelCp = floor((Mxy[2]/vcpc)+2.0+0.5);
  levelAn = floor((Mxy[3]/-vcna)+2.0+0.5);
  levelBn = floor((Mxy[4]/-vcnb)+2.0+0.5);
  levelCn = floor((Mxy[5]/-vcnc)+2.0+0.5);
	if(levelAp>4)levelAp = 4;
	else if(levelAp<0)levelAp=0;
	
	if(levelBp>4)levelBp = 4;
	else if(levelBp<0)levelBp=0;
	
	if(levelCp>4)levelCp = 4;
	else if(levelCp<0)levelCp=0;
	
	if(levelAn>4)levelAn = 4;
	else if(levelAn<0)levelAn=0;
		
	if(levelBn>4)levelBn = 4;
	else if(levelBn<0)levelBn=0;
	
	if(levelCn>4)levelCn = 4;
	else if(levelCn<0)levelCn=0;
	///////corto circuito, evitar este estado, la fpga no lo puede manejar
	if(levelAp==2&&levelAn==2){
		warn_cero++;
		//levelAp=3;//vdc
		//levelAn=1;//-vdc
	}
	if(levelBp==2&&levelBn==2){
		warn_cero++;
		//levelBp=3;//vdc
		//levelBn=1;//-vdc
	}
	if(levelCp==2&&levelCn==2){
		warn_cero++;
		//levelCp=3;//vdc
		//levelCn=1;//-vdc
	}
	//////////////////////////////////
	
}

L3=(levelCn<<1)+signnc;
L2=(levelBn<<1)+signnb;
L1=(levelAn<<1)+signna;

U3=(levelCp<<1)+signpc;
U2=(levelBp<<1)+signpb;
U1=(levelAp<<1)+signpa;

fpga = (L3<<20)+(L2<<16)+(L1<<12)+(U3<<8)+(U2<<4)+U1;
//fpga=0x00_II_LL_UU;
ds1103_bit_io_write(fpga); 

exec_time_A = RTLIB_TIC_READ()*1e6;
RTLIB_SRT_ISR_END(); /* overload check */
}	
/////////////////////////////////////////////////////////////////////////
////																 ////
////				Timer_B lectura de voltage						 ////
////																 ////
/////////////////////////////////////////////////////////////////////////	
/*
void voltage_read(void){
	ds1103_begin_isr_timerB();
	RTLIB_TIC_START();
	host_service(1,0);  Background service 
	int i;
	bitmapB=ds1103_bit_io_read();

	readi =(bitmapB&0xff000000)>>24;

switch(read_state){
	case 0: // espera de inicio de lectura
	{
		counter = 0;
		if(readi==0)
			read_state = 1;
		else read_state =0;
			break;
	}
	case 1:
	{ if(counter<9){
		if(counter >= 7 & readi == 0xff){
			read_state = 0; //lectura correcta
			for(i=0;i<6;i++){
				if(voltage[i]==0 || voltage[i]==255)
					voltage_t[i]=VcxyIn_1[i]; //lectura anterior
				else
					voltage_t[i]=voltage[i]*0.784313;
				VcxyIn[i]=voltage_t[i];
				VcxyO[i]=VcxyIn[i]*B1 +VcxyIn_1[i]*B2+VcxyIn_2[i]*B3- VcxyO_1[i]*A2 - VcxyO_2[i]*A3;
				VcxyIn_2[i]=VcxyIn_1[i];
				VcxyIn_1[i]=VcxyIn[i];
				VcxyO_2[i]=VcxyO_1[i];
				VcxyO_1[i]=VcxyO[i];
				}
			}
		else if(counter==8&readi!=0xff){
			//error de lectura
			read_state = 0;
			error_counter++;
			}
		
		else {
			read_state = 1;
			voltage[counter]=readi;
			counter ++;
			}
		}
		else read_state = 0;
		break;
	}
	
	default:{read_state=0;
		counter = 0;
		for(i=0;i<6;i++)
				voltage_t[i]=0;
	
	}

}
	exec_time_B = RTLIB_TIC_READ();
    ds1103_end_isr_timerB();
	//RTLIB_SRT_ISR_END();  //overload check 
}
	*/
/////////////////////////////////////////////////////////
////												 ////
////					Main						 ////
////												 ////
/////////////////////////////////////////////////////////	
	
void main(void)
{
	init();										/* DS1103 and RTLib1103 initialization */
	constantes();
	ds1103_bit_io_init(DS1103_DIO1_OUT|DS1103_DIO2_OUT|DS1103_DIO3_OUT|DS1103_DIO4_IN);
		/*Initialize digital outputs*/
    ds1103_bit_io_write(0x00000000);
	
	/////////////////////////////////////////	
	ds1103_dac_init(DS1103_DACMODE_LATCHED);
	RTLIB_SRT_START(Ts, controlMMC);				
	//ds1103_start_isr_timerB(0,TsB,voltage_read);   
	
	while(1)									/* Background task */
	{
		RTLIB_BACKGROUND_SERVICE();
	}
}
