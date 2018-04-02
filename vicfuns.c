
#include "vicfuns.h"


#define NR_END 1
#define FREE_ARG char*
#define lrow 2
#define rcol 3

void constantes()
{
ro=2.0*rac+r;
Lo=2.0*Lac+L;
rs=3.0*rac+r;
Ls=3.0*Ldc+L;
//Parametros filtro de voltaje, Se deben modificar
//500 Hz
// B1 =   0.093287236976150156309017802414019;
// B2 =  0.10042883416963054032944313576081;
// B3 = B1;
// A2 =  -1.2010706569493829753980662644608;
// A3 =   0.52309366510790211712844666180899;
//300 Hz
// B1 = 0.039140943738022432318857113386912;
// B2 = -0.032983642108771304979342176011414;
// B3 = B1;
// A2 = -1.733849759140064783480283949757;
// A3 = 0.78467522638964415371276572841452;

//iir orden 6

// B1 = 0.032155054388409695653727737862937;
// B2 =  -0.13033718931674126428887916517851;
// B3 =  0.25943835555117727986740305823332;
// B4= -0.31911558993389427252296286496858;
// B5=  0.25943835555117727986740305823332;
// B6= -0.13033718931674126428887916517851;
// B7=  0.03215505438840968177593993004848;

// A2 =  -4.7870600537391894135907932650298;
// A3 =  9.843800960288552559518393536564;
// A4=-11.072719449563631854971390566789;
// A5=7.1672826109763967394883366068825;
// A6=-2.5267977017353904223284644103842;
// A7= 0.37889439810674097053322384454077;
// ellip orden 3 Fc 20 Hz

// B1 =  0.001667379966876;
// B2 = -0.001664125414121;
// B3 =B2;
// B4=   0.001667379966876;

// A2 =  -2.977661574967806;
// A3 =    2.955963422860258;
// A4= -0.978295338786943;

// ellip orden 3 Fc 80 Hz

// B1 =  0.00939195455596;
// B2 =   -0.008796357408537;
// B3 =   -0.008796357408537;
// B4=   0.009391954555966;

// A2 =  -2.822707027077662;
// A3 =   2.666277979368142;
// A4= -0.842379757995621;

// ellip orden 3 Fc 100 Hz

B1 = 0.018964786512529;
B2 =-0.018432508824266;
B3 =-0.018432508824266;
B4= 0.018964786512529;

A2 = -2.836495701715605;
A3 = 2.690687661199667;
A4= -0.853127404107538;


//Fir orden 5 freq de corte 100 hz b=fir1(5,(100*9*25e-6)*2)

F[0]= 0.028291307304640;
F[1] =0.142595018408818;
F[2] =0.329113674286542;
F[3] =0.329113674286542;
F[4]=0.142595018408818;
F[5]= 0.028291307304640;

//Fir orden 5 freq de corte 50 hz b=fir1(5,(50*9*20e-6)*2)

// F[0]= 0.028699070725939;
// F[1] =0.143029716046051;
// F[2] =0.328271213228010;
// F[3] =0.328271213228010;
// F[4]=0.143029716046051;
// F[5]= 0.028699070725939;

//Fir orden 8 freq de corte 300 hz

// F[0]=   0.017792743734167;
// F[1] = 0.048311483314256;
// F[2] =0.122492362235122;
// F[3] = 0.197240676082559;
// F[4]=  0.228325469267790;
// F[5]=  0.197240676082559;
// F[6]=   0.122492362235122;
// F[7]=  0.048311483314256;
// F[8]= 0.017792743734167;

//Fir orden 8 freq de corte 50 hz

// F[0]=  0.018054421523238;
// F[1] = 0.048642138232548;
// F[2] =  0.122650923949400;
// F[3] = 0.196844380976817;
// F[4]=  0.227616270635995;
// F[5]=  0.196844380976817;
// F[6]=   0.122650923949400;
// F[7]= 0.048642138232548;
// F[8]=  0.018054421523238;

// b=fir1(6,(400*25e-6)*2)
Fii[0]= 0.026407724923238;
Fii[1]=  0.140531362762416;
Fii[2]=  0.333060912314346;
Fii[3]=  0.333060912314346;
Fii[4]=   0.140531362762416;
Fii[5]= 0.026407724923238;
Fii[6]= 0.022257279513162;


//FIltros corriente Ixy
//ellip orden 3 [b,a]=ellip(3,.1,30,(250*20e-6)*2);

B1i= 0.021690785231586;
B2i=-0.014104901191657;
B3i=B2i;
B4i= B1i;

A2i=-2.552077052600769;
A3i= 2.226544566468828;
A4i=-0.659295745788200;

// B1z = 0.002199178731552;
// B2z = -0.008779073382963;
// B3z = 0.010960623495097;
// B4z= 0.0;
// B5z= -0.010960623495097;
// B6z= 0.008779073382963;
// B7z= -0.002199178731552;

// A2z =  -5.954351818165614;
// A3z = 14.780217369405085;
// A4z=-19.577115962071396;
// A5z= 14.593582195531454;
// A6z=-5.804926629691336;
// A7z= 0.962594858984539;
///FIltros promedio 
//[Bp,Ap]=ellip(3,.1,30,(10*Ts*Tv)*2);
Bp[0]=   0.003686868316212;
Bp[1]= -0.003683370618120;
Bp[2]=Bp[1];
Bp[3]= Bp[0];

Ap[0]=-2.970721605692845;
Ap[1]=  2.942031088461970;
Ap[2]=-0.971302487372940;
//Ts = 20 us
// Bp[0]=  0.015379170816530;
// Bp[1]=  -0.015102339737953;
// Bp[2]=Bp[1];
// Bp[3]= Bp[0];

// Ap[0]=-2.869926787561233;
// Ap[1]= 2.751152720264573;
// Ap[2]=-0.880672270546186;



Vx = Vdc*0.5;
Exy_ref = Vc_ref*Vc_ref*C*0.5;
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
//ki_s = 0.0;
Umax_i=Vc_ref*nc;
Umax_v=250.0*Io_max;
Umin_i=-Umax_i;
Umin_v=-Umax_v;
}

void ellip(Float64 *B,Float64 *A,Float64 *In,Float64 *Out,int Order)
{
	int i;
	Float64 Temp = In[0]*B[0];
	for(i=1;i<=Order;i++)
		Temp = Temp+In[i]*B[i]-Out[i]*A[i];
	Out[0]=Temp;
	for(i=Order;i>=1;i--){
		In[i]=In[i-1];
		Out[i]=Out[i-1];
	}
}



void FIR(Float64 *B,Float64 *In,Float64 *Out,int Order)
{
	int i;
	Float64 Temp = 0.0;
	for(i=0;i<=Order;i++)
		Temp = Temp+In[i]*B[i];
	*Out=Temp;
	for(i=Order;i>=1;i--)
		In[i]=In[i-1];
	
}



void PI_tustin(Float64 kp,Float64 ki,Float64 h,Float64 *in_1,Float64 *out_1, Float64 *in,Float64 *out,Float64 umax,Float64 umin)
{
	Float64 out_temp = 0.0;
	out_temp = (ki*h)/2.0*(*in+*in_1)+*out_1;
	if(out_temp-*in*kp>umax)out_temp=umax-*in*kp;
	if(out_temp-*in*kp<umin)out_temp=umin-*in*kp;
	*out_1=out_temp;
	*out = out_temp+*in*kp;
	*in_1=*in;
	
}

void Sinc_ang(Float64 *Vabc,Float64 *Valp,Float64 *Vbeta,Float64 *ang)
{
	*Valp = Vabc[0]-0.5*Vabc[1]-0.5*Vabc[2];
	*Vbeta = sqrt(3.0)*0.5*(Vabc[1]-Vabc[2]);
	*ang = atan2(*Vbeta,*Valp);
}

void PI_tustin_3ph(Float64 kp,Float64 ki,Float64 h,Float64 *in_1,Float64 *out_1, Float64 *in,Float64 *out,Float64 umax,Float64 umin)
{
int i;
Float64 out_temp[3] = {0.0,0.0,0.0};
	for(i=0;i<3;i++){
		out_temp[i] = (ki*h)/2.0*(in[i]+in_1[i])+out_1[i];
		if(out_temp[i]-in[i]*kp>umax)
			out_temp[i]=umax-in[i]*kp;
		if(out_temp[i]-in[i]*kp<umin)
			out_temp[i]=umin-in[i]*kp;
		out_1[i]=out_temp[i];
		in_1[i]=in[i];
		out[i]= out_temp[i]+kp*in[i];
	}
}

void PI_tustin_dq(Float64 kp,Float64 ki,Float64 h,Float64 *in_1,Float64 *out_1, Float64 *in,Float64 *out,Float64 umax,Float64 umin)
{
int i;
Float64 out_temp[2] = {0.0,0.0};
	for(i=0;i<2;i++){
		out_temp[i] = (ki*h)/2.0*(in[i]+in_1[i])+out_1[i];
		if(out_temp[i]-in[i]*kp>umax)
			out_temp[i]=umax-in[i]*kp;
		if(out_temp[i]-in[i]*kp<umin)
			out_temp[i]=umin-in[i]*kp;
		out_1[i]=out_temp[i];
		in_1[i]=in[i];
		out[i]= out_temp[i]+kp*in[i];
	}
}

void PI_tustin_dq_AW(Float64 kp,Float64 ki,Float64 h,Float64 *x_1,Float64 *x,Float64 *out_1, Float64 *in,Float64 *out,Float64 umax,Float64 umin)
{
	int i;
	Float64 out_temp[2] = {0.0,0.0};
	Float64 a1=kp-2.0*ki*h;
	Float64 a2=ki*h-kp;
	Float64 b1=ki*h*ki*h-ki*h*kp;
	Float64 b2=kp*ki*h;	
	
	for(i=0;i<2;i++){
	out_temp[i] =(in[i]-x[i])*ki*h;
	if(out_temp[i]>umax)out_temp[i]=umax;

	if(out_temp[i]<umin)out_temp[i]=umin;
	x[i]=(out_1[i]*a1+out_temp[i]*a2-x_1[i]*b1)/b2;
	x_1[i]=x[i];
	out[i] = out_temp[i];
	out_1[i]=out_temp[i];
	}
}


void PI_forward_AW(Float64 kp,Float64 ki,Float64 h,Float64 *x_1,Float64 *x,Float64 *out_1, Float64 *in,Float64 *out,Float64 umax,Float64 umin)
{
	Float64 out_temp = 0.0;
	Float64 a1=kp-2.0*ki*h;
	Float64 a2=ki*h-kp;
	Float64 b1=ki*h*ki*h-ki*h*kp;
	Float64 b2=kp*ki*h;

	out_temp =(*in-*x)*ki*h;
	if(out_temp>umax)out_temp=umax;

	if(out_temp<umin)out_temp=umin;
	*x=*out_1*a1/b2+out_temp*a2/b2-*x_1*b1/b2;
	*x_1=*x;
	*out = out_temp;
	*out_1=out_temp;
}
void PI_forward_AW_3ph(Float64 kp,Float64 ki,Float64 h,Float64 *x_1,Float64 *x,Float64 *out_1, Float64 *in,Float64 *out,Float64 umax,Float64 umin)
{
	int i;
	Float64 out_temp[3] = {0.0,0.0,0.0};
	Float64 a1=kp-2.0*ki*h;
	Float64 a2=ki*h-kp;
	Float64 b1=ki*h*ki*h-ki*h*kp;
	Float64 b2=kp*ki*h;	
	
	for(i=0;i<3;i++){
	out_temp[i] =(in[i]-x[i])*ki*h;
	if(out_temp[i]>umax)out_temp[i]=umax;

	if(out_temp[i]<umin)out_temp[i]=umin;
	x[i]=out_1[i]*a1/b2+out_temp[i]*a2/b2-x_1[i]*b1/b2;
	x_1[i]=x[i];
	out[i] = out_temp[i];
	out_1[i]=out_temp[i];
	}
}


void PI_tustin_AW_3ph(Float64 kp,Float64 ki,Float64 h,Float64 *x_1,Float64 *x,Float64 *out_1, Float64 *in,Float64 *out,Float64 umax,Float64 umin)
{
	int i;
	Float64 out_temp[3] = {0.0,0.0,0.0};
	Float64 a1=2.0*kp-3.0*ki*h;
	Float64 a2=ki*h-2.0*kp;
	Float64 b1=ki*h*ki*h-2.0*ki*h*kp;
	Float64 b2=ki*h*ki*h+2.0*kp*ki*h;
	
	for(i=0;i<3;i++){
	out_temp[i] =(in[i]-x[i])*ki*h;
	if(out_temp[i]>umax)out_temp[i]=umax;

	if(out_temp[i]<umin)out_temp[i]=umin;
	x[i]=(out_1[i]*a1+out_temp[i]*a2-x_1[i]*b1)/b2;
	x_1[i]=x[i];
	out[i] = out_temp[i];
	out_1[i]=out_temp[i];
	}
}
void PI_tustin_AW(Float64 kp,Float64 ki,Float64 h,Float64 *x_1,Float64 *x,Float64 *out_1, Float64 *in,Float64 *out,Float64 umax,Float64 umin)
{
	Float64 out_temp = 0.0;
	Float64 a1=2.0*kp-3.0*ki*h;
	Float64 a2=ki*h-2.0*kp;
	Float64 b1=ki*h*ki*h-2.0*ki*h*kp;
	Float64 b2=ki*h*ki*h+2.0*kp*ki*h;

	out_temp =((*in)-(*x))*ki*h;
	if(out_temp>umax)out_temp=umax;

	if(out_temp<umin)out_temp=umin;
	*x=((*out_1)*a1+out_temp*a2-(*x_1)*b1)/b2;
	*x_1=*x;
	*out = out_temp;
	*out_1=out_temp;
}

void PI_back(Float64 kp,Float64 ki,Float64 h,Float64 *in_1,Float64 *out_1, Float64 *in,Float64 *out,Float64 umax,Float64 umin)
{
	Float64 out_temp = 0.0;
	out_temp = (ki*h)/2.0*(*in+*in_1)+*in*kp+*out_1;
	if(out_temp>umax)out_temp=umax;

	if(out_temp<umin)out_temp=umin;
	
	*out = out_temp;
	*in_1=*in;
	*out_1=out_temp;
}
void PI_back_dq(Float64 kp,Float64 ki,Float64 h,Float64 *in_1,Float64 *out_1, Float64 *in,Float64 *out,Float64 umax,Float64 umin)
{
	Float64 out_temp[2] = {0.0,0.0};
	int i =0;
	for(i=0;i<2;i++){
	out_temp[i] = (ki*h)/2.0*(in[i]+in_1[i])+in[i]*kp+out_1[i];
	if(out_temp[i]>umax)out_temp[i]=umax;

	if(out_temp[i]<umin)out_temp[i]=umin;
	
	out[i] = out_temp[i];
	in_1[i]=in[i];
	out_1[i]=out_temp[i];
	}
}

void abc2dq_pos(double *abc, double *dq,double sinc)
{
//	double alfa = (float) 2/3*(abc[1]-1/2*abc[2]-1/2*abc[3]);
//	double beta= (float) 2/3*(sqrt(3)/2*abc[2]-sqrt(3)/2*abc[3]);
//printf("%f %f %f\n",abc[0],abc[1],abc[2]);
	dq[0]=0.66666667*(abc[0]-0.5*abc[1]-0.5*abc[2])*cos(sinc)+ 0.66666667*(sqrt(3.0)*0.5*abc[1]-sqrt(3.0)*0.5*abc[2])*sin(sinc);
	dq[1]=-0.66666667*(abc[0]-0.5*abc[1]-0.5*abc[2])*sin(sinc)+(0.66666667*(sqrt(3.0)*0.5*abc[1]-sqrt(3.0)*0.5*abc[2]))*cos(sinc);
	
}

void abc2dq_neg(double *abc, double *dq,double sinc)
{
	//double alfa = (float) 2/3*(abc[1]-1/2*abc[2]-1/2*abc[3]);
	//double beta= (float) 2/3*(-sqrt(3)/2*abc[2]+sqrt(3)/2*abc[3]);
	dq[0]=(double)2/3*(abc[0]-0.5*abc[1]-0.5*abc[2])*cos(sinc)+ (double)2/3*(-sqrt(3.0)*0.5*abc[1]+sqrt(3.0)*0.5*abc[2])*sin(sinc);
	dq[1]=(double)-2/3*(abc[0]-0.5*abc[1]-0.5*abc[2])*sin(sinc)+((double) 2/3*(-sqrt(3.0)*0.5*abc[1]+sqrt(3.0)*0.5*abc[2]))*cos(sinc);
}
void dq2abc_pos(double *dq,double *abc,double sinc)
{
	abc[0]=dq[0]*cos(sinc)-dq[1]*sin(sinc);
	abc[1]=-0.5*(dq[0]*cos(sinc)-dq[1]*sin(sinc))+sqrt(3.0)*0.5*(dq[0]*sin(sinc)+dq[1]*cos(sinc));
	abc[2]=-0.5*(dq[0]*cos(sinc)-dq[1]*sin(sinc))-sqrt(3.0)*0.5*(dq[0]*sin(sinc)+dq[1]*cos(sinc));	
}

void dq2Vm(double *dq,double *m,double sinc)
{
	*m=dq[0]*cos(sinc)-dq[1]*sin(sinc);
}

void dq2abc_neg(double *dq,double *abc,double sinc)
{
	abc[0]=dq[0]*cos(sinc)-dq[1]*sin(sinc);
	abc[1]=-0.5*(dq[0]*cos(sinc)-dq[1]*sin(sinc))-sqrt(3.0)/2.0*(dq[0]*sin(sinc)+dq[1]*cos(sinc));
	abc[2]=-0.5*(dq[0]*cos(sinc)-dq[1]*sin(sinc))+sqrt(3.0)/2.0*(dq[0]*sin(sinc)+dq[1]*cos(sinc));	
}
//Usar en una sola variable

void s_xy2Dec(Float64 *xy,Float64 *o,Float64 *s,Float64 *z, Float64 *m)
{
	o[0]=0.166666667*(2.0*xy[0]+2.0*xy[3]-xy[1]-xy[4]-xy[2]-xy[5]);//1/6(2 ipa+2ina-ipb-inb-ipc-inc)
	o[1]=0.166666667*(-xy[0]-xy[3]+2.0*xy[1]+2.0*xy[4]-xy[2]-xy[5]);
	o[2]=0.166666667*(-xy[0]-xy[3]-xy[1]-xy[4]+2.0*xy[2]+2.0*xy[5]);
	
	*s=0.166666667*(xy[0]-xy[3]+xy[1]-xy[4]+xy[2]-xy[5]);
	
	z[0]=0.166666667*(2.0*xy[0]-2.0*xy[3]-xy[1]+xy[4]-xy[2]+xy[5]);
	z[1]=0.166666667*(-xy[0]+xy[3]+2.0*xy[1]-2.0*xy[4]-xy[2]+xy[5]);
	z[2]=0.166666667*(-xy[0]+xy[3]-xy[1]+xy[4]+2.0*xy[2]-2.0*xy[5]);
	
	*m=0.166666667*(xy[0]+xy[1]+xy[2]+xy[3]+xy[4]+xy[5]);
}
void s_bal(Float64 *o,Float64 *z)
{
	Float64 otemp[3]={0.0,0.0,0.0};
	Float64 ztemp[3]={0.0,0.0,0.0};
	
	otemp[0]=(2.0*o[0]-o[1]-o[2])/3.0;
	otemp[1]=(-o[0]+2.0*o[1]-o[2])/3.0;
	otemp[2]=(-o[0]-o[1]+2.0*o[2])/3.0;
		
	ztemp[0]=(2.0*z[0]-z[1]-z[2])/3.0;
	ztemp[1]=(-z[0]+2.0*z[1]-z[2])/3.0;
	ztemp[2]=(-z[0]-z[1]+2.0*z[2])/3.0;
	
	o[0]=otemp[0];
	o[1]=otemp[1];
	o[2]=otemp[2];
	
	z[0]=ztemp[0];
	z[1]=ztemp[1];
	z[2]=ztemp[2];

}
void s_Dec2xy(Float64 *xy,Float64 *o,Float64 s,Float64 *z, Float64 m)
{
	xy[0]=o[0]+s+z[0]+m; //pa
	xy[1]=o[1]+s+z[1]+m; //pb 
	xy[2]=o[2]+s+z[2]+m; //pc 
	xy[3]=o[0]-s-z[0]+m; //na 
	xy[4]=o[1]-s-z[1]+m; //nb
	xy[5]=o[2]-s-z[2]+m; //nc
	
}


