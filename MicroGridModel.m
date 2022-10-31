close all
clear all
clc

%Induction Machine 1 Parameters**********************************************
rs=0.935;Xls=0.754;Xm=26.13;Xlr=0.754;rr=0.816;
Xss=Xls+Xm;Xrr=Xlr+Xm;
omegab=314;omegae=314;
H=0.705;S0=0.05;
idsoe=5;idroe=5;iqroe=5;iqsoe=5;

E=[Xss 0 Xm 0 0;0 Xss 0 Xm 0;Xm 0 Xrr 0 0;0 Xm 0 Xrr 0;0 0 0 0 -2*H*omegab];
E=(1/omegab).*E;
F=[rs (omegae/omegab)*Xss 0 (omegae/omegab)*Xm 0;-(omegae/omegab)*Xss rs -(omegae/omegab)*Xm 0 0;0 S0*(omegae/omegab)*Xm rr S0*(omegae/omegab)*Xrr -Xm*idsoe-Xrr*idroe;-S0*(omegae/omegab)*Xm 0 -S0*(omegae/omegab)*Xrr rr Xm*iqsoe+Xrr*iqroe;Xm*idroe -Xm*iqroe Xm*idsoe Xm*iqsoe 0];
F=(-1).*F;
A1=inv(E)*F;
AindD1=A1;
BindD1=inv(E);
CindD1=[1 0 0 0 0];
DindD1=[0];
sys1=ss(AindD1,BindD1,CindD1,DindD1);
 AInd1=eig(sys1);
%Induction Machine 2 Parameters********************************************
rs=0.935;Xls=0.754;Xm=26.13;Xlr=0.754;rr=0.816;
Xss=Xls+Xm;Xrr=Xlr+Xm;
omegab=377;omegae=377;
% HStep=HStep+0.2;
H=0.49;S0=0.05;
idsoe=7.5;idroe=7.5;iqroe=7.5;iqsoe=7.5;

E=[Xss 0 Xm 0 0;0 Xss 0 Xm 0;Xm 0 Xrr 0 0;0 Xm 0 Xrr 0;0 0 0 0 -2*H*omegab];
E=(1/omegab).*E;
F=[rs (omegae/omegab)*Xss 0 (omegae/omegab)*Xm 0;-(omegae/omegab)*Xss rs -(omegae/omegab)*Xm 0 0;0 S0*(omegae/omegab)*Xm rr S0*(omegae/omegab)*Xrr -Xm*idsoe-Xrr*idroe;-S0*(omegae/omegab)*Xm 0 -S0*(omegae/omegab)*Xrr rr Xm*iqsoe+Xrr*iqroe;Xm*idroe -Xm*iqroe Xm*idsoe Xm*iqsoe 0];
F=(-1).*F;
A2=inv(E)*F;
AindD2=A2;
BindD2=inv(E);
CindD2=[1 0 0 0 0];
DindD2=[0];
sys2=ss(AindD2,BindD2,CindD2,DindD2);
AInd2=eig(sys2);

%**************************************************************************
%Line 1
R1=0.06;X1=0.314;
%Line 2
R2=0.04;X2=0.471;

%Intial Parameter DG1**********************************************************
r1=0.0052;
Xd1=2.86;Xlkd1=0.0208;Xlfd1=0.6157;Xkd1=2.68;Xfd1=3.2757;rkd1=0.015381;
rfd1=0.0026;H=2.9;Xls=0.2;Xq1=2;Xlk1q1=0.0284;Xlk2q1=0.001;Xk1q1=1.8284;
Xk2q1=1.801;rk1q1=0.0057;rk2q1=0.0015;Da=0;
Xmq1=Xq1-Xls;
Xmd1=Xd1-Xls;

Sb=5;%S base kVA
Vb=13.8; %V base v

%Initial Condition 
omegab=314;
omegae=314;
delta1=0.5; % 30 darje
id10=0.9;
iq10=0.7;
ifd10=1;
Vd10=1;
Vq10=0.9;
idfd10=1;
ifd1=1;

%DG1 Model*****************************************************************
T1=[cos(delta1) -sin(delta1);sin(delta1) cos(delta1)];
Wp=(1/omegab).*([-1.*Xq1 0;0 -1.*Xd1]);
Yp=(1/omegab).*([Xmq1 Xmq1 0 0 0 Xq1.*id10;0 0 Xmd1 Xmd1 0 -1.*Xd1.*iq10]);
Qp=(1/omegab).*([-1.*Xmq1 -1.*Xmq1 0 0 0 0;0 0 -1.*(Xmd1.^2/rfd1) -1.*Xmd1 0 0]);
Qp=Qp';
Sp=(1/omegab).*([Xk1q1 Xmq1 0 0 0 Xmq1.*id10;Xmq1 Xk2q1 0 0 0 Xmq1.*id10;0 0 (Xmd1.*Xfd1/rfd1) (Xmd1.^2/rfd1) 0 -1.*(Xmd1.^2.*iq10/rfd1);0 0 Xmd1 Xkd1 0 -1.*(Xmd1.*iq10);0 0 0 0 2.*H.*omegab Da.*omegab;0 0 0 0 0 omegab]);

Wk=[-1.*r1 -1.*(omegae.*Xd1/omegab);(omegae.*Xq1/omegab) -1.*r1];
Yk=[0 0 (omegae.*Xmd1/omegab) (omegae.*Xmd1/omegab) -1.*Xd1.*id10+Xmd1.*ifd10 r1.*id10-(omegae.*Xd1.*iq10/omegab)+Vd10;-1.*(omegae.*Xmq1/omegab) -1.*(omegae.*Xmq1/omegab) 0 0 Xq1.*iq10 r1.*iq10-(omegae.*Xq1.*id10/omegab)+Vq10];
Qk=[0 0 0 0 Xmq1.*id10-Xmd1.*(id10-idfd10) 0;0 0 0 0 -1.*Xmd1.*iq10+Xmq1.*iq10 0];
Qk=Qk';
Xtit=-1.*id10.*(Xmq1 -Xmd1.*(id10-ifd1))-iq10.*(Xmd1.*iq10-Xmq1.*iq10);
Sk=[rk1q1 0 0 0 0 0;0 rk2q1 0 0 0 0;0 0 Xmd1 0 0 0;0 0 0 rkd1 0 0; -1.*Xmq1.*id10 Xmq1.*id10 Xmd1.*iq10 Xmd1.*iq10 0 Xtit;0 0 0 0 -1.*omegab 0];

E=[(inv(T1))*Wp*T1 (inv(T1))*Yp;Qp*T1 Sp];
F=-1.*([(inv(T1))*Wk*T1 (inv(T1))*Yk;Qk*T1 Sk]);
ADG1=(inv(E))*F;
BVDG1=(inv(E))*[1 1 0 0 0 0 0 0;1 1 0 0 0 0 0 0]';
BUDG1=(inv(E))*[0 0 1 1 1 1 1 1;0 0 1 1 1 1 1 1]';
ADGTot=eig(ADG1)

%Sakht Matrix Bv DG1
FDG1=[X1/omegab 0;0 X1/omegab];
ADotDG1=[R1 omegae*X1/omegab 0 0 0 0 0 0;-omegae*X1/omegab R1 0 0 0 0 0 0];%iq1 & id1
AOwDG1=[1 0;0 1];%Vqpcc & Vdpcc
%*****Modified Bv DG1
FDG1M=BVDG1*FDG1;
ADotDG1M=BVDG1*ADotDG1;
AOwDG1M=BVDG1*AOwDG1;

%**********************PV*************************************************
%Inverter Model

%Parameter*****************************************************************
fs=8000;%kHZ
Lf=1.75;%mH
Cf=50;%micro Farad
rLf=0.1;%ohm
Lc=0.35;%mH
rLc=0.03;%ohm
omegac=0.8;
mp=15e-5;
nq=2.3e-3;
Kpv=0.05;
Kiv=400;
Kpc=10.5;
Kic=15.5e3;
F=0.75;
omegan=300;
%**************************************************************************

%Intial Condition**********************************************************
Vod=[380.3 381.8 380.4];
Iod=[11.4 11.4 11.4];
Ild=[11.4 11.4 11.4];
VbD=[379.5 380.5 379];
omega0=314;
Voq=[0 0 0];
Ioq=[0.4 -1.45 1.25];
Ilq=[-5.5 7.3 4.6];
VbQ=[-6 -6 -5];
delta0=[0 1.9e-3 -0.0113];
%**************************************************************************

%Matrices******************************************************************
AP=[0 -mp 0;0 omegac 0;0 0 -omegac];
CPv=[0 0 -nq;0 0 0];
Bv1=[1 0;0 1];
BC1=[1 0;0 1];
Dv1=[Kpv 0;0 Kpv];
BLCL1=[1/Lf 0;0 1/Lf;0 0;0 0;0 0;0 0];
BLCL2=[0 0;0 0;0 0;0 0;-1/Lc 0;0 -1/Lc];
BLCL3=[Ilq(1,1)-Ild(1,1);Ilq(1,2)-Ild(1,2);Voq(1,1)-Vod(1,1);Voq(1,2)-Vod(1,2);Ioq(1,1)-Iod(1,1);Ioq(1,2)-Iod(1,2)];
DC1=[Kpc 0;0 Kpc];
TVINV=[-VbD(1,2).*sin(delta0(1,2))+VbQ(1,2).*cos(delta0(1,2));-VbD(1,2).*cos(delta0(1,2))-VbQ(1,2).*sin(delta0(1,2))];
CPw=[0 -mp 0];
CV=[Kiv 0;0 Kiv];
CC=[Kic 0;0 Kic];
Bp=[0 0 0 0 0 0;0 0 omegac.*Iod(1,2) omegac.*Ioq(1,2) omegac.*Vod(1,2) omegac.*Voq(1,2);0 0 omegac.*Ioq(1,2) -omegac.*Iod(1,2) -omegac.*Voq(1,2) -omegac.*Vod(1,2)];
Bv2=[0 0 -1 0 0 0;0 0 0 -1 0 0];
Dv2=[0 0 -Kpv -omegan*Cf F 0;0 0 omegan*Cf -Kpv 0 F];
BC2=[-1 0 0 0 0 0;0 -1 0 0 0 0];
DC2=[-Kpc -omegan*Lf 0 0 0 0;omegan*Lf -Kpc 0 0 0 0];
ALCL=[-rLf/Lf omega0 -1/Lf 0 0 0;-omega0 -rLf/Lf 0 -1/Lf 0 0;1/Cf 0 0 omega0 -1/Cf 0;0 1/Cf -omega0 0 0 -1/Cf;0 0 1/Lc 0 -rLf/Lc omega0;0 0 0 1/Lc -omega0 -rLf/Lc];
TS=[cos(delta0(1,2)) -sin(delta0(1,2));sin(delta0(1,2)) cos(delta0(1,2))];
TC=[-Iod(1,2).*sin(delta0(1,2))-Ioq(1,2).*cos(delta0(1,2));Iod(1,2).*cos(delta0(1,2))-Ioq(1,2).*sin(delta0(1,2))];
%**************************************************************************

%INVETER STATE SPACE MODEL*************************************************
AINV1=zeros(13,13);
AINV1(1:3,1:3)=AP;
AINV1(4:5,1:3)=Bv1*CPv;
AINV1(6:7,1:3)=BC1*Dv1*CPv;
Z1=BLCL1*DC1;Z2=Z1*Dv1;Z3=Z2*CPv;
AINV1(8:13,1:3)=Z3+BLCL2*[TVINV(1,1) 0 0;TVINV(2,1) 0 0]+BLCL3*CPw;

AINV1(6:7,4:5)=BC1*CV;
AINV1(8:13,4:5)=BLCL1*DC1*CV;

AINV1(8:13,6:7)=BLCL1*CC;

AINV1(1:3,8:13)=Bp;
AINV1(4:5,8:13)=Bv2;
AINV1(6:7,8:13)=BC1*Dv2+BC2;
AINV1(8:13,8:13)=ALCL+BLCL1*(DC1*Dv2+DC2);

% A=AINV;

BINV1=zeros(13,3);
BINV1(8:13,1:2)=BLCL2*inv(TS);
BINV1(1,3)=-1;
AINVPV=AINV1;
% BUINV1=BINV;

CINV1=zeros(3,13);
CINV1(1,1:3)=CPw;
CINV1(2:3,1)=TC;
CINV1(2:3,12:13)=TS;

% C1=CINV1;

DINV1=[0];
APVTot=eig(AINVPV);


%Sakht Matrix Bv DG1
FPV=[X2/omegab 0;0 X2/omegab];
ADotPV=[R2 omegae*X2 0 0 0 0 0 0 0 0 0 0 0/omegab;-omegae*X2/omegab R2 0 0 0 0 0 0 0 0 0 0 0];%iq2 & id2
AOwPV=[1 0;0 1];%Vqpcc & Vdpcc
%*****Modified Bv DG1
FPVM=FPV;
ADotPVM=ADotPV;
AOwPVM=AOwPV;


%**********************Network Parameter***********************************

%Line 1
R1=0.06;X1=0.314;
%Line 2
R2=0.04;X2=0.471;
Rsl=0.59;Xsl=3.15;
XC=1.803;
omegab=314;omegae=314;
% omegae=omegab;
%Network Small Signal Model
ANET=[-omegab*Rsl/Xsl -omegae omegab/Xsl 0;omegae -omegab*Rsl/Xsl 0 omegab/Xsl;-omegab*XC 0 0 -omegae;0 -omegab*XC omegae 0];
ANETTot=eig(ANET);

ADotNet=[omegab*XC omegab*XC -omegab*XC -omegab*XC -omegab*XC -omegab*XC;omegab*XC omegab*XC -omegab*XC -omegab*XC -omegab*XC -omegab*XC;];
%**********************Thermal Load***********************************

rT=1;%Thermal Load Resistance
syms x1 x2 x3 x4
omegae=314;omegab=314;r=1;T=20;Pl1=0.01*x1;Pl2=0.01*x2;
t1=-omegae*x4;
t2=-omegae*x3;
F1=(-r/T)*(x1^2*(x1-Pl1))/(x3^2)+(2*x1^2/x3^2)*t1;
F2=(-r/T)*(x2^2*(x2-Pl2))/(x4^2)+(2*x2^2/x4^2)*t2;



TH11=diff(F1,x1);
TH12=diff(F1,x2);
TH21=diff(F2,x1);
TH22=diff(F2,x2);

TE11=diff(F1,x3);
TE12=diff(F1,x4);
TE21=diff(F2,x3);
TE22=diff(F2,x4);

%Initial Condition of Thermal State
x10=25;
x20=25;
x30=1;
x40=1;

TH11=subs(TH11,[x1 x2 x3 x4],[x10 x20 x30 x40]);
TH12=subs(TH12,[x1 x2 x3 x4],[x10 x20 x30 x40]);
TH21=subs(TH21,[x1 x2 x3 x4],[x10 x20 x30 x40]);
TH22=subs(TH22,[x1 x2 x3 x4],[x10 x20 x30 x40]);

TE11=subs(TE11,[x1 x2 x3 x4],[x10 x20 x30 x40]);
TE12=subs(TE12,[x1 x2 x3 x4],[x10 x20 x30 x40]);
TE21=subs(TE21,[x1 x2 x3 x4],[x10 x20 x30 x40]);
TE22=subs(TE22,[x1 x2 x3 x4],[x10 x20 x30 x40]);

%Create A TOTAL State Space Model******************************************
a=zeros(37,37);
%*********************************DG1**************************************
a(1,1)=ADG1(1,1)+ADotDG1(1,1);
a(1,2)=ADG1(1,2)+ADotDG1(1,2);
a(1,3)=ADG1(1,3)+ADotDG1(1,3);
a(1,4)=ADG1(1,4)+ADotDG1(1,4);
a(1,5)=ADG1(1,5)+ADotDG1(1,5);
a(1,6)=ADG1(1,6)+ADotDG1(1,6);
a(1,7)=ADG1(1,7)+ADotDG1(1,7);
a(1,8)=ADG1(1,8)+ADotDG1(1,8);

a(2,1)=ADG1(2,1)+ADotDG1(2,1);
a(2,2)=ADG1(2,2)+ADotDG1(2,2);
a(2,3)=ADG1(2,3)+ADotDG1(2,3);
a(2,4)=ADG1(2,4)+ADotDG1(2,4);
a(2,5)=ADG1(2,5)+ADotDG1(2,5);
a(2,6)=ADG1(2,6)+ADotDG1(2,6);
a(2,7)=ADG1(2,7)+ADotDG1(2,7);
a(2,8)=ADG1(2,8)+ADotDG1(2,8);

a(3,1)=ADG1(3,1);
a(3,2)=ADG1(3,2);
a(3,3)=ADG1(3,3);
a(3,4)=ADG1(3,4);
a(3,5)=ADG1(3,5);
a(3,6)=ADG1(3,6);
a(3,7)=ADG1(3,7);
a(3,8)=ADG1(3,8);
 
a(4,1)=ADG1(4,1);
a(4,2)=ADG1(4,2);
a(4,3)=ADG1(4,3);
a(4,4)=ADG1(4,4);
a(4,5)=ADG1(4,5);
a(4,6)=ADG1(4,6);
a(4,7)=ADG1(4,7);
a(4,8)=ADG1(4,8);

a(5,1)=ADG1(5,1);
a(5,2)=ADG1(5,2);
a(5,3)=ADG1(5,3);
a(5,4)=ADG1(5,4);
a(5,5)=ADG1(5,5);
a(5,6)=ADG1(5,6);
a(5,7)=ADG1(5,7);
a(5,8)=ADG1(5,8);

a(6,1)=ADG1(6,1);
a(6,2)=ADG1(6,2);
a(6,3)=ADG1(6,3);
a(6,4)=ADG1(6,4);
a(6,5)=ADG1(6,5);
a(6,6)=ADG1(6,6);
a(6,7)=ADG1(6,7);
a(6,8)=ADG1(6,8);

a(7,1)=ADG1(7,1);
a(7,2)=ADG1(7,2);
a(7,3)=ADG1(7,3);
a(7,4)=ADG1(7,4);
a(7,5)=ADG1(7,5);
a(7,6)=ADG1(7,6);
a(7,7)=ADG1(7,7);
a(7,8)=ADG1(7,8);

a(8,1)=ADG1(8,1);
a(8,2)=ADG1(8,2);
a(8,3)=ADG1(8,3);
a(8,4)=ADG1(8,4);
a(8,5)=ADG1(8,5);
a(8,6)=ADG1(8,6);
a(8,7)=ADG1(8,7);
a(8,8)=ADG1(8,8);
%**************************************************************************


%***********************************PV*************************************
a(9,9)=AINVPV(1,1)+ADotPV(1,1);
a(9,10)=AINVPV(1,2)+ADotPV(1,2);
a(9,11)=AINVPV(1,3)+ADotPV(1,3);
a(9,12)=AINVPV(1,4)+ADotPV(1,4);
a(9,13)=AINVPV(1,5)+ADotPV(1,5);
a(9,14)=AINVPV(1,6)+ADotPV(1,6);
a(9,15)=AINVPV(1,7)+ADotPV(1,7);
a(9,16)=AINVPV(1,8)+ADotPV(1,8);
a(9,17)=AINVPV(1,9)+ADotPV(1,9);
a(9,18)=AINVPV(1,10)+ADotPV(1,10);
a(9,19)=AINVPV(1,11)+ADotPV(1,11);
a(9,20)=AINVPV(1,12)+ADotPV(1,12);
a(9,21)=AINVPV(1,13)+ADotPV(1,13);

a(10,9)=AINVPV(2,1)+ADotPV(2,1);
a(10,10)=AINVPV(2,2)+ADotPV(2,2);
a(10,11)=AINVPV(2,3)+ADotPV(2,3);
a(10,12)=AINVPV(2,4)+ADotPV(2,4);
a(10,13)=AINVPV(2,5)+ADotPV(2,5);
a(10,14)=AINVPV(2,6)+ADotPV(2,6);
a(10,15)=AINVPV(2,7)+ADotPV(2,7);
a(10,16)=AINVPV(2,8)+ADotPV(2,8);
a(10,17)=AINVPV(2,9)+ADotPV(2,9);
a(10,18)=AINVPV(2,10)+ADotPV(2,10);
a(10,19)=AINVPV(2,11)+ADotPV(2,11);
a(10,20)=AINVPV(2,12)+ADotPV(2,12);
a(10,21)=AINVPV(2,13)+ADotPV(2,13);

a(11,9)=AINVPV(3,1);
a(11,10)=AINVPV(3,2);
a(11,11)=AINVPV(3,3);
a(11,12)=AINVPV(3,4);
a(11,13)=AINVPV(3,5);
a(11,14)=AINVPV(3,6);
a(11,15)=AINVPV(3,7);
a(11,16)=AINVPV(3,8);
a(11,17)=AINVPV(3,9);
a(11,18)=AINVPV(3,10);
a(11,19)=AINVPV(3,11);
a(11,20)=AINVPV(3,12);
a(11,21)=AINVPV(3,13);

a(12,9)=AINVPV(4,1);
a(12,10)=AINVPV(4,2);
a(12,11)=AINVPV(4,3);
a(12,12)=AINVPV(4,4);
a(12,13)=AINVPV(4,5);
a(12,14)=AINVPV(4,6);
a(12,15)=AINVPV(4,7);
a(12,16)=AINVPV(4,8);
a(12,17)=AINVPV(4,9);
a(12,18)=AINVPV(4,10);
a(12,19)=AINVPV(4,11);
a(12,20)=AINVPV(4,12);
a(12,21)=AINVPV(4,13);


a(13,9)=AINVPV(5,1);
a(13,10)=AINVPV(5,2);
a(13,11)=AINVPV(5,3);
a(13,12)=AINVPV(5,4);
a(13,13)=AINVPV(5,5);
a(13,14)=AINVPV(5,6);
a(13,15)=AINVPV(5,7);
a(13,16)=AINVPV(5,8);
a(13,17)=AINVPV(5,9);
a(13,18)=AINVPV(5,10);
a(13,19)=AINVPV(5,11);
a(13,20)=AINVPV(5,12);
a(13,21)=AINVPV(5,13);

a(14,9)=AINVPV(6,1);
a(14,10)=AINVPV(6,2);
a(14,11)=AINVPV(6,3);
a(14,12)=AINVPV(6,4);
a(14,13)=AINVPV(6,5);
a(14,14)=AINVPV(6,6);
a(14,15)=AINVPV(6,7);
a(14,16)=AINVPV(6,8);
a(14,17)=AINVPV(6,9);
a(14,18)=AINVPV(6,10);
a(14,19)=AINVPV(6,11);
a(14,20)=AINVPV(6,12);
a(14,21)=AINVPV(6,13);

a(15,9)=AINVPV(7,1);
a(15,10)=AINVPV(7,2);
a(15,11)=AINVPV(7,3);
a(15,12)=AINVPV(7,4);
a(15,13)=AINVPV(7,5);
a(15,14)=AINVPV(7,6);
a(15,15)=AINVPV(7,7);
a(15,16)=AINVPV(7,8);
a(15,17)=AINVPV(7,9);
a(15,18)=AINVPV(7,10);
a(15,19)=AINVPV(7,11);
a(15,20)=AINVPV(7,12);
a(15,21)=AINVPV(7,13);

a(16,9)=AINVPV(8,1);
a(16,10)=AINVPV(8,2);
a(16,11)=AINVPV(8,3);
a(16,12)=AINVPV(8,4);
a(16,13)=AINVPV(8,5);
a(16,14)=AINVPV(8,6);
a(16,15)=AINVPV(8,7);
a(16,16)=AINVPV(8,8);
a(16,17)=AINVPV(8,9);
a(16,18)=AINVPV(8,10);
a(16,19)=AINVPV(8,11);
a(16,20)=AINVPV(8,12);
a(16,21)=AINVPV(8,13);

a(17,9)=AINVPV(9,1);
a(17,10)=AINVPV(9,2);
a(17,11)=AINVPV(9,3);
a(17,12)=AINVPV(9,4);
a(17,13)=AINVPV(9,5);
a(17,14)=AINVPV(9,6);
a(17,15)=AINVPV(9,7);
a(17,16)=AINVPV(9,8);
a(17,17)=AINVPV(9,9);
a(17,18)=AINVPV(9,10);
a(17,19)=AINVPV(9,11);
a(17,20)=AINVPV(9,12);
a(17,21)=AINVPV(9,13);

a(18,9)=AINVPV(10,1);
a(18,10)=AINVPV(10,2);
a(18,11)=AINVPV(10,3);
a(18,12)=AINVPV(10,4);
a(18,13)=AINVPV(10,5);
a(18,14)=AINVPV(10,6);
a(18,15)=AINVPV(10,7);
a(18,16)=AINVPV(10,8);
a(18,17)=AINVPV(10,9);
a(18,18)=AINVPV(10,10);
a(18,19)=AINVPV(10,11);
a(18,20)=AINVPV(10,12);
a(18,21)=AINVPV(10,13);

a(19,9)=AINVPV(11,1);
a(19,10)=AINVPV(11,2);
a(19,11)=AINVPV(11,3);
a(19,12)=AINVPV(11,4);
a(19,13)=AINVPV(11,5);
a(19,14)=AINVPV(11,6);
a(19,15)=AINVPV(11,7);
a(19,16)=AINVPV(11,8);
a(19,17)=AINVPV(11,9);
a(19,18)=AINVPV(11,10);
a(19,19)=AINVPV(11,11);
a(19,20)=AINVPV(11,12);
a(19,21)=AINVPV(11,13);


a(20,9)=AINVPV(12,1);
a(20,10)=AINVPV(12,2);
a(20,11)=AINVPV(12,3);
a(20,12)=AINVPV(12,4);
a(20,13)=AINVPV(12,5);
a(20,14)=AINVPV(12,6);
a(20,15)=AINVPV(12,7);
a(20,16)=AINVPV(12,8);
a(20,17)=AINVPV(12,9);
a(20,18)=AINVPV(12,10);
a(20,19)=AINVPV(12,11);
a(20,20)=AINVPV(12,12);
a(20,21)=AINVPV(12,13);

a(21,9)=AINVPV(13,1);
a(21,10)=AINVPV(13,2);
a(21,11)=AINVPV(13,3);
a(21,12)=AINVPV(13,4);
a(21,13)=AINVPV(13,5);
a(21,14)=AINVPV(13,6);
a(21,15)=AINVPV(13,7);
a(21,16)=AINVPV(13,8);
a(21,17)=AINVPV(13,9);
a(21,18)=AINVPV(13,10);
a(21,19)=AINVPV(13,11);
a(21,20)=AINVPV(13,12);
a(21,21)=AINVPV(13,13);

%**************************************************************************

%*********************************Static Load & Network********************
a(22,22)=ANET(1,1);
a(22,23)=ANET(1,2);
a(22,24)=ANET(1,3);
a(22,25)=ANET(1,4);

a(23,22)=ANET(2,1);
a(23,23)=ANET(2,2);
a(23,24)=ANET(2,3);
a(23,25)=ANET(2,4);

a(24,22)=ANET(3,1);
a(24,23)=ANET(3,2);
a(24,24)=ANET(3,3);
a(24,25)=ANET(3,4);

a(25,22)=ANET(4,1);
a(25,23)=ANET(4,2);
a(25,24)=ANET(4,3);
a(25,25)=ANET(4,4);

%**************************************************************************

%*********************************Induction Motor 1************************
a(26,26)=AindD1(1,1);
a(26,27)=AindD1(1,2);
a(26,28)=AindD1(1,3);
a(26,29)=AindD1(1,4);
a(26,30)=AindD1(1,5);

a(27,26)=AindD1(2,1);
a(27,27)=AindD1(2,2);
a(27,28)=AindD1(2,3);
a(27,29)=AindD1(2,4);
a(27,30)=AindD1(2,5);

a(28,26)=AindD1(3,1);
a(28,27)=AindD1(3,2);
a(28,28)=AindD1(3,3);
a(28,29)=AindD1(3,4);
a(28,30)=AindD1(3,5);

a(29,26)=AindD1(4,1);
a(29,27)=AindD1(4,2);
a(29,28)=AindD1(4,3);
a(29,29)=AindD1(4,4);
a(29,30)=AindD1(4,5);

a(30,26)=AindD1(5,1);
a(30,27)=AindD1(5,2);
a(30,28)=AindD1(5,3);
a(30,29)=AindD1(5,4);
a(30,30)=AindD1(5,5);

%**************************************************************************

%***********************************Induction Motor 2**********************
a(31,31)=AindD2(1,1);
a(31,32)=AindD2(1,2);
a(31,33)=AindD2(1,3);
a(31,34)=AindD2(1,4);
a(31,35)=AindD2(1,5);

a(32,31)=AindD2(2,1);
a(32,32)=AindD2(2,2);
a(32,33)=AindD2(2,3);
a(32,34)=AindD2(2,4);
a(32,35)=AindD2(2,5);

a(33,31)=AindD2(3,1);
a(33,32)=AindD2(3,2);
a(33,33)=AindD2(3,3);
a(33,34)=AindD2(3,4);
a(33,35)=AindD2(3,5);

a(34,31)=AindD2(4,1);
a(34,32)=AindD2(4,2);
a(34,33)=AindD2(4,3);
a(34,34)=AindD2(4,4);
a(34,35)=AindD2(4,5);

a(35,31)=AindD2(5,1);
a(35,32)=AindD2(5,2);
a(35,33)=AindD2(5,3);
a(35,34)=AindD2(5,4);
a(35,35)=AindD2(5,5);

%**************************************************************************

%****************************************Thermal Load***********************
a(36,36)=TH11;
a(36,37)=TH12;

a(37,36)=TH21;
a(37,37)=TH22;
%**************************************************************************

%Elman haie non-ghotrie A
%********************************DG1***************************************
a(1,24)=a(1,24)+AOwDG1M(1,1);
a(1,25)=a(1,25)+AOwDG1M(1,2);

a(2,24)=a(2,24)+AOwDG1M(2,1);
a(2,25)=a(2,25)+AOwDG1M(2,2);

%********************************PV***************************************
a(9,24)=a(9,24)+AOwPV(1,1);
a(9,25)=a(9,25)+AOwPV(1,2);

a(10,24)=a(10,24)+AOwPV(2,1);
a(10,25)=a(10,25)+AOwPV(2,2);

%******************************Network & Static Load***********************
a(24,1)=a(24,1)+ADotNet(1,1);
a(25,2)=a(25,2)+ADotNet(2,1);

a(24,9)=a(24,9)+ADotNet(1,2);
a(25,10)=a(25,10)+ADotNet(2,2);

a(24,23)=a(24,23)+ADotNet(1,3);
a(25,24)=a(25,24)+ADotNet(2,3);

a(24,26)=a(24,26)+ADotNet(1,4);
a(25,27)=a(25,27)+ADotNet(2,4);

a(24,31)=a(24,31)+ADotNet(1,5);
a(25,32)=a(25,32)+ADotNet(2,5);

a(24,36)=a(24,36)+ADotNet(1,6);
a(25,37)=a(25,37)+ADotNet(2,6);

%******************************Therml Load*********************************
a(36,24)=TE11;
a(36,25)=TE12;

a(37,24)=TE21;
a(37,25)=TE22;

%**************************************************************************


%**************** F Matrix Generation
I=eye(37,37);
f=zeros(37,37);
%***************** DG1
f(1,1)=FDG1(1,1);
f(1,2)=FDG1(1,2);

f(2,1)=FDG1(2,1);
f(2,2)=FDG1(2,2);
%**************************************************************************
anew=inv(I-f)*a;
EigenValue=eig(anew)
REAL_z=real(EigenValue);
IMAG_z=imag(EigenValue);

plot(real(EigenValue),imag(EigenValue),'+')
xlabel('Real(1/sec)')
ylabel('Im (rad/sec)')

