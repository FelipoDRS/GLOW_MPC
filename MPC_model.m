function [dydt,W,P,rho,alpha]= MPC_model(t,y,par,u,Dist)
%Felipo Soares, model from Control structure design for stabilizing
%unstable gas-lift oil wells
%Esmaeil Jahanshahi, Sigurd Skogestad and Henrik Hansen,2012
%x1=0.95
%x2+x3 > 1790
%x3<20812 & x3>1789
%Ptt=22e5, rho_gt=12 rho_gt=y(2)/(27.38-y(3)/760) %y(2)=296
Pgs=Dist(1);
GOR=abs(par(1));
Pres=par(2)*1e6;
%wres=18;
PI=par(3)*1e-6;
Kinj=par(4)*1e-4;
u1=u(1);
u2=u(2);

R=8314;
g=9.81;
rho_l=760;
Mg=16.7;
Ta=348;
Va=64.34;
La=2048;
Pgs=140*1.013e5;
Vt=25.03;
Sbh=0.0314;
Lbh=75;
Tt=369.4;
%GOR=0;
%Pres=160*1.013e5;
Dt=0.134;
Lt=2048;
%PI=2.47e-6;
Kgs=9.98e-05;%altered
%Kinj=1.40e-4;
Kpr=2.90e-3;
P0=1.013e5;
u1=u(1);
u2=u(2);
%Pgs=u(3);
%y(2)=max(0,y(2));
dydt=zeros(3,1);
Pat=R*Ta*y(1)/(Mg*Va);
Pab=Pat+y(1)*g*La/Va;
rho_gab=Pab*Mg/(R*Ta);%convert bar to Pa
rho_gin=Pgs*Mg/(R*Ta);%convert bar to Pa
w_gin=Kgs*u2*sqrt(rho_gin*max([Pgs-Pat,0]));

rho_gt=y(2)/(Vt+Sbh*Lbh - y(3)/rho_l);
Ptt = rho_gt*R*Tt/Mg;
rho_mix=(y(2)+y(3)-rho_l*Sbh*Lbh)/Vt;
alpha_l=max([0,(y(3)-rho_l*Sbh*Lbh)/(Vt*rho_l)]);
alpha_mgb=GOR/(GOR+1);

Ptb = Ptt+rho_mix*Lt*g;
w_ginj=Kinj*sqrt(rho_gab*max([Pab-Ptb,0]));
Pbh = Ptb+rho_l*g*Lbh;

w_res=PI*max([Pres-Pbh,0]);
w_lres=(1-alpha_mgb)*w_res;
w_gres=alpha_mgb*w_res;

rho_gtb=Ptb*Mg/(R*Tt);
alpha_lb=(w_lres*rho_gtb)/(w_lres*rho_gtb+(w_ginj+w_gres)*rho_l+1e-20);
alpha_lt=min([max([2*alpha_l-alpha_lb,0]),1]);

rho_mixt=alpha_lt*rho_l+(1-alpha_lt)*rho_gt;
w_out=Kpr*u1*sqrt(rho_mixt*max([Ptt-P0,0]));
Qout=w_out/rho_mixt;
alpha_mgt=(1-alpha_lt)*rho_gt/(alpha_lt*rho_l+(1-alpha_lt)*rho_gt);

w_gout= alpha_mgt*w_out;
w_lout=(1-alpha_mgt)*w_out;

dydt(1)=w_gin-w_ginj;
dydt(2)=w_ginj+w_gres- w_gout;
dydt(3)=w_lres-w_lout;

W=[w_gin,w_ginj,w_res,w_lres,w_gres,w_out,w_gout,w_lout];
P=[Pat, Pab, Ptt,Ptb,Pbh];
rho=[rho_gab,rho_gin,rho_gt,rho_mix,rho_gtb,rho_mixt];
alpha=[alpha_l,alpha_mgb,alpha_lb,alpha_lt,alpha_mgt];


end