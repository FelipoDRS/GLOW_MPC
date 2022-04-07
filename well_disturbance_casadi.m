
%Felipo Soares, model based on Control structure design for stabilizing
%unstable gas-lift oil wells 
%Esmaeil Jahanshahi, Sigurd Skogestad and Henrik Hansen, 2012
%Expanded to include 2 things, Peng-Robinson state equation
%and replacing wres with the past values of w_res

%If you use this code please cite Development of a Nonlinear Model 
%Predictive Control for Stabilization of a Gas-Lift Oil Well
%Felipo Doval Rojas Soares, Argimiro Resende Secchi 
%and Maurício Bezerra de Souza Jr.*, 2022
%https://pubs.acs.org/doi/10.1021/acs.iecr.1c04728

%Recommended initial state=[3000,80,5000,10]
%Recommended disturbances=[140e5,0,160e5,2.47e-6]

% Limits to ensure the functioning of the model:

%x2+x3 > 1790
%x3<20812 & x3>1789
import casadi.*
y = SX.sym('y',4);
Dist=SX.sym('pgs',4);
u=SX.sym('u',2);

%% Peng Robinson EoS calculation
V_PR=SX.sym('V_PR');
VI_PR=SX.sym('VI_PR');
VI_PR_guess=SX.sym('VI_PR',10);
T_PR=SX.sym('T_PR');
P_PR=SX.sym('P_PR');

Tc_PR=190.6;
Pc_PR=46.04*1.013e5;
w_PR=0.011;

R_PR=8.314;
a_PR=0.45724*R_PR^2*Tc_PR^2/Pc_PR;
b_PR=0.07780*R_PR*Tc_PR/Pc_PR;
         
k_PR=0.3764+1.54226*w_PR-0.26992*w_PR^2;

alpha_PR=(1+k_PR*(1-sqrt(T_PR./Tc_PR))).^2;
p_PR=R_PR.*T_PR./(V_PR-b_PR)-alpha_PR*a_PR./(V_PR.^2+2*b_PR.*V_PR-b_PR.^2);
PReosP= Function('P_pengRobinson',{V_PR,T_PR},{p_PR});


f=P_PR-R_PR.*T_PR./(VI_PR-b_PR)+alpha_PR*a_PR./(VI_PR.^2+2*b_PR.*VI_PR-b_PR.^2);
df=R_PR.*T_PR/(VI_PR-b_PR).^2-(alpha_PR*a_PR*(2*VI_PR+2*b_PR))./((VI_PR.^2+2*b_PR.*VI_PR-b_PR.^2).^2);
NR=VI_PR-f/df;
NR_fun= Function('NR_fun',{VI_PR,P_PR,T_PR},{NR});
VI_PR_guess(1)=R_PR.*T_PR./P_PR; %initial guess
ys = {};

for i=1:9
  VI_PR_guess(i+1) = NR_fun(VI_PR_guess(i),P_PR,T_PR);
end
final_V_PR=VI_PR_guess(end);
PReosV= Function('V_pengRobinson',{P_PR,T_PR},{final_V_PR});

%% Gas lift Oil well model

g=9.81;
mi=3.64e-3;
rho_l=760;
Mg=16.7;
Ta=348;
Va=64.34;
La=2048;
%Pgs=140*1.013e5;
Pgs=Dist(1);
Vt=25.03;
Sbh=0.0314;
Lbh=75;
Tt=369.4;
%GOR=0;
GOR=Dist(2);
%Pres=160*1.013e5;
Pres=Dist(3);
wres=y(4);
%wres=18;
Dt=0.134;
Lt=2048;
%PI=2.47e-6;
PI=Dist(4);

Kgs=9.98e-05;%altered
Kinj=1.40e-4;
Kpr=2.90e-3;
epsilon=4e-5;
P0=1.013e5;
u1=u(1);
u2=u(2);
%Pgs=u(3);
Db=Dt;

Pat=PReosP((Mg*Va)/(y(1)*1000),Ta);
Pab=Pat+y(1)*g*La/Va;
rho_gab=Mg/(1000*PReosV(Pab,Ta));%convert bar to Pa
rho_gin=Mg/(1000*PReosV(Pgs,Ta));%convert bar to Pa
w_gin=Kgs*u2*sqrt(rho_gin*max([Pgs-Pat,0]));

rho_gt=y(2)/(Vt+Sbh*Lbh - y(3)/rho_l);
Ptt = PReosP(Mg/(rho_gt*1000),Tt);
rho_mix=(y(2)+y(3)-rho_l*Sbh*Lbh)/Vt;
alpha_l=(y(3)-rho_l*Sbh*Lbh)/(Vt*rho_l);
alpha_mgb=GOR/(GOR+1);

U_slt=(4*(1-alpha_mgb)*wres)/(rho_l*pi*Dt*Dt);
U_sgt=(4*(w_gin+alpha_mgb*wres))/(rho_gab*pi*Dt*Dt);

U_mt=U_slt+U_sgt;

Ret=rho_mix*U_mt*Dt/mi;
lambda_t=(1/(-1.8*log10((epsilon/(Dt*3.7))^1.11+6.9/Ret)))^2;
Ft=(alpha_l*lambda_t*rho_mix*U_mt*U_mt*Lt)/(2*Dt);
Ptb = Ptt+rho_mix*Lt*g+Ft;
w_ginj=Kinj*sqrt(rho_gab*max([Pab-Ptb,0]));
U_lb=wres/(rho_l*Sbh);
Reb=(rho_l*U_lb*Db)/mi;
lambda_b=(1/(-1.8*log10((epsilon/(Db*3.7))^1.11+6.9/Reb)))^2;
Fb=(lambda_b*rho_l*U_lb*U_lb*Lbh)/(2*Db);
Pbh = Ptb+Fb+rho_l*g*Lbh;

w_res=PI*max([Pres-Pbh,0]);
w_lres=(1-alpha_mgb)*w_res;
w_gres=alpha_mgb*w_res;

rho_gtb=Mg/(1000*PReosV(Ptb,Tt));
alpha_lb=(w_lres*rho_gtb)/(w_lres*rho_gtb+(w_ginj+w_gres)*rho_l+1e-20);
alpha_lt=min([max([2*alpha_l-alpha_lb,0]),1]);

rho_mixt=alpha_lt*rho_l+(1-alpha_lt)*rho_gt;
w_out=Kpr*u1*sqrt(rho_mixt*max([Ptt-P0,0]));
Qout=w_out/rho_mixt;
alpha_mgt=(1-alpha_lt)*rho_gt/(alpha_lt*rho_l+(1-alpha_lt)*rho_gt);

w_gout= alpha_mgt*w_out;
w_lout=(1-alpha_mgt)*w_out;

dydt1=w_gin-w_ginj;
dydt2=w_ginj+w_gres- w_gout;
dydt3=w_lres-w_lout;
dydt4=(w_res-y(4))/10;

dydt=vertcat(dydt1, dydt2, dydt3, dydt4); %Derivatives

jac_dydt=jacobian(dydt,y); %Jacobian of the ODE


% This section groups the different variables of the ODE
W=vertcat(w_gin,w_ginj,w_res,w_lres,w_gres,w_out,w_gout,w_lout);
P=vertcat(Pat, Pab, Ptt,Ptb,Pbh);
rho=vertcat(rho_gab,rho_gin,rho_gt,rho_mix,rho_gtb,rho_mixt);
alpha=vertcat(alpha_l,alpha_mgb,alpha_lb,alpha_lt,alpha_mgt);
lambda=vertcat(U_mt,lambda_t,lambda_b,Ft,Fb);

%generate C-file for Mex function
opts = struct('main', true,...
              'mex', true);
fun_dydt=Function('derivatives',{y,u,Dist},{dydt,W,P,rho,alpha,lambda });
fun_jac_dydt=Function('jacobian',{y,u,Dist},{jac_dydt});

fun_dydt.generate('Well_mod_casadi.c',opts);
fun_jac_dydt.generate('Well_mod_casadi_jac.c',opts);

%Function returning the measurements
meas=vertcat(w_gin,w_out,w_gout,w_lout,Pat,  Ptt,rho_gin,rho_gt,rho_mixt,alpha_lt,alpha_mgt);

%Jacobian of the measurements, used for the ideal Extended Kalman Filter
jac_meas=jacobian(meas,y);

fun_meas=Function('measurements',{y,u,Dist},{meas});
fun_jac_meas=Function('jac_meas',{y,u,Dist},{jac_meas});

fun_meas.generate('Well_meas_casadi.c',opts);
fun_jac_meas.generate('Well_meas_casadi_jac.c',opts);

mex Well_mod_casadi.c -R2018a
mex Well_mod_casadi_jac.c -R2018a
mex Well_meas_casadi.c -R2018a
mex Well_meas_casadi_jac.c -R2018a
