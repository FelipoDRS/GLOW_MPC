
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
%Recommended disturbance=140e5
%Recommended parameters=[0,16,2.47,1.4]


% Limits to ensure the functioning of the model:

%x2+x3 > 1790
%x3<20812 & x3>1789


import casadi.*

y = SX.sym('y',3);
par=SX.sym('par',4);
Dist=SX.sym('pgs');
u=SX.sym('u',2);

Pgs=Dist(1);
GOR=par(1);
Pres=par(2)*1e6;
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
Vt=25.03;
Sbh=0.0314;
Lbh=75;
Tt=369.4;
Dt=0.134;
Lt=2048;
Kgs=9.98e-05;
Kpr=2.90e-3;
P0=1.013e5;
u1=u(1);
u2=u(2);

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
alpha_lb=(w_lres*rho_gtb)/(w_lres*rho_gtb+(w_ginj+w_gres)*rho_l);
alpha_lt=min([max([2*alpha_l-alpha_lb,0]),1]);

rho_mixt=alpha_lt*rho_l+(1-alpha_lt)*rho_gt;
w_out=Kpr*u1*sqrt(rho_mixt*max([Ptt-P0,0]));
Qout=w_out/rho_mixt;
alpha_mgt=(1-alpha_lt)*rho_gt/(alpha_lt*rho_l+(1-alpha_lt)*rho_gt);

w_gout= alpha_mgt*w_out;
w_lout=(1-alpha_mgt)*w_out;
%ODE
dydt1=w_gin-w_ginj;
dydt2=w_ginj+w_gres- w_gout;
dydt3=w_lres-w_lout;

dydt=vertcat(dydt1 , dydt2 , dydt3); %ODE

jac_dydt=jacobian(dydt,y); %Jacobian of the ODE

% This section groups the different variables of the ODE
W=vertcat(w_gin,w_ginj,w_res,w_lres,w_gres,w_out,w_gout,w_lout);
P=vertcat(Pat, Pab, Ptt,Ptb,Pbh);
rho=vertcat(rho_gab,rho_gin,rho_gt,rho_mix,rho_gtb,rho_mixt);
alpha=vertcat(alpha_l,alpha_mgb,alpha_lb,alpha_lt,alpha_mgt);

fun_dydt=Function('MPC_ODE',{y,par,u,Dist},{dydt,W,P,rho,alpha});
fun_jac_dydt=Function('jacobian',{y,par,u,Dist},{jac_dydt});

%generate C-file for Mex function
opts = struct('main', true,...
              'mex', true);
fun_dydt.generate('MPC_mod_casadi.c',opts);
fun_jac_dydt.generate('MPC_mod_casadi_jac.c',opts);

meas=vertcat(w_gin,w_out,w_gout,w_lout,Pat,Ptt,rho_gin,rho_gt,rho_mixt,alpha_lt,alpha_mgt);
jac_meas=jacobian(meas,y);

%Function returning the measurements

fun_meas=Function('measurements',{y,par,u,Dist},{meas});
%Jacobian of the measurements, used for the Extended Kalman Filter
fun_jac_meas=Function('jacobian_meas',{y,par,u,Dist},{jac_meas});

opts = struct('main', true,...
              'mex', true);
fun_meas.generate('MPC_meas_casadi.c',opts);
fun_jac_meas.generate('MPC_meas_casadi_jac.c',opts);

mex MPC_mod_casadi.c -R2018a
mex MPC_meas_casadi.c -R2018a
mex MPC_mod_casadi_jac.c -R2018a
mex MPC_meas_casadi_jac.c -R2018a

