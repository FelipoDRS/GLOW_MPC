function [ error ] = MPC_obj( x0,u0,y0,ysp,dt,Hc,Hp,par,U_fut,Pgs )
%Receives [u(t+1);u(t+2) ... ;u(t+Hc)]
%Hc is control horizon
%Hp is prediction horizon
%x0 and u0 are the initial state and control position
%dt is delta t
%ysp are the setpoints
global u
global Dist k
u_fut=[U_fut(1:2:end),U_fut(2:2:end)];

ODS=odeset('MaxStep', 100);
Dist=Pgs;
x=zeros(Hp+1,length(x0));
Wlist=zeros(Hp,2);
x(1,:)=x0;
t=0;
[W1,~,~,~]=MPC_monitor(0,x0,par,u,Dist);
C=W1([1,4])-y0; %corrector
for i=1:Hp
    if i <= Hc
    u=u_fut(i,:);
    else
    u=u_fut(end,:);
    end
    [~,Y_mpc]=ode15s(@(t,y) MPC_model(t,y,par,u,Dist),[t:t+dt] ,x(i,:),ODS);
    [W,~,~,~]=MPC_monitor(0,Y_mpc(end,:),par,u,Dist);
    x(i+1,:)=Y_mpc(end,:);
    Wgin=W(1);
    Woout=W(4);
    Wlist(i,:)=[Wgin, Woout];
end
WlistC=Wlist-ones(size(Wlist))*diag(C);
du=[u0;u_fut(1:end-1,:)]-u_fut;
dy=[y0; Wlist(1:end-1,:)]-Wlist(1:end,:);
errorSP1=(WlistC-ysp).*(WlistC-ysp)*diag([6,5]);
errorSP2=sum(sum(errorSP1(1:end-1,:)));
errorDu=sum(sum(du.*du*diag([1,1])));
errorDy=sum(sum(dy.*dy*diag([1,1])));
error= errorSP2+errorDy+errorDu;

%error =sum(sum((Wlist-ysp).*(Wlist-ysp)*diag([2.5,10])));% +du.*du*diag([0.1,0.1])+dy.*dy*diag([0.1,0.1])));
end

