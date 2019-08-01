function [ S ] = Par_est_obj( X0,par,N )
%This one use local function instead of using a defined here
%I hope this one is faster
%Update: It is. Also variables are now scaled, do not forget it.
D_list=[];
Y=[];
sensor_data=[];
load('data_for_param_est.mat');
if N==0
    N=length(sensor_data);
end
%N= 15000;%length(sensor_data);%number of data points used
sensor_data=sensor_data(1:N,:);
Y_mod=zeros(length(sensor_data),3);
Sensor_estimado=zeros(size(sensor_data));
ODS=odeset('MaxStep', 100);
X0(1)=X0(1)*1e3;
X0(2)=X0(2)*10;
X0(3)=X0(3)*1e3;

global u Dist

u=sensor_data(1,[end-1,end]);
Dist=sensor_data(1,end-2);
[t,ymodelo]=ode15s(@(t,y) MPC_model(t,y,par,u,Dist),[sensor_data(1,1),sensor_data(2,1)],X0,ODS);
Y_mod(1,:)=ymodelo(end,:);
[W,P,rho,alpha]=MPC_monitor(0,ymodelo(end,:),par,u,Dist);
temp1=[W,P,rho,alpha];%to de saco
%cheio de nomear variaveis

Sensor_estimado(1,:)=[sensor_data(1,1),temp1,Dist,u];

for i=2:length(sensor_data)
u=sensor_data(i,[end-1,end]);
Dist=sensor_data(i,end-2);
[~,ymodelo]=ode23t(@(t,y) MPC_model(t,y,par,u,Dist),[sensor_data(i-1,1),sensor_data(i,1)] ,ymodelo(end,:),ODS);
Y_mod(i,:)=ymodelo(end,:);
[W,P,rho,alpha]=MPC_monitor(0,ymodelo(end,:),par,u,Dist);
temp1=[W,P,rho,alpha];%to de saco
%cheio de nomear variaveis

Sensor_estimado(i,:)=[sensor_data(i,1),temp1,Dist,u];
if mod(i,800)==1
i;
end

end
col=size(Sensor_estimado(:,2:end-3),2);
S=sum(sum((Sensor_estimado(:,2:end-3)-sensor_data(:,2:end-3)).^2)./var(sensor_data(:,2:end-3)))/(col*N);
end

