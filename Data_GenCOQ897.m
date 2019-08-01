%this script is an example of how to use the model with disturbance
%I use it in case I forget how to run it

u=[0.55,0.55];
D=[140*1.013e5,0,160*1.013e5 2.47e-6];
A=flip((fullfact([3,3])-0.4)/3.2);%make the valve go to 0.2,0.5 and 8
rng(3720);
dt=5;
TS=9*4*3600;
Y=zeros(ceil(TS/dt),4);
T=zeros(ceil(TS/dt),1);
W=zeros(ceil(TS/dt),4);
P=zeros(ceil(TS/dt),2);
rho=zeros(ceil(TS/dt),3);
alpha= zeros(ceil(TS/dt),2);
u_list=zeros(ceil(TS/dt),2);
D_list=zeros(ceil(TS/dt),4);
Y(1,:)=[3000,50,8000,10];
[w,p,r,al ] = monitor_well_disturbance(0,Y(1,:),u,D) ;
W(1,:)=w;
P(1,:)=p;
rho(1,:)=r;
alpha(1,:)=al;
u_list(1,:)=u;
D_list(1,:)=D;
t=0;
j=1;
u=A(j,:);
for i=2:ceil(TS/dt)
    [t2,y]=ode15s(@(t,y) well_disturbance_der(t,y,u,D),[t,t+dt],Y(i-1,:));
    T(i)=t2(end);
    Y(i,:)=y(end,:);
    [w,p,r,al ] = monitor_well_disturbance(0,Y(i,:),u,D);
    W(i,:)=w;
    P(i,:)=p;
    rho(i,:)=r;
    alpha(i,:)=al;
    u_list(i,:)=u;
    D_list(i,:)=D;
    t=t+dt;
    D(1)=140*1.013e5+randn(1)*1e5;
    D(2)=max([0,randn(1)-2]);
    D(3)=160e5+randn(1)*1e5;
    D(4)=2.47e-6;
    if mod(i,4*3600/dt)==1
        j=j+1
        u=A(j,:);
    end


end


sensor_data=[T,W,P,rho,alpha,D_list,u_list];

%save 'data_for_param_est.mat' sensor_data Y D_list