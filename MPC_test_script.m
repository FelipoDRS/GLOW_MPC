%this script tests the MPC

%at 3000 seconds, step to wgin =0.7 at 6000 seconds step to oil rate =20
% at step 9000 disturbance at Pres, = 150bar
global k
dt=10;
Hc=8;
Hp=8;
u=[0.5,0.5];
D=zeros(4,1);
XI=[3000,50,8000,10]; %initial state
lb=zeros(2*Hc,1);
ub=ones(2*Hc,1);
A=zeros(4*Hc,2*Hc);
b=[0.4*ones(2*Hc,1);-0.4*ones(2*Hc,1)];
P0=1.013e5;%output pressure
LS=27000; %length of simulation
par=[5.83e-3 14.45 3.29 1.20];
f=@(x)( (sim(net,((x-MNX)./STDX)')).*STDY'+MNY');
OPS=optimset('Algorith','sqp','MaxIter',15,'Display','off','tolX',1e-3,'tolF',1e-3);
for i=1:(Hc-1)
A(2*i-1,2*i-1)=1;
A(2*i-1,2*i+1)=-1;
A(2*i,2*i)=1;
A(2*i,2*i+2)=-1;
A(2*Hc+2*i-1,2*i-1)=-1;
A(2*Hc+2*i-1,2*i+1)=1;
A(2*Hc+2*i,2*i)=-1;
A(2*Hc+2*i,2*i+2)=1;
end
D(1)=140e5+randn(1)*1e5;
    D(2)=max([0,randn(1)-2]);
    %D(2)=0;
    D(3)=160e5+randn(1)*1e5;
    t=0;
    D(4)=2.47e-6;

[~,x]=ode15s(@(t,y)(well_disturbance_der(t,y,u,D)),[0,dt],XI);
[w_meas,P_meas,rho_meas,alpha_meas]=monitor_well_disturbance( t,x,u, D );%meas stands for measured
%those are list used to save the results
time=0:dt:LS;
u_list=ones(LS/dt,2)*diag(u);
D_list=ones(LS/dt,4)*diag(D);
x_list=ones(LS/dt,4)*diag(XI);
meas_list=ones(LS/dt,15)*diag([w_meas,w_meas,P_meas,rho_meas,alpha_meas]);
sp_list=ones(LS/dt,2)*diag([0.86, 27.5]);
x0_list=ones(LS/dt,3);
time_list=zeros(LS/dt,1);
obj_list=zeros(LS/dt,1);

for i=1:ceil(LS/dt)
    k=i;
    winsp= 0.5;
    wooilsp=21;
    D(1)=140e5+randn(1)*1e5;
    D(2)=max([0,randn(1)-2]);
    %D(2)=0;
    D(3)=160e5+randn(1)*1e5;
    D(4)=2.47e-6;

    if i>9000/dt && i<12000/dt
    winsp=0.8;
    end
    
    if i>13000/dt && i<16000/dt
    wooilsp=25;
    end
    
    if i>17000/dt && i<20000/dt
    D(3)=150e5+randn(1)*1e4;
    end
    
    if i>21000/dt && i<24000/dt
    D(4)=3e-6;
    end
    
    ysp=ones(Hp,2)*diag([winsp, wooilsp]);    
    x1=P_meas(1)*16.7*64.34/(8314*348);

    u0=ones(Hc,2)*diag([wooilsp/((1-alpha_meas(end))*2.9e-3*sqrt(rho_meas(3)*(P_meas(2)-P0))),winsp/(9.98e-5*sqrt(rho_meas(1)*(D(1)-P_meas(1))))]);%initial guess is the solution assuming constant pressure
    u0=reshape(u0',Hc*2,1);
    x23=f([meas_list(i,:),D(1),u_list(i,:)]);
    x0=[x1,x23'];
    lb(1)=max([0,u_list(i,1)-0.4]);
    ub(1)=min([1,u_list(i,1)+0.4]);
    lb(2)=max([0,u_list(i,2)-0.4]);
    ub(2)=min([1,u_list(i,2)+0.4]);
    
    mpc_eq=@(U)(MPC_obj(x0,u_list(i,:),meas_list(i,[5,8]),ysp,dt,Hc,Hp,par,U,D(1)));
    t1=tic;
    [U,obj]=fmincon(mpc_eq,u0,A,b,[],[],lb,ub,[],OPS);
    t2=toc(t1);
    u=U([1,2]);

    [~,x]=ode15s(@(t,y)(well_disturbance_der(t,y,u,D)),[0,dt],x_list(i,:));    
    [w_meas,P_meas,rho_meas,alpha_meas]=monitor_well_disturbance( t,x(end,:),u, D );%meas stands for measured

    [~,w_check,P_check,rho_check,alpha_check,~]=well_disturbance(t,x(end,:),u,D);
u_list(i+1,:)=u;
D_list(i+1,:)=D;
x_list(i+1,:)=x(end,:);
meas_list(i+1,5:end)=[w_meas,P_meas,rho_meas,alpha_meas];
sp_list(i+1,:)=[winsp, wooilsp]; 
x0_list(i+1,:)=x0; 
time_list(i+1)=t2;
obj_list(i+1)=obj;

if i>6
        meas_list(i+1,1:4)=mean(meas_list((i-5):i,5:8));
 else
     meas_list(i+1,1:4)=w_meas;
 end
   if mod(i,350) ==0
       i
   end
end