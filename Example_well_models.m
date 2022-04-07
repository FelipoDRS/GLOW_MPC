XI=[3000,80,5000,10];
par=[0,16,2.47,1.4];
Dist_well=[140e5,0,160e5,2.47e-6];
Dist_mpc=140e5;
dt=40;

u=[1,0.8];
n_points=3000;
MPC_meas_list=zeros(n_points,11);
MPC_state_list=zeros(n_points,3);
Well_meas_list=zeros(n_points,11);
Well_state_list=zeros(n_points,4);
MPC_meas_list(1,:)=MPC_meas_casadi(XI(1:3),par,u,Dist_mpc);
MPC_state_list(1,:)=XI(1:3);
Well_meas_list(1,:)=Well_meas_casadi(XI,u,Dist_well);
Well_state_list(1,:)=XI;

for i=2:n_points

    if i>1000 & i <2000
    u=[1,0.5];
    end
    if i>2000 
    u=[1,0.2];
    end
    
    [~,x_mpc]=ode45(@(t,y) MPC_mod_casadi(y,par,u,Dist_mpc),[0,dt],MPC_state_list(i-1,:));
    [~,x_well]=ode45(@(t,y) Well_mod_casadi(y,u,Dist_well),[0,dt],Well_state_list(i-1,:));
MPC_state_list(i,:)=x_mpc(end,:);
Well_state_list(i,:)=x_well(end,:);
    
MPC_meas_list(i,:)=MPC_meas_casadi(x_mpc(end,:),par,u,Dist_mpc);
Well_meas_list(i,:)=Well_meas_casadi(x_well(end,:),u,Dist_well);
    
end

figure()
subplot(2,1,1)
plot((1:length(MPC_meas_list))*dt/3600,MPC_meas_list(:,1))
xlabel('Time (h)')
ylabel('Gas inlet flow rate (kg/s)')
hold on
plot((1:length(Well_meas_list))*dt/3600,Well_meas_list(:,1))
hold off
%xlim([0,lim_x])
legend({'MPC model','Well model'})
subplot(2,1,2)
plot((1:length(MPC_meas_list))*dt/3600,MPC_meas_list(:,4))
xlabel('Time (h)')
ylabel('Oil outlet flow rate (kg/s)')
hold on
plot((1:length(Well_meas_list))*dt/3600,Well_meas_list(:,4))
hold off
legend({'MPC model','Well model'})

