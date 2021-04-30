

%Example script to show how to use the model

%First the initial states, gas mass in the annulus,gas mass in the tube,oil mass in the tube,velocity of the gas liquid mixture
%in that order
XI=[3000,100,8000,20];
% The valve positions, 80% gas production valve opening, 40% gas injection
% valve opening
u=[0.8,0.4];

%Finally the disturbances
%The disturbances are, source gas pressure, Gas Oil Ratio, Pressure in the
%reservoir and productivity index
D=[140e5;0;160e5;2.47e-6];

% empty variable t
t=0;

%total integration time
t_total=3600;
%time step
dt=10;

steps=ceil(t_total/dt);
%preallocating the matrix that will hold the data
state_array=zeros(steps+1,4);
flow_array=zeros(steps+1,4);
pressure_array=zeros(steps+1,2);
t_vector=0:dt:t_total;

%monitor_well_disturbance returns the measurements on the top side of the well 
[w_meas,P_meas,rho_meas,alpha_meas]=monitor_well_disturbance( t,XI,u, D );

state_array(1,:)=XI;
flow_array(1,:)=w_meas;
pressure_array(1,:)=P_meas;

%  @(t,y)(well_disturbance_der(t,y,u,D)) creates an anonymous function that
%  will feed the derivatives to the integrator, in this case ode15s

%the integration shoud be run inside a loop
for i= 1:steps 
% that returns the derivatives used in the integration
[t,x]=ode15s(@(t,y)(well_disturbance_der(t,y,u,D)),[0,dt],XI);

[w_meas,P_meas,rho_meas,alpha_meas]=monitor_well_disturbance( t,x,u, D );
XI=x(end,:);
state_array(i+1,:)=XI;
flow_array(i+1,:)=w_meas;
pressure_array(i+1,:)=P_meas;
end

figure()
plot(t_vector,state_array)
legend({'gas mass in the annulus','gas mass in the tube','oil mass in the tube','velocity of the gas liquid mixture'})

figure()
plot(t_vector,flow_array(:,1),t_vector,flow_array(:,4))
legend({'gas inlet flow','oil outlet flow'})
figure()
plot(t_vector,pressure_array)
legend({'gas inlet pressure','oil outlet pressure'})


%Now doing a step change on the valves openings 60% gas production valve opening, 25% gas injection
% valve opening
u=[0.5,0.20];
t_total=3600*9;
%time step
dt=5;

steps=ceil(t_total/dt);

state_array2=zeros(steps,4);
flow_array2=zeros(steps,4);
pressure_array2=zeros(steps,2);
t_vector2=t_vector(end):dt:(t_vector(end)+t_total);


state_array2(1,:)=XI;
flow_array2(1,:)=w_meas;
pressure_array2(1,:)=P_meas;

%  @(t,y)(well_disturbance_der(t,y,u,D)) creates an anonymous function that
%  will feed the derivatives to the integrator, in this case ode15s

%this shoud be run inside a loop
for i= 1:steps 
% that returns the derivatives used in the integration
[t,x]=ode15s(@(t,y)(well_disturbance_der(t,y,u,D)),[0,dt],XI);

%monitor_well_disturbance returns the measurements on the top side of the well 
[w_meas,P_meas,rho_meas,alpha_meas]=monitor_well_disturbance( t,x,u, D );
XI=x(end,:);
state_array2(i+1,:)=XI;
flow_array2(i+1,:)=w_meas;
pressure_array2(i+1,:)=P_meas;
end
state_array3=[state_array;state_array2];
flow_array3=[flow_array;flow_array2];
pressure_array3=[pressure_array;pressure_array2];
t_vector3=[t_vector,t_vector2];


figure()
plot(t_vector3,state_array3)
legend({'gas mass in the annulus','gas mass in the tube','oil mass in the tube','velocity of the gas liquid mixture'})

figure()
plot(t_vector3,flow_array3(:,1),t_vector3,flow_array3(:,4))
legend({'gas inlet flow','oil outlet flow'})
figure()
plot(t_vector3,pressure_array3)
legend({'gas inlet pressure','oil outlet pressure'})
