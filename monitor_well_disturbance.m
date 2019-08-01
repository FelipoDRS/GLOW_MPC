function [  W,P,rho,alpha ] = monitor_well_Disturbance( t,y,u, D ) 

[~,W,P,rho,alpha,~]= well_disturbance(t,y,u, D);

W=max(W+randn(size(W))*0.1,0);
P=P+randn(size(P))*1e4;
rho=rho+randn(size(rho))*0.05;
alpha=max(0,alpha+randn(size(alpha))*0.01);


W=W([1,6,7,8]);
P=P([1,3]);
rho=rho([2,3,6]);
alpha=alpha([4,5]);
end