function [ W,P,rho,alpha ] = MPC_monitor(t, y,par,u,Dist )

[~,W,P,rho,alpha]=MPC_model(t,y,par,u,Dist);

W=W([1,6,7,8]);
P=P([1,3]);
rho=rho([2,3,6]);
alpha=alpha([4,5]);



end
