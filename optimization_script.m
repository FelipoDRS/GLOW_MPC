%optimization script
N=3000;
f=@(x) Par_est_obj2(x(1:3),x(4:7),N);
OPS=optimset('Display','iter','MaxIter',120);
[a,b,c,d,e]=fminsearch3(f,[1.7,5,10,0.01,16,2.5,10],OPS)
%resposta é a=[3,5,8,0,16,2.47,1.4]
f2=@(x) Par_est_obj_LS(x(1:3),x(4:7));
opt = optimoptions(@lsqnonlin,'Display','iter-detailed','MaxIter',20,'Algorithm','levenberg-marquardt','InitDamping',0);
[x,resnorm,residual,exitflag,output,lambda,jacobian]  = lsqnonlin(f2,[1,5,10,0.01,16,2.5,10],[],[],opt)
a=zeros(1,7);
load data_for_param_est.mat
a(1)=sensor_data(1,6)*16.7*64.34/(8314*348*1000);
a(2:7)=[5,10,0.01,16,2.5,10]+randn(1,6);
a(4)=4.4e-3;
E=[];

for i=1:4
    N=i*6000;
   f=@(x) Par_est_obj2(x(1:3),x(4:7),N);

    OPS=optimset('Display','iter','MaxIter',i*4);
 
   [a,b,~,~,e]=fminsearch3(f,a,OPS);
   e=e(sum(e')~=0,:);
   E=[E;e];
end
OPS=optimset('Display','iter','MaxIter',600);
 N=00;
 f=@(x) Par_est_obj2(x(1:3),x(4:7),N);
   [a,b,~,~,e]=fminsearch3(f,a,OPS);
   E=[E;e];
%this calculate the numerical hessian and jacobian

a1=[3.34	8.31	10.08 5.82e-3	14.45	3.29	1.20];
dh=1e-4;
H=zeros(7);
jac=zeros(7,1);
for i=1:7
daF=a1;
daB=a1;
daF(i)=daF(i)+dh;
daB(i)=daB(i)-dh;
jac(i)=(Par_est_obj2(daF(1:3),daF(4:7),0)-Par_est_obj2(daB(1:3),daB(4:7),00))/(2*dh);
end

for i=1:6
    for j=i:6
        da1=a1;
da2=a1;
da3=a1;
da4=a1;

da1(i)=da1(i)+dh;
da1(j+1)=da1(j+1)+dh;
da2(i)=da2(i)-dh;
da2(j+1)=da2(j+1)+dh;
da3(i)=da3(i)+dh;
da3(j+1)=da3(j+1)-dh;
da4(i)=da4(i)-dh;
da4(j+1)=da4(j+1)-dh;


        H(i,j+1)=(Par_est_obj2(da1(1:3),da1(4:7),000)+Par_est_obj2(da4(1:3),da4(4:7),000)...
            -Par_est_obj2(da2(1:3),da2(4:7),00)-Par_est_obj2(da3(1:3),da3(4:7),000))/(4*dh*dh);
        
        H(j+1,i)=H(i,j+1);
    end
    i
end
F0=Par_est_obj2(a1(1:3),a1(4:7),00);
for i=1:7
    da1=a1;
    da2=a1;
    da1(i)=da1(i)+dh;
    da2(i)=da2(i)-dh;
    H(i,i)=(Par_est_obj2(da1(1:3),da1(4:7),00)+Par_est_obj2(da2(1:3),da2(4:7),000)-2*F0)/(4*dh*dh);
end

a1=[2.97,7.94,7.47,1.9e-4,17.99,5,3.36];
par=a1(4:7);
dh=5e-4;
H=zeros(3);
jac=zeros(4,1);
for i=1:3
daF=a1(1:3);
daB=a1(1:3);
daF(i)=daF(i)+dh;
daB(i)=daB(i)-dh;
jac(i)=(Par_est_obj2(daF,par,0)-Par_est_obj2(daB,par,00))/(2*dh);
end

for i=1:2
    for j=i:2
        da1=a1;
da2=a1;
da3=a1;
da4=a1;

da1(i)=da1(i)+dh;
da1(j+1)=da1(j+1)+dh;
da2(i)=da2(i)-dh;
da2(j+1)=da2(j+1)+dh;
da3(i)=da3(i)+dh;
da3(j+1)=da3(j+1)-dh;
da4(i)=da4(i)-dh;
da4(j+1)=da4(j+1)-dh;


        H(i,j+1)=(Par_est_obj2(da1,par,1000)+Par_est_obj2(da4,par,1000)-Par_est_obj2(da2,par,1000)-Par_est_obj2(da3,par,1000))/(4*dh*dh);
        
        H(j+1,i)=H(i,j+1);
    end
    i
end
F0=Par_est_obj2(a1(1:3),par,1000);
for i=1:3
    da1=a1(1:3);
    da2=a1(1:3);
    da1(i)=da1(i)+dh;
    da2(i)=da2(i)-dh;
    H(i,i)=(Par_est_obj2(da1,par,1000)+Par_est_obj2(da2,par,1000)-2*F0)/(4*dh*dh);
end
