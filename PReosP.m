function p = PReosP(V,T)
Tc=190.6;
Pc=46.04*1.013e5;
w=0.011;

R=8.314;
a=0.45724*R^2*Tc^2/Pc;
b=0.07780*R*Tc/Pc;
         
k=0.3764+1.54226*w-0.26992*w^2;

alpha=(1+k*(1-sqrt(T./Tc))).^2;
p=R.*T./(V-b)-alpha*a./(V.^2+2*b.*V-b.^2);
end