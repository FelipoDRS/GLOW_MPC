function V = PReosV(P,T)
Tc=190.6;
Pc=46.04*1.013e5;
w=0.011;

R=8.314;
a=0.45724*R^2*Tc^2/Pc;
b=0.07780*R*Tc/Pc;
         
k=0.3764+1.54226*w-0.26992*w^2;

alpha=(1+k*(1-sqrt(T./Tc))).^2;
VI=R.*T./P; %initial guess
e=1;
while e>1e-6
    f=P-R.*T./(VI-b)+alpha*a./(VI.^2+2*b.*VI-b.^2);
    df=R.*T/(VI-b).^2-(alpha*a*(2*VI+2*b))./((VI.^2+2*b.*VI-b.^2).^2);
    V=VI-f/df;
    e=abs((V-VI))/V;
    VI=V;
end
end