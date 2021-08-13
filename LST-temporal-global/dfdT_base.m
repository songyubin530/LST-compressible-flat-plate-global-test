function [dAmudT_stretch,dkdT_stretch,dlambdadT_stretch,ddAmudT_stretch,ddkdT_stretch]=dfdT_base(y,y_stretch,T_reference,T_eta)
%根据Satherland公式，粘度的表达式为Amu/Amu_ref=(T/T_ref)^3/2*(1+S/T_ref)/(T/T_ref+S/T_ref),其中S=110.4;
%（1）那么dAmudT=((1+S/T_ref)*3/2*(T/T_ref)^1/2*(T/T_ref+S/T_ref)-(T/T_ref)^3/2*(1+S/T_ref))/(T/T_ref+S/T_ref)^2;
%考虑到Pr=CpAmu_\infty/k_\infty=CpAmu_\infty/k_\infty——Amu/Amu_\infty=k/k_\infty
%（2）从而dkdT=dAmudT=((1+S/T_ref)*3/2*(T/T_ref)^1/2*(T/T_ref+S/T_ref)-(T/T_ref)^3/2*(1+S/T_ref))/(T/T_ref+S/T_ref)^2;
%再考虑到lambda=-2/3Amu，lambda/Amu_\infty=-2/3Amu/Amu_\infty，lambda/(-3/2lambda_\infty)=-2/3Amu/Amu_\infty，lambda/lambda_\infty=Amu/Amu_\infty
%（3）从而dlambdadT=-2/3*dAmudT=-2/3*((1+S/T_ref)*3/2*(T/T_ref)^1/2*(T/T_ref+S/T_ref)-(T/T_ref)^3/2*(1+S/T_ref))/(T/T_ref+S/T_ref)^2;
%考虑到雷诺数的定义，其实lambda是用Amu无量纲的
%（4）ddAmudT=(1+S/T_ref)*(3/2*1/2*(T/T_ref)^(-1/2)*(T/T_ref+S/T_ref)-3/2*(T/T_ref)^1/2)/(T/T_ref+S/T_ref)^2-(1+S/T_ref)*(3/2*(T/T_ref)^1/2*(T/T_ref+S/T_ref)^2-2*(T/T_ref+S/T_ref)*(T/T_ref)^3/2)
dAmudT_init=((1+110.4/T_reference)*1.5.*T_eta.^0.5.*(T_eta+110.4/T_reference)-T_eta.^1.5*(1+110.4/T_reference))./((T_eta+110.4/T_reference)).^2;
dkdT_init=dAmudT_init;
dlambdadT_init=-2/3*dAmudT_init;
tem1=((1+110.4/T_reference)*1.5*0.5.*T_eta.^-0.5.*(T_eta+110.4/T_reference)-(1+110.4/T_reference)*1.5.*T_eta.^0.5)./(T_eta+110.4/T_reference).^2;
tem2=(1+110.4/T_reference).*(1.5.*T_eta.^0.5.*(T_eta+110.4/T_reference).^2-2.*(T_eta+110.4/T_reference).*T_eta.^1.5)./(T_eta+110.4/T_reference).^4;
ddAmudT_init=tem1-tem2;
ddkdT_init=ddAmudT_init;

dAmudT_stretch=interp1(y,dAmudT_init,y_stretch);
dkdT_stretch=interp1(y,dkdT_init,y_stretch);
dlambdadT_stretch=interp1(y,dlambdadT_init,y_stretch);
ddAmudT_stretch=interp1(y,ddAmudT_init,y_stretch);
ddkdT_stretch=interp1(y,ddkdT_init,y_stretch);
end