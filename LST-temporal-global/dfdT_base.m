function [dAmudT_stretch,dkdT_stretch,dlambdadT_stretch,ddAmudT_stretch,ddkdT_stretch]=dfdT_base(y,y_stretch,T_reference,T_eta)
%����Satherland��ʽ��ճ�ȵı��ʽΪAmu/Amu_ref=(T/T_ref)^3/2*(1+S/T_ref)/(T/T_ref+S/T_ref),����S=110.4;
%��1����ôdAmudT=((1+S/T_ref)*3/2*(T/T_ref)^1/2*(T/T_ref+S/T_ref)-(T/T_ref)^3/2*(1+S/T_ref))/(T/T_ref+S/T_ref)^2;
%���ǵ�Pr=CpAmu_\infty/k_\infty=CpAmu_\infty/k_\infty����Amu/Amu_\infty=k/k_\infty
%��2���Ӷ�dkdT=dAmudT=((1+S/T_ref)*3/2*(T/T_ref)^1/2*(T/T_ref+S/T_ref)-(T/T_ref)^3/2*(1+S/T_ref))/(T/T_ref+S/T_ref)^2;
%�ٿ��ǵ�lambda=-2/3Amu��lambda/Amu_\infty=-2/3Amu/Amu_\infty��lambda/(-3/2lambda_\infty)=-2/3Amu/Amu_\infty��lambda/lambda_\infty=Amu/Amu_\infty
%��3���Ӷ�dlambdadT=-2/3*dAmudT=-2/3*((1+S/T_ref)*3/2*(T/T_ref)^1/2*(T/T_ref+S/T_ref)-(T/T_ref)^3/2*(1+S/T_ref))/(T/T_ref+S/T_ref)^2;
%���ǵ���ŵ���Ķ��壬��ʵlambda����Amu�����ٵ�
%��4��ddAmudT=(1+S/T_ref)*(3/2*1/2*(T/T_ref)^(-1/2)*(T/T_ref+S/T_ref)-3/2*(T/T_ref)^1/2)/(T/T_ref+S/T_ref)^2-(1+S/T_ref)*(3/2*(T/T_ref)^1/2*(T/T_ref+S/T_ref)^2-2*(T/T_ref+S/T_ref)*(T/T_ref)^3/2)
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