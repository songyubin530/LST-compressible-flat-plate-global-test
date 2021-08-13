function [A,B,C,C_omega]=getABC(ny,alpha,beta,Re_x,gamma,Ma,Pr,u_stretch,T_stretch,Amu_stretch,dudy_stretch,dTdy_stretch,ddudy_stretch,ddTdy_stretch,dAmudT_stretch,dkdT_stretch,dlambdadT_stretch,ddAmudT_stretch,ddkdT_stretch)
A=zeros(5,5,ny);
B=zeros(5,5,ny);
C=zeros(5,5,ny);
C_omega=zeros(5,5,ny);
for j=1:ny
    A(1,1,j)=1;
    A(2,2,j)=1;
    A(4,4,j)=1;
    A(5,5,j)=1;
end

for j=1:ny
    B(1,1,j)=1/Amu_stretch(j)*dAmudT_stretch(j)*dTdy_stretch(j);
    B(1,2,j)=1i*alpha/Amu_stretch(j)*(Amu_stretch(j)-2/3*Amu_stretch(j));
    B(1,4,j)=1/Amu_stretch(j)*dAmudT_stretch(j)*dudy_stretch(j);
    
    B(2,1,j)=1i*alpha*(Amu_stretch(j)-2/3*Amu_stretch(j))/(2*Amu_stretch(j)-2/3*Amu_stretch(j));
    B(2,2,j)=1/(2*Amu_stretch(j)-2/3*Amu_stretch(j))*(2*dAmudT_stretch(j)+dlambdadT_stretch(j))*dTdy_stretch(j);
    B(2,3,j)=-Re_x/(2*Amu_stretch(j)-2/3*Amu_stretch(j));
    B(2,5,j)=1i*beta*(Amu_stretch(j)-2/3*Amu_stretch(j))/(2*Amu_stretch(j)-2/3*Amu_stretch(j));
            
    B(3,2,j)=1;
    
    B(4,1,j)=2*Amu_stretch(j)/Amu_stretch(j)*(gamma-1)*Ma^2*Pr*dudy_stretch(j);
    B(4,4,j)=2/Amu_stretch(j)*dkdT_stretch(j)*dTdy_stretch(j);
    
    B(5,2,j)=1i*beta*(Amu_stretch(j)-2/3*Amu_stretch(j))/Amu_stretch(j);
    B(5,5,j)=1/Amu_stretch(j)*dAmudT_stretch(j)*dTdy_stretch(j);
end

for j=1:ny
    C(1,1,j)=-1i*u_stretch(j)*alpha*Re_x/Amu_stretch(j)/T_stretch(j)-1/Amu_stretch(j)*(2*Amu_stretch(j)-2/3*Amu_stretch(j))*alpha^2-beta^2;
    C(1,2,j)=1i*alpha/Amu_stretch(j)*dAmudT_stretch(j)*dTdy_stretch(j)-Re_x/Amu_stretch(j)/T_stretch(j)*dudy_stretch(j);
    C(1,3,j)=-1i*alpha*Re_x/Amu_stretch(j);
    C(1,4,j)=1/Amu_stretch(j)*(ddAmudT_stretch(j)*dTdy_stretch(j)*dudy_stretch(j)+dAmudT_stretch(j)*ddudy_stretch(j));
    C(1,5,j)=-1/Amu_stretch(j)*(-2/3*Amu_stretch(j)+Amu_stretch(j))*alpha*beta;
    
    C(2,1,j)=1i*alpha/(2*Amu_stretch(j)-2/3*Amu_stretch(j))*dlambdadT_stretch(j)*dTdy_stretch(j);
    C(2,2,j)=Re_x*1i*(-u_stretch(j)*alpha)/(2*Amu_stretch(j)-2/3*Amu_stretch(j))/T_stretch(j)-Amu_stretch(j)*(alpha^2+beta^2)/(2*Amu_stretch(j)-2/3*Amu_stretch(j));
    C(2,4,j)=1/(2*Amu_stretch(j)-2/3*Amu_stretch(j))*dAmudT_stretch(j)*(1i*alpha*dudy_stretch(j));
    C(2,5,j)=1i*beta/(2*Amu_stretch(j)-2/3*Amu_stretch(j))*dlambdadT_stretch(j)*dTdy_stretch(j);
    
    C(3,1,j)=1i*alpha;
    C(3,2,j)=-dTdy_stretch(j)/T_stretch(j);
    C(3,3,j)=1i*(u_stretch(j)*alpha)*gamma*Ma^2;
    C(3,4,j)=1i*(-u_stretch(j)*alpha)/T_stretch(j);
    C(3,5,j)=1i*beta;
    
    C(4,2,j)=2*Amu_stretch(j)/Amu_stretch(j)*(gamma-1)*Ma^2*Pr*(1i*alpha*dudy_stretch(j))-Pr*Re_x/Amu_stretch(j)/T_stretch(j)*dTdy_stretch(j);
    C(4,3,j)=-1/Amu_stretch(j)*(gamma-1)*Ma^2*Pr*Re_x*1i*(-u_stretch(j)*alpha);
    C(4,4,j)=1/Amu_stretch(j)*ddkdT_stretch(j)*(dTdy_stretch(j))^2+1/Amu_stretch(j)*dkdT_stretch(j)*ddTdy_stretch(j)-(alpha^2+beta^2)+Pr*Re_x/Amu_stretch(j)/T_stretch(j)*1i*(-u_stretch(j)*alpha)+(gamma-1)*Pr*Ma^2/Amu_stretch(j)*dAmudT_stretch(j)*(dudy_stretch(j))^2;
    
    C(5,1,j)=-(Amu_stretch(j)-2/3*Amu_stretch(j))/Amu_stretch(j)*alpha*beta;
    C(5,2,j)=1i*beta/Amu_stretch(j)*dAmudT_stretch(j)*dTdy_stretch(j);
    C(5,3,j)=-1i*beta*Re_x/Amu_stretch(j);
    C(5,5,j)=1i*(-u_stretch(j)*alpha)*Re_x/Amu_stretch(j)/T_stretch(j)-alpha^2-(2*Amu_stretch(j)-2/3*Amu_stretch(j))/Amu_stretch(j)*beta^2;
end
for j=1:ny
    C_omega(1,1,j)=1i*Re_x/Amu_stretch(j)/T_stretch(j);
    C_omega(2,2,j)=Re_x*1i/(2*Amu_stretch(j)-2/3*Amu_stretch(j))/T_stretch(j);
    C_omega(3,3,j)=-1i*gamma*Ma^2;
    C_omega(3,4,j)=1i/T_stretch(j);
    C_omega(4,3,j)=-1/Amu_stretch(j)*(gamma-1)*Ma^2*Pr*Re_x*1i;
    C_omega(4,4,j)=Pr*Re_x/Amu_stretch(j)/T_stretch(j)*1i;
    C_omega(5,5,j)=1i*Re_x/Amu_stretch(j)/T_stretch(j);
end
end