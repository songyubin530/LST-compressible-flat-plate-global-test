function [y,d_eta,u_eta,v_eta,T_eta,Amu_eta]=shooting_method_isothermy(Ma,Pr,gamma,T_reference,nx,SLX,Re_x,Tw)
eta=linspace(0,SLX,nx);
s1=0.5;
s2=15;
epsilon1=1;
epsilon2=1;
RK=4;%采用四阶龙格库塔
%定义初始数组
f=zeros(nx,9);
g=zeros(nx,6);
fn=zeros(1,9);%传递量
gn=zeros(1,6);%传递量
kf=zeros(9,4);%传递量
kg=zeros(6,4);%传递量
% k1=zeros(9,1);%最终量
% g1=zeros(6,1);%最终量
step=0;%迭代步数

while abs(epsilon1)>10e-5 || abs(epsilon2)>10e-5
%赋初值
f(1,1)=0;f(1,2)=0;f(1,3)=s1;f(1,4)=0;f(1,5)=0;
f(1,6)=1;f(1,7)=0;f(1,8)=0;f(1,9)=0;
g(1,1)=Tw;g(1,2)=s2;g(1,3)=0;
g(1,4)=0;g(1,5)=0;g(1,6)=1;
for i=1:nx-1
    h=eta(i+1)-eta(i);
    for j=1:RK
            switch j
                case 1
                fn=f(i,:);
                gn=g(i,:);
                case 2
                fn=f(i,:)+0.5*h*kf(:,j-1)';
                gn=g(i,:)+0.5*h*kg(:,j-1)';
                case 3
                fn=f(i,:)+0.5*h*kf(:,j-1)';
                gn=g(i,:)+0.5*h*kg(:,j-1)';
                case 4
                fn=f(i,:)+h*kf(:,j-1)';
                gn=g(i,:)+h*kg(:,j-1)'; 
            end  
[C1,dC1]=get_C1(gn(1),gn(2),T_reference);
[A1,A2,A3,A4]=get_A(fn(1),fn(3),fn(4),fn(6),fn(7),fn(9),gn(1),gn(2),gn(3),gn(4),gn(5),gn(6),C1,dC1,T_reference,Ma,Pr,gamma);
kf(1,j)=fn(2);
kf(2,j)=fn(3);
kf(3,j)=-(dC1*fn(3)+fn(1)*fn(3))/C1;
kf(4,j)=fn(5);
kf(5,j)=fn(6);
kf(6,j)=A1;
kf(7,j)=fn(8);
kf(8,j)=fn(9);
kf(9,j)=A2;

kg(1,j)=gn(2);
temp=(1-gamma)*Ma^2*Pr;
kg(2,j)=(temp*C1*fn(3)^2-dC1*gn(2)-Pr*fn(1)*gn(2))/C1;
kg(3,j)=gn(4);
kg(4,j)=A3;
kg(5,j)=gn(6);
kg(6,j)=A4;
    end
  k1=kf(:,1)+2*kf(:,2)+2*kf(:,3)+kf(:,4);
  g1=kg(:,1)+2*kg(:,2)+2*kg(:,3)+kg(:,4);
  
  f(i+1,:)=f(i,:)+h/6*k1';
  g(i+1,:)=g(i,:)+h/6*g1';
end
ds1=((1-f(nx,2))/f(nx,8)-(1-g(nx,1))/g(nx,5))/(f(nx,5)/f(nx,8)-g(nx,3)/g(nx,5));
ds2=((1-f(nx,2))/f(nx,5)-(1-g(nx,1))/g(nx,3))/(f(nx,8)/f(nx,5)-g(nx,5)/g(nx,3));
s1=s1+ds1;
s2=s2+ds2;
epsilon1=1-f(nx,2);
epsilon2=1-g(nx,1);
step=step+1;
% fprintf('step= %d \n',step)
end

%从eta转换到坐标y-----------------------------------------------------------
y_tem=zeros(nx,1);
for i=2:nx
    y_tem(i,1)=y_tem(i-1,1)+(g(i,1)+g(i-1,1))*(eta(i)-eta(i-1))/2;
end
y=sqrt(2)*y_tem;
%--------------------------------------------------------------------------

d_eta=1./g(:,1);
u_eta=f(:,2);
v_eta=-1/Re_x/sqrt(2)*(f(:,1).*g(:,1)-f(:,2).*y_tem);
T_eta=g(:,1);
Amu_eta=zeros(nx,1);
for i=1:nx
    Amu_eta(i,1)=T_eta(i)^(3/2)*(1+110.4/T_reference)/(T_eta(i)+110.4/T_reference);
end
end

function [C1,dC1]=get_C1(g1,g2,T_eta)
C1=sqrt(g1)*(1+110.4/T_eta)/(g1+110.4/T_eta);
temp1=0.5/sqrt(g1)*(1+110.4/T_eta)/(g1+110.4/T_eta);
temp2=sqrt(g1)*(1+110.4/T_eta)/(g1+110.4/T_eta)^2;
dC1=(temp1-temp2)*g2;
end

function [A1,A2,A3,A4]=get_A(f1,f3,f4,f6,f7,f9,g1,g2,g3,g4,g5,g6,C1,dC1,T_eta,Ma,Pr,gamma)
temp1=(1+110.4/T_eta)/(g1+110.4/T_eta);
temp2=(1+110.4/T_eta)/(g1+110.4/T_eta)^2;
temp3=(1+110.4/T_eta)/(g1+110.4/T_eta)^3;
ddC1ds1=g4*(0.5/sqrt(g1)*temp1-sqrt(g1)*temp2)+g2*((-0.25*g1^(-3/2))*g3*temp1-0.5/sqrt(g1)*g3*temp2)-g2*(0.5/sqrt(g1)*g3*temp2-2*sqrt(g1)*g3*temp3);
dC1ds1=0.5/sqrt(g1)*g3*temp1-sqrt(g1)*g3*temp2;
A1=(-ddC1ds1*f3*C1-f6*dC1*C1-f4*f3*C1-f6*f1*C1+dC1ds1*(dC1*f3+f1*f3))/C1^2;
ddC1ds2=g6*(0.5/sqrt(g1)*temp1-sqrt(g1)*temp2)+g2*((-0.25*g1^(-3/2))*g5*temp1-0.5/sqrt(g1)*g5*temp2)-g2*(0.5/sqrt(g1)*g5*temp2-2*sqrt(g1)*g5*temp3);
dC1ds2=0.5/sqrt(g1)*g5*temp1-sqrt(g1)*g5*temp2;
A2=(-ddC1ds2*f3*C1-f9*dC1*C1-f7*f3*C1-f9*f1*C1+dC1ds2*(dC1*f3+f1*f3))/C1^2;
temp4=(1-gamma)*Ma^2*Pr;
A3=temp4/C1*(dC1ds1*f3^2+C1*2*f3*f6)-1/C1*(g2*ddC1ds1+dC1*g4+Pr*g2*f4+Pr*f1*g4)-(temp4*C1*f3^2-dC1*g2-Pr*f1*g2)/C1^2*dC1ds1;
A4=temp4/C1*(dC1ds2*f3^2+C1*2*f3*f9)-1/C1*(g2*ddC1ds2+dC1*g6+Pr*g2*f7+Pr*f1*g6)-(temp4*C1*f3^2-dC1*g2-Pr*f1*g2)/C1^2*dC1ds2;
end