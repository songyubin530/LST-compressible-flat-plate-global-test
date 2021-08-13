function [dudy_stretch,dTdy_stretch,ddudy_stretch,ddTdy_stretch]=derivation_base(nx_base,y,y_stretch,u_eta,T_eta)
hx=1/(nx_base-1);
sx=getsx(y,nx_base,hx);
sx=1./sx;
%dudy_init为基本流坐标下的导数（因为拉伸后的坐标，网格点比较少，直接用拉伸后的网格误差较大；使用初始网格求得导数后，再映射到拉伸后的网格，导数更加精确）
%dudy_stretch为拉伸后的网格
dudy_init=derivation(u_eta,nx_base,sx,1);
dTdy_init=derivation(T_eta,nx_base,sx,1);
ddudy_init=derivation(dudy_init,nx_base,sx,1);
ddTdy_init=derivation(dTdy_init,nx_base,sx,1);

dudy_stretch=interp1(y,dudy_init,y_stretch);
dTdy_stretch=interp1(y,dTdy_init,y_stretch);
ddudy_stretch=interp1(y,ddudy_init,y_stretch);
ddTdy_stretch=interp1(y,ddTdy_init,y_stretch);
end

function df=derivation(f,ny,sy,SLY)
hy=SLY/(ny-1);
df(1,1)=(-11*f(1)+18*f(2)-9*f(3)+2*f(4))/6/hy;
df(2,1)=(-2*f(1)-3*f(2)+6*f(3)-f(4))/6/hy; 
df(3,1)=(f(1)-6*f(2)+3*f(3)+2*f(4))/6/hy;
% df(3,1)=(8*(f(4)-f(2))-(f(5)-f(1)))/12/hy;
df(ny-2,1)=-(f(ny)-6*f(ny-1)+3*f(ny-2)+2*f(ny-3))/6/hy;
df(ny-1,1)=-(-2*f(ny)-3*f(ny-1)+6*f(ny-2)-f(ny-3))/6/hy; 
df(ny,1)=-(-11*f(ny)+18*f(ny-1)-9*f(ny-2)+2*f(ny-3))/6/hy;
for i=4:ny-3
df(i,1)=(45*(f(i+1)-f(i-1))-9*(f(i+2)-f(i-2))+(f(i+3)-f(i-3)))/(hy*60);
end
    df=df.*sy;
end


function fy=getsx(f,ny,hy)
    b1=8/(12*hy);
    b2=1/(12*hy);
    a1=1/(60*hy);
    a2=-3/(20*hy);
    a3=3/(4*hy);
    for j=4:ny-3
        fy(j,1)=a1*(f(j+3)-f(j-3))+a2*(f(j+2)-f(j-2))+a3*(f(j+1)-f(j-1));
    end
    fy(1,1)=(-11*f(1)+18*f(2)-9*f(3)+2*f(4))/(6*hy);
    fy(2,1)=(-2*f(1)-3*f(2)+6*f(3)-f(4))/(6*hy);  
    fy(3,1)=(f(1)-6*f(2)+3*f(3)+2*f(4))/6/hy;
%     fy(3,1)=b1*(f(4)-f(2))-b2*(f(5)-f(1));
    fy(ny-2,1)=b1*(f(ny-1)-f(ny-3))-b2*(f(ny)-f(ny-4));
    fy(ny-1,1)=(f(ny-3)-6*f(ny-2)+3*f(ny-1) +2*f(ny))/(6*hy);
    fy(ny,1)=(f(ny-2)-4*f(ny-1)+3*f(ny))/(2*hy);
end

function fy=getsx_new(yl,y_max,ny,SLY)
a=yl*y_max/(y_max-2*yl);
b=1+2*a/y_max;
dy_cal=SLY/(ny-1);
for i=1:ny
% temp1=-sin(pi*(i-1)/(ny-1));
% temp2=1+cos(pi*(i-1)/(ny-1));
% temp3=b-cos(pi*(i-1)/(ny-1));
% dy_phy(i,1)=a*(temp1*temp3*pi/(ny-1)+temp1*temp2*pi/(ny-1))/(temp3)^2;    
temp1=-sin(pi*(ny-i)/(ny-1));
temp2=1+cos(pi*(ny-i)/(ny-1));
temp3=b-cos(pi*(ny-i)/(ny-1));
dy_phy(i,1)=-a*(temp1*temp3*pi/(ny-1)+temp1*temp2*pi/(ny-1))/(temp3)^2;  
end
fy=dy_phy./dy_cal;
end