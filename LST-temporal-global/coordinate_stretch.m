function [xi_com,dxidy,y_stretch]=coordinate_stretch(ny,y_max,yi)
xi_com=zeros(ny,1);
dxidy=zeros(ny,1);
y_stretch=zeros(ny,1);
a=yi*y_max/(y_max-2*yi);
b=1+2*a/y_max;
for i=1:ny
xi_com(i,1)=cos(pi*(i-1)/(ny-1));
y_stretch(i,1)=a*(1+xi_com(i))/(b-xi_com(i));
dxidy(i,1)=(b/a*(y_stretch(i,1)/a+1)-1/a*(b*y_stretch(i,1)/a-1))/(y_stretch(i,1)/a+1)^2;
end
% y_stretch=flipud(y_stretch);
end