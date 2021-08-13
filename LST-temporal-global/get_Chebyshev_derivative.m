function [F,G]=get_Chebyshev_derivative(N,ny,xi,dxidy)
ck=ones(ny,1);
E=zeros(ny,ny);
ck(1)=2;ck(ny)=2;
for j=1:ny
    for k=1:ny
        if (j==k)
            if (j==1)
                E(j,k)=(2*N^2+1)/6;
            elseif (j==ny)
                E(j,k)=-(2*N^2+1)/6;
            else
                E(j,k)=-xi(j)/2/(1-xi(j)^2);
            end              
        else
            E(j,k)=ck(j)/ck(k)*(-1)^(j+k)/(xi(j)-xi(k));
        end
    end
end
F=zeros(ny,ny);
for j=1:ny
    for k=1:ny
    F(j,k)=dxidy(j)*E(j,k);
    end
end
G=F*F;
end