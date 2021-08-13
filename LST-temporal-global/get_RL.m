function [R,R_omega,L]=get_RL(ny,N,F,G,A,B,C,C_omega)
R=zeros(5*N-3,5*N-3);
for j=2:ny-1
    for k=1:5
        for l=1:5
        m=5*(j-1-1)+l;
        n=5*(j-1-1)+k;
%推导时j=0-N，我们使用为1-(N-1)；而程序中，j=1,N+1对应0-N，所以程序中需要需要减一
        R(m,n)=C(l,k,j);
        end
    end
end
R_omega=zeros(5*N-3,5*N-3);
for j=2:ny-1
    for k=1:5
        for l=1:5
        m=5*(j-1-1)+l;
        n=5*(j-1-1)+k;
        R_omega(m,n)=C_omega(l,k,j);
        end
    end
end
L=zeros(5*N-3,5*N-3);
for j=2:ny-1
    m=1+5*(j-1-1);
    for k=2:ny-1
        L(m,5*(k-1-1)+1)=A(1,1,j)*G(j,k)+B(1,1,j)*F(j,k);
        L(m,5*(k-1-1)+2)=B(1,2,j)*F(j,k);
        L(m,5*(k-1-1)+4)=B(1,4,j)*F(j,k);
    end
end
for j=2:ny-1
    m=2+5*(j-1-1);
    for k=2:ny-1
        L(m,5*(k-1-1)+1)=B(2,1,j)*F(j,k);
        L(m,5*(k-1-1)+2)=A(2,2,j)*G(j,k)+B(2,2,j)*F(j,k);
        L(m,5*(k-1-1)+3)=B(2,3,j)*F(j,k);
        L(m,5*(k-1-1)+5)=B(2,5,j)*F(j,k);
    end
        L(m,5*N-4)=B(2,3,j)*F(j,1);
        L(m,5*N-3)=B(2,3,j)*F(j,ny);
end
for j=2:ny-1
    m=5+5*(j-1-1);
    for k=2:ny-1
        L(m,5*(k-1-1)+2)=B(5,2,j)*F(j,k);
        L(m,5*(k-1-1)+4)=B(5,4,j)*F(j,k);
        L(m,5*(k-1-1)+5)=A(5,5,j)*G(j,k)+B(5,5,j)*F(j,k);
    end
end
for j=2:ny-1
    m=4+5*(j-1-1);
    for k=2:ny-1
        L(m,5*(k-1-1)+1)=B(4,1,j)*F(j,k);
        L(m,5*(k-1-1)+4)=A(4,4,j)*G(j,k)+B(4,4,j)*F(j,k);
        L(m,5*(k-1-1)+5)=B(4,5,j)*F(j,k);
    end
end
for j=2:ny-1
    m=3+5*(j-1-1);
    for k=2:ny-1
        L(m,5*(k-1-1)+2)=B(3,2,j)*F(j,k);
    end
end
    for k=2:ny-1
        L(5*N-4,5*(k-1-1)+1)=B(2,1,1)*F(1,k);
        L(5*N-4,5*(k-1-1)+2)=A(2,2,1)*G(1,k)+B(2,2,1)*F(1,k);
        L(5*N-4,5*(k-1-1)+3)=B(2,3,1)*F(1,k);
        L(5*N-4,5*(k-1-1)+5)=B(2,5,1)*F(1,k);
    end
        L(5*N-4,5*N-4)=B(2,3,1)*F(1,1);
        L(5*N-4,5*N-3)=B(2,3,1)*F(1,ny);
    for k=2:ny-1
        L(5*N-3,5*(k-1-1)+1)=B(2,1,ny)*F(ny,k);
        L(5*N-3,5*(k-1-1)+2)=A(2,2,ny)*G(ny,k)+B(2,2,ny)*F(ny,k);
        L(5*N-3,5*(k-1-1)+3)=B(2,3,ny)*F(ny,k);
        L(5*N-3,5*(k-1-1)+5)=B(2,5,ny)*F(ny,k);
    end
        L(5*N-3,5*N-4)=B(2,3,ny)*F(ny,1);
        L(5*N-3,5*N-3)=B(2,3,ny)*F(ny,ny);
end