function [A2,B2,phi]=get_Elimination(N,A_overline,B_overline)
A1=zeros(5*N-5,5*N-3);
B1=zeros(5*N-5,5*N-3);
% A2=zeros(5*N-5,5*N-5);
% B2=zeros(5*N-5,5*N-5);
% phi=zeros(5*N-5,5*N-5);
for m=1:5*N-5
tem1=(A_overline(m,5*N-4)*A_overline(5*N-3,5*N-3)/A_overline(5*N-3,5*N-4)-A_overline(m,5*N-3))/(A_overline(5*N-4,5*N-4)*A_overline(5*N-3,5*N-3)/A_overline(5*N-3,5*N-4)-A_overline(5*N-4,5*N-3));
tem2=(A_overline(m,5*N-3)*A_overline(5*N-4,5*N-4)/A_overline(5*N-4,5*N-3)-A_overline(m,5*N-4))/(A_overline(5*N-3,5*N-3)*A_overline(5*N-4,5*N-4)/A_overline(5*N-4,5*N-3)-A_overline(5*N-3,5*N-4));
A1(m,:)=A_overline(m,:)-(tem1*A_overline(5*N-4,:)+tem2*A_overline(5*N-3,:));
B1(m,:)=B_overline(m,:)-(tem1*B_overline(5*N-4,:)+tem2*B_overline(5*N-3,:));
end
A2=A1(1:5*N-5,1:5*N-5);
B2=B1(1:5*N-5,1:5*N-5);
phi=B2^(-1)*A2;
end