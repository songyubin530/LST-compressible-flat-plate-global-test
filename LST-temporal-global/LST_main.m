%--------------------------------------------------------------------------
%���㽭��ѧ���պ���ѧԺ���񻪸����ڿ����鿪��
%--------------------------------------------------------------------------
%�����ȶ��Լ������-ʱ���ȶ���-ȫ�ַ�������������ģ�飺
%��1�����������㣻
%��2������ֲ�����������ӳ�䵽��������
%��3�����������������
%��4���������A��B��C��
%��5�������б�ѩ�����ʽ��������Ԫ�أ�
%��6���������L������R1������R2���Լ�����ֵ����A_overline��B_overline��
%��7�����ڵ㷽������ȥ�߽��
%��8���������ֵ���⡣

clear all
%��������-------------------------------------------------------------------
%����������ʱ��boundary_layer_condition=0����ʾ���ȱ��棻boundary_layer_condition=1����ʾ���±��档
boundary_layer_condition=0;
Ma=10.0;
Pr=0.7;
Re_x=1000;
gamma=1.4;
%����ο��¶ȣ�Malik������פ���¶ȣ����¶ȵ�λΪR��K=R/1.8;
tem_refer_T=1+(gamma-1)/2*Ma^2;
T_reference=4200/1.8/tem_refer_T;
%���±ڱ߽���������Ҫ��������
% Tw=(1+sqrt(Pr)*(gamma-1)/2*Ma^2)*0.1;
%������������
nx_base=10001;
%������������
SLX_base=120;
%������
alpha=0.12;
%չ����
beta=0.0;
%�׷����������
N=300;
ny=N+1;
%�������
yi=32;
%�ȶ��Լ�����
y_max=100;
%��1�����������㣻----------------------------------------------------------
if boundary_layer_condition==0 
[y,d_eta,u_eta,v_eta,T_eta,Amu_eta]=shooting_method_adiabat(Ma,Pr,gamma,T_reference,nx_base,SLX_base,Re_x);
else
[y,d_eta,u_eta,v_eta,T_eta,Amu_eta]=shooting_method_isothermy(Ma,Pr,gamma,T_reference,nx_base,SLX_base,Re_x,Tw);   
end
p_eta=d_eta.*T_eta/gamma/Ma^2;
%��2������ֲ�����������ӳ�䵽��������--------------------------------------
[xi,dxidy,y_stretch]=coordinate_stretch(ny,y_max,yi);
u_stretch=interp1(y,u_eta,y_stretch);
v_stretch=interp1(y,v_eta,y_stretch);
p_stretch=interp1(y,p_eta,y_stretch);
T_stretch=interp1(y,T_eta,y_stretch);
Amu_stretch=interp1(y,Amu_eta,y_stretch);
%��3�����������������------------------------------------------------------
[dudy_stretch,dTdy_stretch,ddudy_stretch,ddTdy_stretch]=derivation_base(nx_base,y,y_stretch,u_eta,T_eta);
[dAmudT_stretch,dkdT_stretch,dlambdadT_stretch,ddAmudT_stretch,ddkdT_stretch]=dfdT_base(y,y_stretch,T_reference,T_eta);
%��4���������A��B��C��-----------------------------------------------------
%�����ָ����������Ӧ��ϵΪ��j=1��ny����0��N������0��ʾy_max;N��ʾ���棬��y=0��
[A,B,C,C_omega]=getABC(ny,alpha,beta,Re_x,gamma,Ma,Pr,...
    u_stretch,T_stretch,Amu_stretch,dudy_stretch,dTdy_stretch,...
    ddudy_stretch,ddTdy_stretch,dAmudT_stretch,dkdT_stretch,dlambdadT_stretch,ddAmudT_stretch,ddkdT_stretch);
%��5�������б�ѩ�����ʽ��������Ԫ�أ�----------------------------------------
[F,G]=get_Chebyshev_derivative(N,ny,xi,dxidy);
%��6���������L������R1������R2���Լ�����ֵ����A_overline��B_overline��-------
[R,R_omega,L]=get_RL(ny,N,F,G,A,B,C,C_omega);
A_overline=L+R;B_overline=-R_omega;
%��7�����ڵ㷽������ȥ�߽��-----------------------------------------------
[A1,B1,phi]=get_Elimination(N,A_overline,B_overline);
%��8���������ֵ���⡣-------------------------------------------------------
[vector,omega]=eig(phi);
omega=diag(omega);
%ע������ʱ����ںܶ���ٵĲ��ȶ�����ֵ����Щ��ٵ�����ֵ�����������ı仯���仯�����Ըı����������Ժܺõ�������ٲ��ȶ�����ֵ����ʵ���ȶ�����ֵ��
%�����Ŷ���������-----------------------------------------------------------
[omega_row,omega_col]=find(imag(omega)>0);
eigenfunction=vector(:,omega_row(end));
eigenfunction_U=eigenfunction(1:5:end);
eigenfunction_V=eigenfunction(2:5:end);
eigenfunction_p=eigenfunction(3:5:end);
eigenfunction_Temperature=eigenfunction(4:5:end);
eigenfunction_W=eigenfunction(5:5:end);
%���¶��Ŷ����������ֵ������й�һ��-----------------------------------------
eigenfunction_Temperature_real=real(eigenfunction_Temperature)/max(abs(eigenfunction_Temperature));
eigenfunction_Temperature_imag=imag(eigenfunction_Temperature)/max(abs(eigenfunction_Temperature));
%���ݶԱ�-Malik�¶����������ط���ֲ�----------------------------------------
fid1=('T_real_Malik.dat');
[yy_real_Malik,T_real_Malik]=textread(fid1,'%f%f',1000,'delimiter', ',','headerlines',3);
fid1=('T_imag_Malik.dat');
[yy_imag_Malik,T_imag_Malik]=textread(fid1,'%f%f',1000,'delimiter', ',','headerlines',3);
%�ԱȽ��-------------------------------------------------------------------
% ���� figure
figure1 = figure('Name','compare');
% ���� axes
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.60 0.775 0.33]);
hold(axes1,'on');

% ���� plot
plot(real(omega),imag(omega),'Parent',axes1,'MarkerSize',8,'Marker','*','LineStyle','none',...
    'Color',[0 0 0]);
plot(real(omega(omega_row(end))),imag(omega(omega_row(end))),'ro');
text(real(omega(omega_row(end))),imag(omega(omega_row(end))),sprintf('(%.8f,%.8f)',real(omega(omega_row(2))),imag(omega(omega_row(2)))),...
    'VerticalAlignment','bottom','FontSize',14)

% ���� ylabel
ylabel({'$\omega_i$'},'Interpreter','latex');

% ���� xlabel
xlabel({'$\omega_r$'},'Interpreter','latex');

% ���� title
title({'����ֵ�ֲ�'});

% ȡ�������е�ע���Ա����������� X ��Χ
xlim(axes1,[0 0.2]);
% ȡ�������е�ע���Ա����������� Y ��Χ
ylim(axes1,[-0.05 0.05]);
box(axes1,'on');
% ������������������
set(axes1,'FontSize',30,'LineWidth',3);
% ���� subplot
subplot1 = subplot(2,1,2,'Parent',figure1);
hold(subplot1,'on');

% ʹ�� plot �ľ������봴������
plot1 = plot([y_stretch(2:end-1),y_stretch(2:end-1)],[eigenfunction_Temperature_real,eigenfunction_Temperature_imag],'Parent',subplot1,'LineWidth',2);
set(plot1(1),'DisplayName','Temperature-real','Color',[1 0 0]);
set(plot1(2),'DisplayName','Temperature-imag','Color',[0 0 1]);

% ���� plot
plot(yy_real_Malik,T_real_Malik,'Parent',subplot1,'DisplayName','Temperature-real-Malik',...
    'MarkerSize',8,...
    'Marker','*',...
    'LineStyle','none',...
    'Color',[1 0 0]);

% ���� plot
plot(yy_imag_Malik,T_imag_Malik,'Parent',subplot1,'DisplayName','Temperature-imag-Malik',...
    'MarkerSize',8,...
    'Marker','*',...
    'LineStyle','none',...
    'Color',[0 0 1]);

% ���� ylabel
ylabel({'Temperature'},'Interpreter','latex');

% ���� xlabel
xlabel({'y'},'Interpreter','latex');

% ���� title
title({'�¶����������ط���ֲ�'});

% ȡ�������е�ע���Ա����������� X ��Χ
xlim(subplot1,[0 50]);
% ȡ�������е�ע���Ա����������� Y ��Χ
ylim(subplot1,[-0.6 1.2]);
box(subplot1,'on');
% ������������������
set(subplot1,'FontSize',30,'LineWidth',3);
% ���� legend
legend1 = legend(subplot1,'show');
set(legend1,...
    'Position',[0.736770460312534 0.26958872861114 0.147916663050031 0.124451750856743],...
    'FontSize',16);
set(figure1, 'position', get(0,'ScreenSize'));