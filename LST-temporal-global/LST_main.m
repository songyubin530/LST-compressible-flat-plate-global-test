%--------------------------------------------------------------------------
%由浙江大学航空航天学院夏振华副教授课题组开发
%--------------------------------------------------------------------------
%线性稳定性计算程序-时间稳定性-全局方法，包含以下模块：
%（1）基本流计算；
%（2）网格分布、将基本流映射到计算网格；
%（3）计算基本流导数；
%（4）计算矩阵A、B、C；
%（5）计算切比雪夫多项式导数矩阵元素；
%（6）计算矩阵L、矩阵R1、矩阵R2，以及特征值矩阵A_overline，B_overline；
%（7）从内点方程中消去边界项；
%（8）求解特征值问题。

clear all
%基本参数-------------------------------------------------------------------
%基本流计算时：boundary_layer_condition=0，表示绝热壁面；boundary_layer_condition=1，表示等温壁面。
boundary_layer_condition=0;
Ma=10.0;
Pr=0.7;
Re_x=1000;
gamma=1.4;
%计算参考温度，Malik给出了驻点温度，且温度单位为R，K=R/1.8;
tem_refer_T=1+(gamma-1)/2*Ma^2;
T_reference=4200/1.8/tem_refer_T;
%等温壁边界条件还需要给定壁温
% Tw=(1+sqrt(Pr)*(gamma-1)/2*Ma^2)*0.1;
%基本流网格数
nx_base=10001;
%基本流计算域
SLX_base=120;
%流向波数
alpha=0.12;
%展向波数
beta=0.0;
%谱方法网格点数
N=300;
ny=N+1;
%拉伸参数
yi=32;
%稳定性计算域
y_max=100;
%（1）基本流计算；----------------------------------------------------------
if boundary_layer_condition==0 
[y,d_eta,u_eta,v_eta,T_eta,Amu_eta]=shooting_method_adiabat(Ma,Pr,gamma,T_reference,nx_base,SLX_base,Re_x);
else
[y,d_eta,u_eta,v_eta,T_eta,Amu_eta]=shooting_method_isothermy(Ma,Pr,gamma,T_reference,nx_base,SLX_base,Re_x,Tw);   
end
p_eta=d_eta.*T_eta/gamma/Ma^2;
%（2）网格分布、将基本流映射到计算网格；--------------------------------------
[xi,dxidy,y_stretch]=coordinate_stretch(ny,y_max,yi);
u_stretch=interp1(y,u_eta,y_stretch);
v_stretch=interp1(y,v_eta,y_stretch);
p_stretch=interp1(y,p_eta,y_stretch);
T_stretch=interp1(y,T_eta,y_stretch);
Amu_stretch=interp1(y,Amu_eta,y_stretch);
%（3）计算基本流导数；------------------------------------------------------
[dudy_stretch,dTdy_stretch,ddudy_stretch,ddTdy_stretch]=derivation_base(nx_base,y,y_stretch,u_eta,T_eta);
[dAmudT_stretch,dkdT_stretch,dlambdadT_stretch,ddAmudT_stretch,ddkdT_stretch]=dfdT_base(y,y_stretch,T_reference,T_eta);
%（4）计算矩阵A、B、C；-----------------------------------------------------
%矩阵的指标与坐标点对应关系为：j=1：ny——0：N，但是0表示y_max;N表示壁面，即y=0；
[A,B,C,C_omega]=getABC(ny,alpha,beta,Re_x,gamma,Ma,Pr,...
    u_stretch,T_stretch,Amu_stretch,dudy_stretch,dTdy_stretch,...
    ddudy_stretch,ddTdy_stretch,dAmudT_stretch,dkdT_stretch,dlambdadT_stretch,ddAmudT_stretch,ddkdT_stretch);
%（5）计算切比雪夫多项式导数矩阵元素；----------------------------------------
[F,G]=get_Chebyshev_derivative(N,ny,xi,dxidy);
%（6）计算矩阵L、矩阵R1、矩阵R2，以及特征值矩阵A_overline，B_overline；-------
[R,R_omega,L]=get_RL(ny,N,F,G,A,B,C,C_omega);
A_overline=L+R;B_overline=-R_omega;
%（7）从内点方程中消去边界项；-----------------------------------------------
[A1,B1,phi]=get_Elimination(N,A_overline,B_overline);
%（8）求解特征值问题。-------------------------------------------------------
[vector,omega]=eig(phi);
omega=diag(omega);
%注：计算时会存在很多虚假的不稳定特征值，这些虚假的特征值会随网格数的变化而变化，所以改变网格数可以很好地区分虚假不稳定特征值和真实不稳定特征值；
%计算扰动特征函数-----------------------------------------------------------
[omega_row,omega_col]=find(imag(omega)>0);
eigenfunction=vector(:,omega_row(end));
eigenfunction_U=eigenfunction(1:5:end);
eigenfunction_V=eigenfunction(2:5:end);
eigenfunction_p=eigenfunction(3:5:end);
eigenfunction_Temperature=eigenfunction(4:5:end);
eigenfunction_W=eigenfunction(5:5:end);
%用温度扰动函数的最大值对其进行归一化-----------------------------------------
eigenfunction_Temperature_real=real(eigenfunction_Temperature)/max(abs(eigenfunction_Temperature));
eigenfunction_Temperature_imag=imag(eigenfunction_Temperature)/max(abs(eigenfunction_Temperature));
%数据对比-Malik温度特征函数沿法向分布----------------------------------------
fid1=('T_real_Malik.dat');
[yy_real_Malik,T_real_Malik]=textread(fid1,'%f%f',1000,'delimiter', ',','headerlines',3);
fid1=('T_imag_Malik.dat');
[yy_imag_Malik,T_imag_Malik]=textread(fid1,'%f%f',1000,'delimiter', ',','headerlines',3);
%对比结果-------------------------------------------------------------------
% 创建 figure
figure1 = figure('Name','compare');
% 创建 axes
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.60 0.775 0.33]);
hold(axes1,'on');

% 创建 plot
plot(real(omega),imag(omega),'Parent',axes1,'MarkerSize',8,'Marker','*','LineStyle','none',...
    'Color',[0 0 0]);
plot(real(omega(omega_row(end))),imag(omega(omega_row(end))),'ro');
text(real(omega(omega_row(end))),imag(omega(omega_row(end))),sprintf('(%.8f,%.8f)',real(omega(omega_row(2))),imag(omega(omega_row(2)))),...
    'VerticalAlignment','bottom','FontSize',14)

% 创建 ylabel
ylabel({'$\omega_i$'},'Interpreter','latex');

% 创建 xlabel
xlabel({'$\omega_r$'},'Interpreter','latex');

% 创建 title
title({'特征值分布'});

% 取消以下行的注释以保留坐标区的 X 范围
xlim(axes1,[0 0.2]);
% 取消以下行的注释以保留坐标区的 Y 范围
ylim(axes1,[-0.05 0.05]);
box(axes1,'on');
% 设置其余坐标区属性
set(axes1,'FontSize',30,'LineWidth',3);
% 创建 subplot
subplot1 = subplot(2,1,2,'Parent',figure1);
hold(subplot1,'on');

% 使用 plot 的矩阵输入创建多行
plot1 = plot([y_stretch(2:end-1),y_stretch(2:end-1)],[eigenfunction_Temperature_real,eigenfunction_Temperature_imag],'Parent',subplot1,'LineWidth',2);
set(plot1(1),'DisplayName','Temperature-real','Color',[1 0 0]);
set(plot1(2),'DisplayName','Temperature-imag','Color',[0 0 1]);

% 创建 plot
plot(yy_real_Malik,T_real_Malik,'Parent',subplot1,'DisplayName','Temperature-real-Malik',...
    'MarkerSize',8,...
    'Marker','*',...
    'LineStyle','none',...
    'Color',[1 0 0]);

% 创建 plot
plot(yy_imag_Malik,T_imag_Malik,'Parent',subplot1,'DisplayName','Temperature-imag-Malik',...
    'MarkerSize',8,...
    'Marker','*',...
    'LineStyle','none',...
    'Color',[0 0 1]);

% 创建 ylabel
ylabel({'Temperature'},'Interpreter','latex');

% 创建 xlabel
xlabel({'y'},'Interpreter','latex');

% 创建 title
title({'温度特征函数沿法向分布'});

% 取消以下行的注释以保留坐标区的 X 范围
xlim(subplot1,[0 50]);
% 取消以下行的注释以保留坐标区的 Y 范围
ylim(subplot1,[-0.6 1.2]);
box(subplot1,'on');
% 设置其余坐标区属性
set(subplot1,'FontSize',30,'LineWidth',3);
% 创建 legend
legend1 = legend(subplot1,'show');
set(legend1,...
    'Position',[0.736770460312534 0.26958872861114 0.147916663050031 0.124451750856743],...
    'FontSize',16);
set(figure1, 'position', get(0,'ScreenSize'));