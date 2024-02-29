function [flow]=NPSE_Baseflow(MESH,NPSE)

global parameter
 Pr=parameter.Pr;              %普朗特数
 r=parameter. r;                  %气体常数
 Ma=parameter. Ma;                 %马赫数  
 Te=parameter. Te;              
 Re0=parameter.Re0;
 
 ym=MESH.y;
 Nx=MESH.Nx;
 Ny=MESH.Ny;
 a1=110.4/Te;                
 a2=194/Te;                  
 
 X0=Re0;
 Xmax=NPSE.Xmax;
 Xall=NPSE.X;
 ReAll=sqrt(Re0*Xall);
 
%存储基本流用的结构体
flow.U0=zeros(Ny,Nx);
flow.Ux0=zeros(Ny,Nx);
flow.Uy0=zeros(Ny,Nx);
flow.Uxy0=zeros(Ny,Nx);
flow.Uyy0=zeros(Ny,Nx);

flow.V0=zeros(Ny,Nx);
flow.Vx0=zeros(Ny,Nx);
flow.Vy0=zeros(Ny,Nx);
flow.Vxy0=zeros(Ny,Nx);
flow.Vyy0=zeros(Ny,Nx);

flow.T0=zeros(Ny,Nx);
flow.Tx0=zeros(Ny,Nx);
flow.Ty0=zeros(Ny,Nx);
flow.Tyy0=zeros(Ny,Nx);

flow.Den0=zeros(Ny,Nx);
flow.Denx0=zeros(Ny,Nx);
flow.Deny0=zeros(Ny,Nx);
flow.Denyy0=zeros(Ny,Nx);

 
 
solinit = bvpinit(linspace(0,10,100),@NPSE_guess);
sol=bvp4c(@NPSE_bvpeq,@NPSE_bvpbc,solinit);
t0=linspace(0,10,20000);                    
y=deval(sol,t0);                                 


f=y(1,:);
f1=y(2,:);   %f一阶导
f2=y(3,:);   %f二阶导
g=y(4,:);
g1=y(5,:);   %g一阶导
G=y(6,:);    %g积分，即η

%U0=f1;
%T0=g;
%plot(G,U0,'c')
%plot(G,T0,'m')
%hold on

C1=sqrt(g)*(1+a1)./(g+a1);                      %C1
C2=sqrt(g)*(1+a2)./(g+a2);                      %C2
C11=(1+a1)*g1.*(0.5*g.^(-0.5).*(g+a1)-sqrt(g))./((g+a1).^2); %C1的一阶导
f3=(-0.5*f.*f2-C11.*f2)./C1;                                                 %f的三阶导
g2=(-Pr*(r-1)*Ma^2*C1.*f2.^2 -0.5*Pr*f.*g1-C11.*g1)./C2;       %g的两阶导


for i=1:Nx
   R=ReAll(i);
   Y=G;
  %速度U各分量场      
  U=f1;
  Ux=-0.5/R*f2.*G./g;
  Uy=f2./g;
  Uxy=0.5/R*(f2.*g1.*G./g.^3-f2./g-f3.*G./g.^2);
  Uyy=(f3.*g-f2.*g1)./g.^3;
  
  %速度V各分量场
  V=0.5/R*(G.*f1-f.*g);
  Vx=0.25/R^2*(f.*g1.*G./g-G.^2.*f2./g+f.*g-G.*f1);
  Vy=0.5/R*(G.*f2-f.*g1)./g;
  Vxy=0.25/R^2*((f1.*g1.*G+f.*g2.*G+f.*g1.*g-2*f2.*G.*g-G.^2.*f3)./g+(f2.*g1.*G.^2-f.*g1.^2.*G)./g.^2+f.*g1-G.*f2)./g;
  Vyy=0.5/R*(f2.*g.^2+f3.*g.*G-f1.*g1.*g-f.*g2.*g-f2.*g1.*G+f.*g1.^2);
  
  %温度T各分量场
  T=g;
  Tx=-0.5/R*(g1.*G./g);
  Ty=g1./g;
  Tyy=(g2.*g-g1.^2)./g.^3;
  
  %密度场ρ各分量
  Den=1./g;
  Denx=0.5/R*(G.*g1./g.^3);
  Deny=-g1./g.^3;
  Denyy=(3*g1.^2.*g.^2-g2.*g.^3)./g.^7;
  
%对各物理场进行扩展
[m,n]=size(Y);
ex=10000;
Y_e=linspace(30,300,ex+1);
Y0=[Y Y_e(2:end)]*R/Re0;

U_e=U(n)*ones(1,ex);
Ux_e=Ux(n)*ones(1,ex);
Uy_e=Uy(n)*ones(1,ex);
Uxy_e=Uxy(n)*ones(1,ex);
Uyy_e=Uyy(n)*ones(1,ex);

V_e=V(n)*ones(1,ex);
Vx_e=Vx(n)*ones(1,ex);
Vy_e=Vy(n)*ones(1,ex);
Vxy_e=Vxy(n)*ones(1,ex);
Vyy_e=Vyy(n)*ones(1,ex);

T_e=T(n)*ones(1,ex);
Tx_e=Tx(n)*ones(1,ex);
Ty_e=Ty(n)*ones(1,ex);
Tyy_e=Tyy(n)*ones(1,ex);

Den_e=Den(n)*ones(1,ex);
Denx_e=Denx(n)*ones(1,ex);
Deny_e=Deny(n)*ones(1,ex);
Denyy_e=Denyy(n)*ones(1,ex);

U=[U U_e];
Ux=[Ux Ux_e]*Re0/R;
Uy=[Uy Uy_e]*Re0/R;
Uxy=[Uxy Uxy_e]*(Re0/R)^2;
Uyy=[Uyy Uyy_e]*(Re0/R)^2;

V=[V V_e];
Vx=[Vx Vx_e]*Re0/R;
Vy=[Vy Vy_e]*Re0/R;
Vxy=[Vxy Vxy_e]*(Re0/R)^2;
Vyy=[Vyy Vyy_e]*(Re0/R)^2;

T=[T T_e];
Tx=[Tx Tx_e]*Re0/R;
Ty=[Ty Ty_e]*Re0/R;
Tyy=[Tyy Tyy_e]*(Re0/R)^2;

Den=[Den Den_e];
Denx=[Denx Denx_e]*Re0/R;
Deny=[Deny Deny_e]*Re0/R;
Denyy=[Denyy Denyy_e]*(Re0/R)^2;

%插值
flow.U0(:,i)=interp1(Y0,U,ym,'spline');
flow.Ux0(:,i)=interp1(Y0,Ux,ym,'spline');
flow.Uy0(:,i)=interp1(Y0,Uy,ym,'spline');
flow.Uxy0(:,i)=interp1(Y0,Uxy,ym,'spline');
flow.Uyy0(:,i)=interp1(Y0,Uyy,ym,'spline');

flow.V0(:,i)=interp1(Y0,V,ym,'spline');
flow.Vx0(:,i)=interp1(Y0,Vx,ym,'spline');
flow.Vy0(:,i)=interp1(Y0,Vy,ym,'spline');
flow.Vxy0(:,i)=interp1(Y0,Vxy,ym,'spline');
flow.Vyy0(:,i)=interp1(Y0,Vyy,ym,'spline');

flow.T0(:,i)=interp1(Y0,T,ym,'spline');
flow.Tx0(:,i)=interp1(Y0,Tx,ym,'spline');
flow.Ty0(:,i)=interp1(Y0,Ty,ym,'spline');
flow.Tyy0(:,i)=interp1(Y0,Tyy,ym,'spline');

flow.Den0(:,i)=interp1(Y0,Den,ym,'spline');
flow.Denx0(:,i)=interp1(Y0,Denx,ym,'spline');
flow.Deny0(:,i)=interp1(Y0,Deny,ym,'spline');
flow.Denyy0(:,i)=interp1(Y0,Denyy,ym,'spline');

end
flow.R=ReAll;
end
