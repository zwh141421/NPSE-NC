function [flow0]=NPSE_Baseflow(y,Nx,Ny,X,NPSE)

parameter=NPSE_SetupParameter;
 Pr  = parameter(1);              
 r   = parameter(2);                
 Ma  = parameter(3);                  
 Te  = parameter(4);              
 Re0 = parameter(7);
 
 ym = y;
 a1 = 110.4/Te;                
 a2 = 194/Te;                  
 
 X0   = NPSE(6);
 Xmax = NPSE(7);
 Xall=X;
 ReAll=sqrt(Re0*Xall);
 

flow0.U0=zeros(Ny,Nx);
flow0.Ux0=zeros(Ny,Nx);
flow0.Uy0=zeros(Ny,Nx);
flow0.Uxy0=zeros(Ny,Nx);
flow0.Uyy0=zeros(Ny,Nx);

flow0.V0=zeros(Ny,Nx);
flow0.Vx0=zeros(Ny,Nx);
flow0.Vy0=zeros(Ny,Nx);
flow0.Vxy0=zeros(Ny,Nx);
flow0.Vyy0=zeros(Ny,Nx);

flow0.T0=zeros(Ny,Nx);
flow0.Tx0=zeros(Ny,Nx);
flow0.Ty0=zeros(Ny,Nx);
flow0.Tyy0=zeros(Ny,Nx);

flow0.Den0=zeros(Ny,Nx);
flow0.Denx0=zeros(Ny,Nx);
flow0.Deny0=zeros(Ny,Nx);
flow0.Denyy0=zeros(Ny,Nx);

 
 
solinit = bvpinit(linspace(0,10,100),@NPSE_guess);
sol=bvp4c(@NPSE_bvpeq,@NPSE_bvpbc,solinit);
t0=linspace(0,10,20000);                    
y=deval(sol,t0);                                 


f=y(1,:);
f1=y(2,:);   %first derivative of f
f2=y(3,:);   %second derivative of f
g=y(4,:);
g1=y(5,:);   %first derivative of g
G=y(6,:);    %Integral of g, η

%U0=f1;
%T0=g;
%plot(G,U0,'c')
%plot(G,T0,'m')
%hold on

C1=sqrt(g)*(1+a1)./(g+a1);                      %C1
C2=sqrt(g)*(1+a2)./(g+a2);                      %C2
C11=(1+a1)*g1.*(0.5*g.^(-0.5).*(g+a1)-sqrt(g))./((g+a1).^2); %first derivative of C1
f3=(-0.5*f.*f2-C11.*f2)./C1;                                                 %third derivative of f
g2=(-Pr*(r-1)*Ma^2*C1.*f2.^2 -0.5*Pr*f.*g1-C11.*g1)./C2;       %second derivative of g


for i=1:Nx
   R=ReAll(i);
   Y=G;
  %The component of the velocity U      
  U=f1;
  Ux=-0.5/R*f2.*G./g;
  Uy=f2./g;
  Uxy=0.5/R*(f2.*g1.*G./g.^3-f2./g-f3.*G./g.^2);
  Uyy=(f3.*g-f2.*g1)./g.^3;
  
  %The component of the velocity V
  V=0.5/R*(G.*f1-f.*g);
  Vx=0.25/R^2*(f.*g1.*G./g-G.^2.*f2./g+f.*g-G.*f1);
  Vy=0.5/R*(G.*f2-f.*g1)./g;
  Vxy=0.25/R^2*((f1.*g1.*G+f.*g2.*G+f.*g1.*g-2*f2.*G.*g-G.^2.*f3)./g+(f2.*g1.*G.^2-f.*g1.^2.*G)./g.^2+f.*g1-G.*f2)./g;
  Vyy=0.5/R*(f2.*g.^2+f3.*g.*G-f1.*g1.*g-f.*g2.*g-f2.*g1.*G+f.*g1.^2);
  
  %The component of the temperature T
  T=g;
  Tx=-0.5/R*(g1.*G./g);
  Ty=g1./g;
  Tyy=(g2.*g-g1.^2)./g.^3;
  
  %The component of the density ρ
  Den=1./g;
  Denx=0.5/R*(G.*g1./g.^3);
  Deny=-g1./g.^3;
  Denyy=(3*g1.^2.*g.^2-g2.*g.^3)./g.^7;
  
%Expand the physical fields
[~,n]=size(Y);
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

%interpolation
flow0.U0(:,i)=interp1(Y0,U,ym,'spline');
flow0.Ux0(:,i)=interp1(Y0,Ux,ym,'spline');
flow0.Uy0(:,i)=interp1(Y0,Uy,ym,'spline');
flow0.Uxy0(:,i)=interp1(Y0,Uxy,ym,'spline');
flow0.Uyy0(:,i)=interp1(Y0,Uyy,ym,'spline');

flow0.V0(:,i)=interp1(Y0,V,ym,'spline');
flow0.Vx0(:,i)=interp1(Y0,Vx,ym,'spline');
flow0.Vy0(:,i)=interp1(Y0,Vy,ym,'spline');
flow0.Vxy0(:,i)=interp1(Y0,Vxy,ym,'spline');
flow0.Vyy0(:,i)=interp1(Y0,Vyy,ym,'spline');

flow0.T0(:,i)=interp1(Y0,T,ym,'spline');
flow0.Tx0(:,i)=interp1(Y0,Tx,ym,'spline');
flow0.Ty0(:,i)=interp1(Y0,Ty,ym,'spline');
flow0.Tyy0(:,i)=interp1(Y0,Tyy,ym,'spline');

flow0.Den0(:,i)=interp1(Y0,Den,ym,'spline');
flow0.Denx0(:,i)=interp1(Y0,Denx,ym,'spline');
flow0.Deny0(:,i)=interp1(Y0,Deny,ym,'spline');
flow0.Denyy0(:,i)=interp1(Y0,Denyy,ym,'spline');

end
flow0.R=ReAll;
end
