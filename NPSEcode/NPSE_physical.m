function[nonlin]=NPSE_physical(flow,F,Ft,Fx,Fy,Fz,Fxx,Fyy,Fzz,Fxy,Fxz,Fyz,MESH,xi)
parameter=NPSE_SetupParameter;
r=parameter. r;
Rg=parameter.Rg;
RE=parameter.Re0;
%K=parameter.K;
%h1=1+K*y(i);
%h11=K
h1=1;
h11=0;
PR=parameter.Pr;
EC=parameter.Ec;

%基本流部分
Ny=MESH.Ny;
U=flow.U0(:,xi);   
Ux=flow.Ux0(:,xi);
%Ux=zeros(Ny,1);
Uy=flow.Uy0(:,xi);
%Uxy=flow.Uxy0(:,xi);
%Uxy=zeros(Ny,1);
%Uyy=flow.Uyy0(:,xi);
%Uxx=zeros(Ny,1);

V=flow.V0(:,xi);
%V=zeros(Ny,1);
Vx=flow.Vx0(:,xi);
%Vx=zeros(Ny,1);
Vy=flow.Vy0(:,xi);
%Vy=zeros(Ny,1);
%Vxy=flow.Vxy0(:,xi);
%Vxy=zeros(Ny,1);
%Vyy=flow.Vyy0(:,xi);
%Vyy=zeros(Ny,1);
%Vxx=zeros(Ny,1);

T=flow.T0(:,xi);
Tx=flow.Tx0(:,xi);
%Tx=zeros(Ny,1);
Ty=flow.Ty0(:,xi);
%Tyy=flow.Tyy0(:,xi);
%Txx=zeros(Ny,1);
%Txy=zeros(l,1);

Den=flow.Den0(:,xi);
Denx=flow.Denx0(:,xi);
%Denx=zeros(Ny,1);
Deny=flow.Deny0(:,xi);
%Denyy=flow.Denyy0(:,xi);
%Denxx=zeros(l,1);
%Denxy=zeros(l,1);

        e= Rg*T/(r-1);
        et=Rg/(r-1)*ones(Ny,1);
        ett=zeros(Ny,1);
        er=zeros(Ny,1);
        err=zeros(Ny,1);
        ert=zeros(Ny,1);
        
        %p=Den*Rg.*T;
        pt=Den*Rg;
        %ptt=zeros(Ny,1);
        pr=Rg*T;
        %prr=zeros(Ny,1);
        prt=Rg*ones(Ny,1);
%sutherland公式        
       a1=110.4/parameter.Te;
       miu=(1.0+a1)*T.^1.5./(T+a1);
       miut=(1.0+a1)*sqrt(T).*(0.5*T+1.5*a1)./((T+a1).^2.0);
       miutt=(0.75*(T.^0.5+a1*T.^(-0.5)).*(T+a1).^2.0-(T.^1.5+3.0*a1*T.^0.5).*(T+a1)).*(1.0+a1)./(T+a1).^4.0;
       miur=zeros(Ny,1);
       miurr=zeros(Ny,1);
       miurt=zeros(Ny,1);
       
       a2=194/parameter.Te;
       %ka=(1.0+a2)*T.^1.5./(T+a2);
       kat=(1.0+a2)*sqrt(T).*(0.5*T+1.5*a2)./((T+a2).^2.0);
       katt=(0.75*(T.^0.5+a2*T.^(-0.5)).*(T+a2).^2.0-(T.^1.5+3.0*a2*T.^0.5).*(T+a2)).*(1.0+a2)./(T+a2).^4.0;
       kar=zeros(Ny,1);
       karr=zeros(Ny,1);
       kart=zeros(Ny,1);
       
       
%扰动项部分
Den2=zeros(Ny,1);u2=zeros(Ny,1);V2=zeros(Ny,1);W2=zeros(Ny,1);T2=zeros(Ny,1);
Denx2=zeros(Ny,1);Ux2=zeros(Ny,1);Vx2=zeros(Ny,1);Wx2=zeros(Ny,1);Tx2=zeros(Ny,1);
Deny2=zeros(Ny,1);Uy2=zeros(Ny,1);Vy2=zeros(Ny,1);Wy2=zeros(Ny,1);Ty2=zeros(Ny,1);
Denz2=zeros(Ny,1);Uz2=zeros(Ny,1);Vz2=zeros(Ny,1);Wz2=zeros(Ny,1);Tz2=zeros(Ny,1);
Dent2=zeros(Ny,1);Ut2=zeros(Ny,1);Vt2=zeros(Ny,1);Wt2=zeros(Ny,1);Tt2=zeros(Ny,1);
Denxx2=zeros(Ny,1);Uxx2=zeros(Ny,1);Vxx2=zeros(Ny,1);Wxx2=zeros(Ny,1);Txx2=zeros(Ny,1);
Denyy2=zeros(Ny,1);Uyy2=zeros(Ny,1);Vyy2=zeros(Ny,1);Wyy2=zeros(Ny,1);Tyy2=zeros(Ny,1);
Denzz2=zeros(Ny,1);Uzz2=zeros(Ny,1);Vzz2=zeros(Ny,1);Wzz2=zeros(Ny,1);Tzz2=zeros(Ny,1);
Denxy2=zeros(Ny,1);Uxy2=zeros(Ny,1);Vxy2=zeros(Ny,1);Wxy2=zeros(Ny,1);Txy2=zeros(Ny,1);
Denxz2=zeros(Ny,1);Uxz2=zeros(Ny,1);Vxz2=zeros(Ny,1);Wxz2=zeros(Ny,1);Txz2=zeros(Ny,1);
Denyz2=zeros(Ny,1);Uyz2=zeros(Ny,1);Vyz2=zeros(Ny,1);Wyz2=zeros(Ny,1);Tyz2=zeros(Ny,1);

for j=1:Ny
Den2(j)=F(5*j-4);u2(j)=F(5*j-3);V2(j)=F(5*j-2);W2(j)=F(5*j-1);T2(j)=F(5*j);
Denx2(j)=Fx(5*j-4);Ux2(j)=Fx(5*j-3);Vx2(j)=Fx(5*j-2);Wx2(j)=Fx(5*j-1);Tx2(j)=Fx(5*j);
Deny2(j)=Fy(5*j-4);Uy2(j)=Fy(5*j-3);Vy2(j)=Fy(5*j-2);Wy2(j)=Fy(5*j-1);Ty2(j)=Fy(5*j);
Denz2(j)=Fz(5*j-4);Uz2(j)=Fz(5*j-3);Vz2(j)=Fz(5*j-2);Wz2(j)=Fz(5*j-1);Tz2(j)=Fz(5*j);
Dent2(j)=Ft(5*j-4);Ut2(j)=Ft(5*j-3);Vt2(j)=Ft(5*j-2);Wt2(j)=Ft(5*j-1);Tt2(j)=Ft(5*j);
Denxx2(j)=Fxx(5*j-4);Uxx2(j)=Fxx(5*j-3);Vxx2(j)=Fxx(5*j-2);Wxx2(j)=Fxx(5*j-1);Txx2(j)=Fxx(5*j);
Denyy2(j)=Fyy(5*j-4);Uyy2(j)=Fyy(5*j-3);Vyy2(j)=Fyy(5*j-2);Wyy2(j)=Fyy(5*j-1);Tyy2(j)=Fyy(5*j);
Denzz2(j)=Fzz(5*j-4);Uzz2(j)=Fzz(5*j-3);Vzz2(j)=Fzz(5*j-2);Wzz2(j)=Fzz(5*j-1);Tzz2(j)=Fzz(5*j);
Denxy2(j)=Fxy(5*j-4);Uxy2(j)=Fxy(5*j-3);Vxy2(j)=Fxy(5*j-2);Wxy2(j)=Fxy(5*j-1);Txy2(j)=Fxy(5*j);
Denxz2(j)=Fxz(5*j-4);Uxz2(j)=Fxz(5*j-3);Vxz2(j)=Fxz(5*j-2);Wxz2(j)=Fxz(5*j-1);Txz2(j)=Fxz(5*j);
Denyz2(j)=Fyz(5*j-4);Uyz2(j)=Fyz(5*j-3);Vyz2(j)=Fyz(5*j-2);Wyz2(j)=Fyz(5*j-1);Tyz2(j)=Fyz(5*j);    
end

nonlin=zeros(Ny,5);

nonlin(:,1)=-(Denx2.*h1.^(-1).*u2+Den2.*h1.^(-1).*Ux2+Deny2.*V2+Den2.*h1.^(-1).*h11.*V2+Den2.*Vy2+Denz2.*W2+Den2.*Wz2);

nonlin(:,2)=-(Denx2.*h1.^(-1).*prt.*T2+Den2.*h1.^(-1).*prt.*Tx2+Dent2.*u2+Deny2.*h1.^(-1).*h11.*miur.*RE.^(-1).*u2+Den2.*h1.^(-2).*h11.^2.*miur.*RE.^(-1).*u2+Den2.* ...
  Deny.*h1.^(-1).*h11.*miurr.*RE.^(-1).*u2+Deny.*h1.^(-1).*h11.*miurt.*RE.^(-1).*T2.*u2+h1.^(-2).*h11.^2.*miut.*RE.^(-1).*T2.*u2+Den2.*h1.^(-1).*h11.* ...
  miurt.*RE.^(-1).*Ty.*u2+h1.^(-1).*h11.*miutt.*RE.^(-1).*T2.*Ty.*u2+h1.^(-1).*h11.*miut.*RE.^(-1).*Ty2.*u2+2.*Denx2.*h1.^(-1).*U.*u2+Denx.*h1.^(-1).* ...
  u2.^2+Denx2.*h1.^(-1).*u2.^2+Den2.*Ut2+2.*Den2.*h1.^(-1).*u2.*Ux+(-4/3).*Denx2.*h1.^(-2).*miur.*RE.^(-1).*Ux2+(-4/3).*Den2.*Denx.*h1.^(-2).*miurr.* ...
  RE.^(-1).*Ux2+(-4/3).*Denx.*h1.^(-2).*miurt.*RE.^(-1).*T2.*Ux2+(-4/3).*Den2.*h1.^(-2).*miurt.*RE.^(-1).*Tx.*Ux2+(-4/3).*h1.^(-2).*miutt.*RE.^(-1).* ...
  T2.*Tx.*Ux2+(-4/3).*h1.^(-2).*miut.*RE.^(-1).*Tx2.*Ux2+2.*Den2.*h1.^(-1).*U.*Ux2+2.*Den.*h1.^(-1).*u2.*Ux2+2.*Den2.*h1.^(-1).*u2.*Ux2+(-4/3).*Den2.* ...
  h1.^(-2).*miur.*RE.^(-1).*Uxx2+(-4/3).*h1.^(-2).*miut.*RE.^(-1).*T2.*Uxx2+(-1).*Deny2.*miur.*RE.^(-1).*Uy2+(-1).*Den2.*h1.^(-1).*h11.*miur.*RE.^(-1).* ...
  Uy2+(-1).*Den2.*Deny.*miurr.*RE.^(-1).*Uy2+(-1).*Deny.*miurt.*RE.^(-1).*T2.*Uy2+(-1).*h1.^(-1).*h11.*miut.*RE.^(-1).*T2.*Uy2+(-1).*Den2.*miurt.*RE.^( ...
  -1).*Ty.*Uy2+(-1).*miutt.*RE.^(-1).*T2.*Ty.*Uy2+(-1).*miut.*RE.^(-1).*Ty2.*Uy2+(-1).*Den2.*miur.*RE.^(-1).*Uyy2+(-1).*miut.*RE.^(-1).*T2.*Uyy2+(-1).* ...
  Denz2.*miur.*RE.^(-1).*Uz2+(-1).*miut.*RE.^(-1).*Tz2.*Uz2+(-1).*Den2.*miur.*RE.^(-1).*Uzz2+(-1).*miut.*RE.^(-1).*T2.*Uzz2+Deny2.*u2.*V+2.*Den2.*h1.^( ...
  -1).*h11.*u2.*V+Den2.*Uy2.*V+(-4/3).*Denx2.*h1.^(-2).*h11.*miur.*RE.^(-1).*V2+(-4/3).*Den2.*Denx.*h1.^(-2).*h11.*miurr.*RE.^(-1).*V2+(-4/3).*Denx.* ...
  h1.^(-2).*h11.*miurt.*RE.^(-1).*T2.*V2+(-4/3).*Den2.*h1.^(-2).*h11.*miurt.*RE.^(-1).*Tx.*V2+(-4/3).*h1.^(-2).*h11.*miutt.*RE.^(-1).*T2.*Tx.*V2+(-4/3) ...
  .*h1.^(-2).*h11.*miut.*RE.^(-1).*Tx2.*V2+Deny2.*U.*V2+2.*Den2.*h1.^(-1).*h11.*U.*V2+Deny.*u2.*V2+Deny2.*u2.*V2+2.*Den.*h1.^(-1).*h11.*u2.*V2+2.*Den2.* ...
  h1.^(-1).*h11.*u2.*V2+Den2.*Uy.*V2+Den.*Uy2.*V2+Den2.*Uy2.*V2+(-1).*Deny2.*h1.^(-1).*miur.*RE.^(-1).*Vx2+(-7/3).*Den2.*h1.^(-2).*h11.*miur.*RE.^(-1).* ...
  Vx2+(-1).*Den2.*Deny.*h1.^(-1).*miurr.*RE.^(-1).*Vx2+(-1).*Deny.*h1.^(-1).*miurt.*RE.^(-1).*T2.*Vx2+(-7/3).*h1.^(-2).*h11.*miut.*RE.^(-1).*T2.*Vx2+( ...
  -1).*Den2.*h1.^(-1).*miurt.*RE.^(-1).*Ty.*Vx2+(-1).*h1.^(-1).*miutt.*RE.^(-1).*T2.*Ty.*Vx2+(-1).*h1.^(-1).*miut.*RE.^(-1).*Ty2.*Vx2+(-1/3).*Den2.* ...
  h1.^(-1).*miur.*RE.^(-1).*Vxy2+(-1/3).*h1.^(-1).*miut.*RE.^(-1).*T2.*Vxy2+Den2.*u2.*Vy+(2/3).*Denx2.*h1.^(-1).*miur.*RE.^(-1).*Vy2+(2/3).*Den2.*Denx.* ...
  h1.^(-1).*miurr.*RE.^(-1).*Vy2+(2/3).*Denx.*h1.^(-1).*miurt.*RE.^(-1).*T2.*Vy2+(2/3).*Den2.*h1.^(-1).*miurt.*RE.^(-1).*Tx.*Vy2+(2/3).*h1.^(-1).* ...
  miutt.*RE.^(-1).*T2.*Tx.*Vy2+(2/3).*h1.^(-1).*miut.*RE.^(-1).*Tx2.*Vy2+Den2.*U.*Vy2+Den.*u2.*Vy2+Den2.*u2.*Vy2+Denz2.*U.*W2+Denz2.*u2.*W2+Den.*Uz2.* ...
  W2+Den2.*Uz2.*W2+(-1).*Denz2.*h1.^(-1).*miur.*RE.^(-1).*Wx2+(-1).*h1.^(-1).*miut.*RE.^(-1).*Tz2.*Wx2+(-1/3).*Den2.*h1.^(-1).*miur.*RE.^(-1).*Wxz2+( ...
  -1/3).*h1.^(-1).*miut.*RE.^(-1).*T2.*Wxz2+(2/3).*Denx2.*h1.^(-1).*miur.*RE.^(-1).*Wz2+(2/3).*Den2.*Denx.*h1.^(-1).*miurr.*RE.^(-1).*Wz2+(2/3).*Denx.* ...
  h1.^(-1).*miurt.*RE.^(-1).*T2.*Wz2+(2/3).*Den2.*h1.^(-1).*miurt.*RE.^(-1).*Tx.*Wz2+(2/3).*h1.^(-1).*miutt.*RE.^(-1).*T2.*Tx.*Wz2+(2/3).*h1.^(-1).* ...
  miut.*RE.^(-1).*Tx2.*Wz2+Den2.*U.*Wz2+Den.*u2.*Wz2+Den2.*u2.*Wz2);

nonlin(:,3)=-(Deny2.*prt.*T2+Den2.*prt.*Ty2+Denx2.*h1.^(-2).*h11.*miur.*RE.^(-1).*u2+Den2.*Denx.*h1.^(-2).*h11.*miurr.*RE.^(-1).*u2+Denx.*h1.^(-2).*h11.*miurt.* ...
  RE.^(-1).*T2.*u2+Den2.*h1.^(-2).*h11.*miurt.*RE.^(-1).*Tx.*u2+h1.^(-2).*h11.*miutt.*RE.^(-1).*T2.*Tx.*u2+h1.^(-2).*h11.*miut.*RE.^(-1).*Tx2.*u2+(-2).* ...
  Den2.*h1.^(-1).*h11.*U.*u2+(-1).*Den.*h1.^(-1).*h11.*u2.^2+(-1).*Den2.*h1.^(-1).*h11.*u2.^2+(2/3).*Deny2.*h1.^(-1).*miur.*RE.^(-1).*Ux2+(7/3).*Den2.* ...
  h1.^(-2).*h11.*miur.*RE.^(-1).*Ux2+(2/3).*Den2.*Deny.*h1.^(-1).*miurr.*RE.^(-1).*Ux2+(2/3).*Deny.*h1.^(-1).*miurt.*RE.^(-1).*T2.*Ux2+(7/3).*h1.^(-2).* ...
  h11.*miut.*RE.^(-1).*T2.*Ux2+(2/3).*Den2.*h1.^(-1).*miurt.*RE.^(-1).*Ty.*Ux2+(2/3).*h1.^(-1).*miutt.*RE.^(-1).*T2.*Ty.*Ux2+(2/3).*h1.^(-1).*miut.* ...
  RE.^(-1).*Ty2.*Ux2+(-1/3).*Den2.*h1.^(-1).*miur.*RE.^(-1).*Uxy2+(-1/3).*h1.^(-1).*miut.*RE.^(-1).*T2.*Uxy2+(-1).*Denx2.*h1.^(-1).*miur.*RE.^(-1).*Uy2+ ...
  (-1).*Den2.*Denx.*h1.^(-1).*miurr.*RE.^(-1).*Uy2+(-1).*Denx.*h1.^(-1).*miurt.*RE.^(-1).*T2.*Uy2+(-1).*Den2.*h1.^(-1).*miurt.*RE.^(-1).*Tx.*Uy2+(-1).* ...
  h1.^(-1).*miutt.*RE.^(-1).*T2.*Tx.*Uy2+(-1).*h1.^(-1).*miut.*RE.^(-1).*Tx2.*Uy2+Denx2.*h1.^(-1).*u2.*V+Den2.*h1.^(-1).*Ux2.*V+Dent2.*V2+(2/3).*Deny2.* ...
  h1.^(-1).*h11.*miur.*RE.^(-1).*V2+(4/3).*Den2.*h1.^(-2).*h11.^2.*miur.*RE.^(-1).*V2+(2/3).*Den2.*Deny.*h1.^(-1).*h11.*miurr.*RE.^(-1).*V2+(2/3).* ...
  Deny.*h1.^(-1).*h11.*miurt.*RE.^(-1).*T2.*V2+(4/3).*h1.^(-2).*h11.^2.*miut.*RE.^(-1).*T2.*V2+(2/3).*Den2.*h1.^(-1).*h11.*miurt.*RE.^(-1).*Ty.*V2+(2/3) ...
  .*h1.^(-1).*h11.*miutt.*RE.^(-1).*T2.*Ty.*V2+(2/3).*h1.^(-1).*h11.*miut.*RE.^(-1).*Ty2.*V2+Denx2.*h1.^(-1).*U.*V2+Denx.*h1.^(-1).*u2.*V2+Denx2.*h1.^( ...
  -1).*u2.*V2+Den2.*h1.^(-1).*Ux.*V2+Den.*h1.^(-1).*Ux2.*V2+Den2.*h1.^(-1).*Ux2.*V2+2.*Deny2.*V.*V2+2.*Den2.*h1.^(-1).*h11.*V.*V2+Deny.*V2.^2+Deny2.* ...
  V2.^2+Den.*h1.^(-1).*h11.*V2.^2+Den2.*h1.^(-1).*h11.*V2.^2+Den2.*Vt2+Den2.*h1.^(-1).*u2.*Vx+(-1).*Denx2.*h1.^(-2).*miur.*RE.^(-1).*Vx2+(-1).*Den2.* ...
  Denx.*h1.^(-2).*miurr.*RE.^(-1).*Vx2+(-1).*Denx.*h1.^(-2).*miurt.*RE.^(-1).*T2.*Vx2+(-1).*Den2.*h1.^(-2).*miurt.*RE.^(-1).*Tx.*Vx2+(-1).*h1.^(-2).* ...
  miutt.*RE.^(-1).*T2.*Tx.*Vx2+(-1).*h1.^(-2).*miut.*RE.^(-1).*Tx2.*Vx2+Den2.*h1.^(-1).*U.*Vx2+Den.*h1.^(-1).*u2.*Vx2+Den2.*h1.^(-1).*u2.*Vx2+(-1).* ...
  Den2.*h1.^(-2).*miur.*RE.^(-1).*Vxx2+(-1).*h1.^(-2).*miut.*RE.^(-1).*T2.*Vxx2+2.*Den2.*V2.*Vy+(-4/3).*Deny2.*miur.*RE.^(-1).*Vy2+(-4/3).*Den2.*h1.^( ...
  -1).*h11.*miur.*RE.^(-1).*Vy2+(-4/3).*Den2.*Deny.*miurr.*RE.^(-1).*Vy2+(-4/3).*Deny.*miurt.*RE.^(-1).*T2.*Vy2+(-4/3).*h1.^(-1).*h11.*miut.*RE.^(-1).* ...
  T2.*Vy2+(-4/3).*Den2.*miurt.*RE.^(-1).*Ty.*Vy2+(-4/3).*miutt.*RE.^(-1).*T2.*Ty.*Vy2+(-4/3).*miut.*RE.^(-1).*Ty2.*Vy2+2.*Den2.*V.*Vy2+2.*Den.*V2.*Vy2+ ...
  2.*Den2.*V2.*Vy2+(-4/3).*Den2.*miur.*RE.^(-1).*Vyy2+(-4/3).*miut.*RE.^(-1).*T2.*Vyy2+(-1).*Denz2.*miur.*RE.^(-1).*Vz2+(-1).*miut.*RE.^(-1).*Tz2.*Vz2+( ...
  -1).*Den2.*miur.*RE.^(-1).*Vzz2+(-1).*miut.*RE.^(-1).*T2.*Vzz2+Denz2.*V.*W2+Denz2.*V2.*W2+Den.*Vz2.*W2+Den2.*Vz2.*W2+(-1).*Denz2.*miur.*RE.^(-1).*Wy2+ ...
  (-1).*miut.*RE.^(-1).*Tz2.*Wy2+(-1/3).*Den2.*miur.*RE.^(-1).*Wyz2+(-1/3).*miut.*RE.^(-1).*T2.*Wyz2+(2/3).*Deny2.*miur.*RE.^(-1).*Wz2+(2/3).*Den2.* ...
  Deny.*miurr.*RE.^(-1).*Wz2+(2/3).*Deny.*miurt.*RE.^(-1).*T2.*Wz2+(2/3).*Den2.*miurt.*RE.^(-1).*Ty.*Wz2+(2/3).*miutt.*RE.^(-1).*T2.*Ty.*Wz2+(2/3).* ...
  miut.*RE.^(-1).*Ty2.*Wz2+Den2.*V.*Wz2+Den.*V2.*Wz2+Den2.*V2.*Wz2);


nonlin(:,4)=-(Denz2.*prt.*T2+Den2.*prt.*Tz2+(2/3).*Denz2.*h1.^(-1).*miur.*RE.^(-1).*Ux2+(2/3).*h1.^(-1).*miut.*RE.^(-1).*Tz2.*Ux2+(-1/3).*Den2.*h1.^(-1).*miur.* ...
  RE.^(-1).*Uxz2+(-1/3).*h1.^(-1).*miut.*RE.^(-1).*T2.*Uxz2+(-1).*Denx2.*h1.^(-1).*miur.*RE.^(-1).*Uz2+(-1).*Den2.*Denx.*h1.^(-1).*miurr.*RE.^(-1).*Uz2+ ...
  (-1).*Denx.*h1.^(-1).*miurt.*RE.^(-1).*T2.*Uz2+(-1).*Den2.*h1.^(-1).*miurt.*RE.^(-1).*Tx.*Uz2+(-1).*h1.^(-1).*miutt.*RE.^(-1).*T2.*Tx.*Uz2+(-1).*h1.^( ...
  -1).*miut.*RE.^(-1).*Tx2.*Uz2+(2/3).*Denz2.*h1.^(-1).*h11.*miur.*RE.^(-1).*V2+(2/3).*h1.^(-1).*h11.*miut.*RE.^(-1).*Tz2.*V2+(2/3).*Denz2.*miur.*RE.^( ...
  -1).*Vy2+(2/3).*miut.*RE.^(-1).*Tz2.*Vy2+(-1/3).*Den2.*miur.*RE.^(-1).*Vyz2+(-1/3).*miut.*RE.^(-1).*T2.*Vyz2+(-1).*Deny2.*miur.*RE.^(-1).*Vz2+(-1/3).* ...
  Den2.*h1.^(-1).*h11.*miur.*RE.^(-1).*Vz2+(-1).*Den2.*Deny.*miurr.*RE.^(-1).*Vz2+(-1).*Deny.*miurt.*RE.^(-1).*T2.*Vz2+(-1/3).*h1.^(-1).*h11.*miut.* ...
  RE.^(-1).*T2.*Vz2+(-1).*Den2.*miurt.*RE.^(-1).*Ty.*Vz2+(-1).*miutt.*RE.^(-1).*T2.*Ty.*Vz2+(-1).*miut.*RE.^(-1).*Ty2.*Vz2+Dent2.*W2+Denx2.*h1.^(-1).* ...
  U.*W2+Denx.*h1.^(-1).*u2.*W2+Denx2.*h1.^(-1).*u2.*W2+Den2.*h1.^(-1).*Ux.*W2+Den.*h1.^(-1).*Ux2.*W2+Den2.*h1.^(-1).*Ux2.*W2+Deny2.*V.*W2+Den2.*h1.^(-1) ...
  .*h11.*V.*W2+Deny.*V2.*W2+Deny2.*V2.*W2+Den.*h1.^(-1).*h11.*V2.*W2+Den2.*h1.^(-1).*h11.*V2.*W2+Den2.*Vy.*W2+Den.*Vy2.*W2+Den2.*Vy2.*W2+Denz2.*W2.^2+ ...
  Den2.*Wt2+(-1).*Denx2.*h1.^(-2).*miur.*RE.^(-1).*Wx2+(-1).*Den2.*Denx.*h1.^(-2).*miurr.*RE.^(-1).*Wx2+(-1).*Denx.*h1.^(-2).*miurt.*RE.^(-1).*T2.*Wx2+( ...
  -1).*Den2.*h1.^(-2).*miurt.*RE.^(-1).*Tx.*Wx2+(-1).*h1.^(-2).*miutt.*RE.^(-1).*T2.*Tx.*Wx2+(-1).*h1.^(-2).*miut.*RE.^(-1).*Tx2.*Wx2+Den2.*h1.^(-1).* ...
  U.*Wx2+Den.*h1.^(-1).*u2.*Wx2+Den2.*h1.^(-1).*u2.*Wx2+(-1).*Den2.*h1.^(-2).*miur.*RE.^(-1).*Wxx2+(-1).*h1.^(-2).*miut.*RE.^(-1).*T2.*Wxx2+(-1).* ...
  Deny2.*miur.*RE.^(-1).*Wy2+(-1).*Den2.*h1.^(-1).*h11.*miur.*RE.^(-1).*Wy2+(-1).*Den2.*Deny.*miurr.*RE.^(-1).*Wy2+(-1).*Deny.*miurt.*RE.^(-1).*T2.*Wy2+ ...
  (-1).*h1.^(-1).*h11.*miut.*RE.^(-1).*T2.*Wy2+(-1).*Den2.*miurt.*RE.^(-1).*Ty.*Wy2+(-1).*miutt.*RE.^(-1).*T2.*Ty.*Wy2+(-1).*miut.*RE.^(-1).*Ty2.*Wy2+ ...
  Den2.*V.*Wy2+Den.*V2.*Wy2+Den2.*V2.*Wy2+(-1).*Den2.*miur.*RE.^(-1).*Wyy2+(-1).*miut.*RE.^(-1).*T2.*Wyy2+(-4/3).*Denz2.*miur.*RE.^(-1).*Wz2+(-4/3).* ...
  miut.*RE.^(-1).*Tz2.*Wz2+2.*Den.*W2.*Wz2+2.*Den2.*W2.*Wz2+(-4/3).*Den2.*miur.*RE.^(-1).*Wzz2+(-4/3).*miut.*RE.^(-1).*T2.*Wzz2);

nonlin(:,5)=-(2.*Den2.*Dent2.*er+Dent2.*et.*T2+Den2.*et.*Tt2+(-1).*Denx2.*EC.^(-1).*h1.^(-2).*kar.*PR.^(-1).*RE.^(-1).*Tx2+(-1).*Den2.*Denx.*EC.^(-1).*h1.^(-2).* ...
  karr.*PR.^(-1).*RE.^(-1).*Tx2+(-1).*Denx.*EC.^(-1).*h1.^(-2).*kart.*PR.^(-1).*RE.^(-1).*T2.*Tx2+(-1).*Den2.*EC.^(-1).*h1.^(-2).*kart.*PR.^(-1).*RE.^( ...
  -1).*Tx.*Tx2+(-1).*EC.^(-1).*h1.^(-2).*katt.*PR.^(-1).*RE.^(-1).*T2.*Tx.*Tx2+(-1).*EC.^(-1).*h1.^(-2).*kat.*PR.^(-1).*RE.^(-1).*Tx2.^2+(-1).*Den2.* ...
  EC.^(-1).*h1.^(-2).*kar.*PR.^(-1).*RE.^(-1).*Txx2+(-1).*EC.^(-1).*h1.^(-2).*kat.*PR.^(-1).*RE.^(-1).*T2.*Txx2+(-1).*Deny2.*EC.^(-1).*kar.*PR.^(-1).* ...
  RE.^(-1).*Ty2+(-1).*Den2.*EC.^(-1).*h1.^(-1).*h11.*kar.*PR.^(-1).*RE.^(-1).*Ty2+(-1).*Den2.*Deny.*EC.^(-1).*karr.*PR.^(-1).*RE.^(-1).*Ty2+(-1).*Deny.* ...
  EC.^(-1).*kart.*PR.^(-1).*RE.^(-1).*T2.*Ty2+(-1).*EC.^(-1).*h1.^(-1).*h11.*kat.*PR.^(-1).*RE.^(-1).*T2.*Ty2+(-1).*Den2.*EC.^(-1).*kart.*PR.^(-1).* ...
  RE.^(-1).*Ty.*Ty2+(-1).*EC.^(-1).*katt.*PR.^(-1).*RE.^(-1).*T2.*Ty.*Ty2+(-1).*EC.^(-1).*kat.*PR.^(-1).*RE.^(-1).*Ty2.^2+(-1).*Den2.*EC.^(-1).*kar.* ...
  PR.^(-1).*RE.^(-1).*Tyy2+(-1).*EC.^(-1).*kat.*PR.^(-1).*RE.^(-1).*T2.*Tyy2+(-1).*Denz2.*EC.^(-1).*kar.*PR.^(-1).*RE.^(-1).*Tz2+(-1).*EC.^(-1).*kat.* ...
  PR.^(-1).*RE.^(-1).*Tz2.^2+(-1).*Den2.*EC.^(-1).*kar.*PR.^(-1).*RE.^(-1).*Tzz2+(-1).*EC.^(-1).*kat.*PR.^(-1).*RE.^(-1).*T2.*Tzz2+2.*Den2.*Denx2.*er.* ...
  h1.^(-1).*U+Den2.^2.*Denx.*err.*h1.^(-1).*U+Den2.*Denx.*ert.*h1.^(-1).*T2.*U+Denx2.*et.*h1.^(-1).*T2.*U+Den2.^2.*ert.*h1.^(-1).*Tx.*U+Den2.*ett.*h1.^( ...
  -1).*T2.*Tx.*U+Den2.*et.*h1.^(-1).*Tx2.*U+Denx2.*e.*h1.^(-1).*u2+2.*Den2.*Denx.*er.*h1.^(-1).*u2+Den.*Denx2.*er.*h1.^(-1).*u2+2.*Den2.*Denx2.*er.* ...
  h1.^(-1).*u2+Den.*Den2.*Denx.*err.*h1.^(-1).*u2+Den2.^2.*Denx.*err.*h1.^(-1).*u2+Den.*Denx.*ert.*h1.^(-1).*T2.*u2+Den2.*Denx.*ert.*h1.^(-1).*T2.*u2+ ...
  Denx.*et.*h1.^(-1).*T2.*u2+Denx2.*et.*h1.^(-1).*T2.*u2+Den.*Den2.*ert.*h1.^(-1).*Tx.*u2+Den2.^2.*ert.*h1.^(-1).*Tx.*u2+Den2.*et.*h1.^(-1).*Tx.*u2+ ...
  Den.*ett.*h1.^(-1).*T2.*Tx.*u2+Den2.*ett.*h1.^(-1).*T2.*Tx.*u2+Den.*et.*h1.^(-1).*Tx2.*u2+Den2.*et.*h1.^(-1).*Tx2.*u2+(-2).*Den2.*h1.^(-2).*h11.^2.* ...
  miur.*RE.^(-1).*U.*u2+(-2).*h1.^(-2).*h11.^2.*miut.*RE.^(-1).*T2.*U.*u2+(-1).*h1.^(-2).*h11.^2.*miu.*RE.^(-1).*u2.^2+(-1).*Den2.*h1.^(-2).*h11.^2.* ...
  miur.*RE.^(-1).*u2.^2+(-1).*h1.^(-2).*h11.^2.*miut.*RE.^(-1).*T2.*u2.^2+Den2.^2.*er.*h1.^(-1).*Ux+Den2.*et.*h1.^(-1).*T2.*Ux+Den2.*h1.^(-1).*prt.*T2.* ...
  Ux+Den2.*e.*h1.^(-1).*Ux2+Den.*Den2.*er.*h1.^(-1).*Ux2+Den2.^2.*er.*h1.^(-1).*Ux2+Den2.*h1.^(-1).*pr.*Ux2+Den.*et.*h1.^(-1).*T2.*Ux2+Den2.*et.*h1.^( ...
  -1).*T2.*Ux2+Den2.*h1.^(-1).*prt.*T2.*Ux2+h1.^(-1).*pt.*T2.*Ux2+(-8/3).*Den2.*h1.^(-2).*miur.*RE.^(-1).*Ux.*Ux2+(-8/3).*h1.^(-2).*miut.*RE.^(-1).*T2.* ...
  Ux.*Ux2+(-4/3).*h1.^(-2).*miu.*RE.^(-1).*Ux2.^2+(-4/3).*Den2.*h1.^(-2).*miur.*RE.^(-1).*Ux2.^2+(-4/3).*h1.^(-2).*miut.*RE.^(-1).*T2.*Ux2.^2+2.*Den2.* ...
  h1.^(-1).*h11.*miur.*RE.^(-1).*u2.*Uy+2.*h1.^(-1).*h11.*miut.*RE.^(-1).*T2.*u2.*Uy+2.*Den2.*h1.^(-1).*h11.*miur.*RE.^(-1).*U.*Uy2+2.*h1.^(-1).*h11.* ...
  miut.*RE.^(-1).*T2.*U.*Uy2+2.*h1.^(-1).*h11.*miu.*RE.^(-1).*u2.*Uy2+2.*Den2.*h1.^(-1).*h11.*miur.*RE.^(-1).*u2.*Uy2+2.*h1.^(-1).*h11.*miut.*RE.^(-1).* ...
  T2.*u2.*Uy2+(-2).*Den2.*miur.*RE.^(-1).*Uy.*Uy2+(-2).*miut.*RE.^(-1).*T2.*Uy.*Uy2+(-1).*miu.*RE.^(-1).*Uy2.^2+(-1).*Den2.*miur.*RE.^(-1).*Uy2.^2+(-1) ...
  .*miut.*RE.^(-1).*T2.*Uy2.^2+(-1).*miu.*RE.^(-1).*Uz2.^2+(-1).*Den2.*miur.*RE.^(-1).*Uz2.^2+(-1).*miut.*RE.^(-1).*T2.*Uz2.^2+2.*Den2.*Deny2.*er.*V+ ...
  Den2.^2.*Deny.*err.*V+Den2.^2.*er.*h1.^(-1).*h11.*V+Den2.*Deny.*ert.*T2.*V+Deny2.*et.*T2.*V+Den2.*et.*h1.^(-1).*h11.*T2.*V+Den2.*h1.^(-1).*h11.*prt.* ...
  T2.*V+Den2.^2.*ert.*Ty.*V+Den2.*ett.*T2.*Ty.*V+Den2.*et.*Ty2.*V+(-8/3).*Den2.*h1.^(-2).*h11.*miur.*RE.^(-1).*Ux2.*V+(-8/3).*h1.^(-2).*h11.*miut.*RE.^( ...
  -1).*T2.*Ux2.*V+Deny2.*e.*V2+2.*Den2.*Deny.*er.*V2+Den.*Deny2.*er.*V2+2.*Den2.*Deny2.*er.*V2+Den.*Den2.*Deny.*err.*V2+Den2.^2.*Deny.*err.*V2+Den2.*e.* ...
  h1.^(-1).*h11.*V2+Den.*Den2.*er.*h1.^(-1).*h11.*V2+Den2.^2.*er.*h1.^(-1).*h11.*V2+Den2.*h1.^(-1).*h11.*pr.*V2+Den.*Deny.*ert.*T2.*V2+Den2.*Deny.*ert.* ...
  T2.*V2+Deny.*et.*T2.*V2+Deny2.*et.*T2.*V2+Den.*et.*h1.^(-1).*h11.*T2.*V2+Den2.*et.*h1.^(-1).*h11.*T2.*V2+Den2.*h1.^(-1).*h11.*prt.*T2.*V2+h1.^(-1).* ...
  h11.*pt.*T2.*V2+Den.*Den2.*ert.*Ty.*V2+Den2.^2.*ert.*Ty.*V2+Den2.*et.*Ty.*V2+Den.*ett.*T2.*Ty.*V2+Den2.*ett.*T2.*Ty.*V2+Den.*et.*Ty2.*V2+Den2.*et.* ...
  Ty2.*V2+(-8/3).*Den2.*h1.^(-2).*h11.*miur.*RE.^(-1).*Ux.*V2+(-8/3).*h1.^(-2).*h11.*miut.*RE.^(-1).*T2.*Ux.*V2+(-8/3).*h1.^(-2).*h11.*miu.*RE.^(-1).* ...
  Ux2.*V2+(-8/3).*Den2.*h1.^(-2).*h11.*miur.*RE.^(-1).*Ux2.*V2+(-8/3).*h1.^(-2).*h11.*miut.*RE.^(-1).*T2.*Ux2.*V2+(-8/3).*Den2.*h1.^(-2).*h11.^2.*miur.* ...
  RE.^(-1).*V.*V2+(-8/3).*h1.^(-2).*h11.^2.*miut.*RE.^(-1).*T2.*V.*V2+(-4/3).*h1.^(-2).*h11.^2.*miu.*RE.^(-1).*V2.^2+(-4/3).*Den2.*h1.^(-2).*h11.^2.* ...
  miur.*RE.^(-1).*V2.^2+(-4/3).*h1.^(-2).*h11.^2.*miut.*RE.^(-1).*T2.*V2.^2+2.*Den2.*h1.^(-2).*h11.*miur.*RE.^(-1).*u2.*Vx+2.*h1.^(-2).*h11.*miut.*RE.^( ...
  -1).*T2.*u2.*Vx+(-2).*Den2.*h1.^(-1).*miur.*RE.^(-1).*Uy2.*Vx+(-2).*h1.^(-1).*miut.*RE.^(-1).*T2.*Uy2.*Vx+2.*Den2.*h1.^(-2).*h11.*miur.*RE.^(-1).*U.* ...
  Vx2+2.*h1.^(-2).*h11.*miut.*RE.^(-1).*T2.*U.*Vx2+2.*h1.^(-2).*h11.*miu.*RE.^(-1).*u2.*Vx2+2.*Den2.*h1.^(-2).*h11.*miur.*RE.^(-1).*u2.*Vx2+2.*h1.^(-2) ...
  .*h11.*miut.*RE.^(-1).*T2.*u2.*Vx2+(-2).*Den2.*h1.^(-1).*miur.*RE.^(-1).*Uy.*Vx2+(-2).*h1.^(-1).*miut.*RE.^(-1).*T2.*Uy.*Vx2+(-2).*h1.^(-1).*miu.* ...
  RE.^(-1).*Uy2.*Vx2+(-2).*Den2.*h1.^(-1).*miur.*RE.^(-1).*Uy2.*Vx2+(-2).*h1.^(-1).*miut.*RE.^(-1).*T2.*Uy2.*Vx2+(-2).*Den2.*h1.^(-2).*miur.*RE.^(-1).* ...
  Vx.*Vx2+(-2).*h1.^(-2).*miut.*RE.^(-1).*T2.*Vx.*Vx2+(-1).*h1.^(-2).*miu.*RE.^(-1).*Vx2.^2+(-1).*Den2.*h1.^(-2).*miur.*RE.^(-1).*Vx2.^2+(-1).*h1.^(-2) ...
  .*miut.*RE.^(-1).*T2.*Vx2.^2+Den2.^2.*er.*Vy+Den2.*et.*T2.*Vy+Den2.*prt.*T2.*Vy+(4/3).*Den2.*h1.^(-1).*miur.*RE.^(-1).*Ux2.*Vy+(4/3).*h1.^(-1).*miut.* ...
  RE.^(-1).*T2.*Ux2.*Vy+(4/3).*Den2.*h1.^(-1).*h11.*miur.*RE.^(-1).*V2.*Vy+(4/3).*h1.^(-1).*h11.*miut.*RE.^(-1).*T2.*V2.*Vy+Den2.*e.*Vy2+Den.*Den2.*er.* ...
  Vy2+Den2.^2.*er.*Vy2+Den2.*pr.*Vy2+Den.*et.*T2.*Vy2+Den2.*et.*T2.*Vy2+Den2.*prt.*T2.*Vy2+pt.*T2.*Vy2+(4/3).*Den2.*h1.^(-1).*miur.*RE.^(-1).*Ux.*Vy2+( ...
  4/3).*h1.^(-1).*miut.*RE.^(-1).*T2.*Ux.*Vy2+(4/3).*h1.^(-1).*miu.*RE.^(-1).*Ux2.*Vy2+(4/3).*Den2.*h1.^(-1).*miur.*RE.^(-1).*Ux2.*Vy2+(4/3).*h1.^(-1).* ...
  miut.*RE.^(-1).*T2.*Ux2.*Vy2+(4/3).*Den2.*h1.^(-1).*h11.*miur.*RE.^(-1).*V.*Vy2+(4/3).*h1.^(-1).*h11.*miut.*RE.^(-1).*T2.*V.*Vy2+(4/3).*h1.^(-1).* ...
  h11.*miu.*RE.^(-1).*V2.*Vy2+(4/3).*Den2.*h1.^(-1).*h11.*miur.*RE.^(-1).*V2.*Vy2+(4/3).*h1.^(-1).*h11.*miut.*RE.^(-1).*T2.*V2.*Vy2+(-8/3).*Den2.*miur.* ...
  RE.^(-1).*Vy.*Vy2+(-8/3).*miut.*RE.^(-1).*T2.*Vy.*Vy2+(-4/3).*miu.*RE.^(-1).*Vy2.^2+(-4/3).*Den2.*miur.*RE.^(-1).*Vy2.^2+(-4/3).*miut.*RE.^(-1).*T2.* ...
  Vy2.^2+(-1).*miu.*RE.^(-1).*Vz2.^2+(-1).*Den2.*miur.*RE.^(-1).*Vz2.^2+(-1).*miut.*RE.^(-1).*T2.*Vz2.^2+Denz2.*e.*W2+Den.*Denz2.*er.*W2+2.*Den2.* ...
  Denz2.*er.*W2+Denz2.*et.*T2.*W2+Den.*et.*Tz2.*W2+Den2.*et.*Tz2.*W2+(-2).*h1.^(-1).*miu.*RE.^(-1).*Uz2.*Wx2+(-2).*Den2.*h1.^(-1).*miur.*RE.^(-1).*Uz2.* ...
  Wx2+(-2).*h1.^(-1).*miut.*RE.^(-1).*T2.*Uz2.*Wx2+(-1).*h1.^(-2).*miu.*RE.^(-1).*Wx2.^2+(-1).*Den2.*h1.^(-2).*miur.*RE.^(-1).*Wx2.^2+(-1).*h1.^(-2).* ...
  miut.*RE.^(-1).*T2.*Wx2.^2+(-2).*miu.*RE.^(-1).*Vz2.*Wy2+(-2).*Den2.*miur.*RE.^(-1).*Vz2.*Wy2+(-2).*miut.*RE.^(-1).*T2.*Vz2.*Wy2+(-1).*miu.*RE.^(-1).* ...
  Wy2.^2+(-1).*Den2.*miur.*RE.^(-1).*Wy2.^2+(-1).*miut.*RE.^(-1).*T2.*Wy2.^2+Den2.*e.*Wz2+Den.*Den2.*er.*Wz2+Den2.^2.*er.*Wz2+Den2.*pr.*Wz2+Den.*et.* ...
  T2.*Wz2+Den2.*et.*T2.*Wz2+Den2.*prt.*T2.*Wz2+pt.*T2.*Wz2+(4/3).*Den2.*h1.^(-1).*miur.*RE.^(-1).*Ux.*Wz2+(4/3).*h1.^(-1).*miut.*RE.^(-1).*T2.*Ux.*Wz2+( ...
  4/3).*h1.^(-1).*miu.*RE.^(-1).*Ux2.*Wz2+(4/3).*Den2.*h1.^(-1).*miur.*RE.^(-1).*Ux2.*Wz2+(4/3).*h1.^(-1).*miut.*RE.^(-1).*T2.*Ux2.*Wz2+(4/3).*Den2.* ...
  h1.^(-1).*h11.*miur.*RE.^(-1).*V.*Wz2+(4/3).*h1.^(-1).*h11.*miut.*RE.^(-1).*T2.*V.*Wz2+(4/3).*h1.^(-1).*h11.*miu.*RE.^(-1).*V2.*Wz2+(4/3).*Den2.*h1.^( ...
  -1).*h11.*miur.*RE.^(-1).*V2.*Wz2+(4/3).*h1.^(-1).*h11.*miut.*RE.^(-1).*T2.*V2.*Wz2+(4/3).*Den2.*miur.*RE.^(-1).*Vy.*Wz2+(4/3).*miut.*RE.^(-1).*T2.* ...
  Vy.*Wz2+(4/3).*miu.*RE.^(-1).*Vy2.*Wz2+(4/3).*Den2.*miur.*RE.^(-1).*Vy2.*Wz2+(4/3).*miut.*RE.^(-1).*T2.*Vy2.*Wz2+(-4/3).*miu.*RE.^(-1).*Wz2.^2+(-4/3) ...
  .*Den2.*miur.*RE.^(-1).*Wz2.^2+(-4/3).*miut.*RE.^(-1).*T2.*Wz2.^2);

end