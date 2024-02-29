function[gamma,A,A2,B,C,D,Hxx,Hyy,Hzz,Hxy,Hxz,Hyz]=NPSE_matrix_PSE(i,xi,MESH,flow)

global parameter
           
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


Ny=MESH.Ny;

U=flow.U0(:,xi);   
Ux=flow.Ux0(:,xi);
Uy=flow.Uy0(:,xi);
Uxy=flow.Uxy0(:,xi);
Uyy=flow.Uyy0(:,xi);
Uxx=zeros(Ny,1);

V=flow.V0(:,xi);
Vx=flow.Vx0(:,xi);
Vy=flow.Vy0(:,xi);
Vxy=flow.Vxy0(:,xi);
Vyy=flow.Vyy0(:,xi);
Vxx=zeros(Ny,1);

T=flow.T0(:,xi);
Tx=flow.Tx0(:,xi);
Ty=flow.Ty0(:,xi);
Tyy=flow.Tyy0(:,xi);
Txx=zeros(Ny,1);

Den=flow.Den0(:,xi);
Denx=flow.Denx0(:,xi);
Deny=flow.Deny0(:,xi);
%Denyy=flow.Denyy0(:,xi);
%Denxx=zeros(Ny,1);
%Denxy=zeros(Ny,1);

%理想气体
        e= Rg*T/(r-1);
        et=Rg/(r-1)*ones(Ny,1);
        ett=zeros(Ny,1);
        er=zeros(Ny,1);
        err=zeros(Ny,1);
        ert=zeros(Ny,1);
        
        p=Den*Rg.*T;
        pt=Den*Rg;
        ptt=zeros(Ny,1);
        pr=Rg*T;
        prr=zeros(Ny,1);
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
       ka=(1.0+a2)*T.^1.5./(T+a2);
       kat=(1.0+a2)*sqrt(T).*(0.5*T+1.5*a2)./((T+a2).^2.0);
       katt=(0.75*(T.^0.5+a2*T.^(-0.5)).*(T+a2).^2.0-(T.^1.5+3.0*a2*T.^0.5).*(T+a2)).*(1.0+a2)./(T+a2).^4.0;
       kar=zeros(Ny,1);
       karr=zeros(Ny,1);
       kart=zeros(Ny,1);



gamma=[1,0,0,0,0;U(i),Den(i),0,0,0;V(i),0,Den(i),0,0;0,0,0,Den(i),0;e(i)+Den(i).*er(i),0,0,0,Den(i).*et(i)];

A=[h1.^(-1).*U(i),h1.^(-1).*Den(i),0,0,0;h1.^(-1).*pr(i)+h1.^(-1).*U(i).^2+(-4/3).*h1.^(-2).*RE.^(-1).*miur(i).*Ux(i)+(-4/3).*h1.^(-2).*h11.*RE.^(-1).* ...
  miur(i).*V(i)+(2/3).*h1.^(-1).*RE.^(-1).*miur(i).*Vy(i),(-4/3).*h1.^(-2).*RE.^(-1).*Denx(i).*miur(i)+(-4/3).*h1.^(-2).*RE.^(-1).*miut(i).*Tx(i)+2.* ...
  h1.^(-1).*Den(i).*U(i),(-7/3).*h1.^(-2).*h11.*RE.^(-1).*miu(i)+(-1).*h1.^(-1).*RE.^(-1).*Deny(i).*miur(i)+(-1).*h1.^(-1).*RE.^(-1).*miut(i).*Ty(i),0, ...
  h1.^(-1).*pt(i)+(-4/3).*h1.^(-2).*RE.^(-1).*miut(i).*Ux(i)+(-4/3).*h1.^(-2).*h11.*RE.^(-1).*miut(i).*V(i)+(2/3).*h1.^(-1).*RE.^(-1).*miut(i).*Vy(i); ...
  h1.^(-2).*h11.*RE.^(-1).*miur(i).*U(i)+(-1).*h1.^(-1).*RE.^(-1).*miur(i).*Uy(i)+h1.^(-1).*U(i).*V(i)+(-1).*h1.^(-2).*RE.^(-1).*miur(i).*Vx(i),(7/3).* ...
  h1.^(-2).*h11.*RE.^(-1).*miu(i)+(2/3).*h1.^(-1).*RE.^(-1).*Deny(i).*miur(i)+(2/3).*h1.^(-1).*RE.^(-1).*miut(i).*Ty(i)+h1.^(-1).*Den(i).*V(i),(-1).* ...
  h1.^(-2).*RE.^(-1).*Denx(i).*miur(i)+(-1).*h1.^(-2).*RE.^(-1).*miut(i).*Tx(i)+h1.^(-1).*Den(i).*U(i),0,h1.^(-2).*h11.*RE.^(-1).*miut(i).*U(i)+(-1).* ...
  h1.^(-1).*RE.^(-1).*miut(i).*Uy(i)+(-1).*h1.^(-2).*RE.^(-1).*miut(i).*Vx(i);0,0,0,(-1).*h1.^(-2).*RE.^(-1).*Denx(i).*miur(i)+(-1).*h1.^(-2).*RE.^(-1) ...
  .*miut(i).*Tx(i)+h1.^(-1).*Den(i).*U(i),0;(-1).*EC.^(-1).*h1.^(-2).*PR.^(-1).*RE.^(-1).*kar(i).*Tx(i)+h1.^(-1).*e(i).*U(i)+h1.^(-1).*Den(i).*er(i).*U( ...
  i),h1.^(-1).*Den(i).*e(i)+h1.^(-1).*p(i)+(-8/3).*h1.^(-2).*RE.^(-1).*miu(i).*Ux(i)+(-8/3).*h1.^(-2).*h11.*RE.^(-1).*miu(i).*V(i)+(4/3).*h1.^(-1).* ...
  RE.^(-1).*miu(i).*Vy(i),2.*h1.^(-2).*h11.*RE.^(-1).*miu(i).*U(i)+(-2).*h1.^(-1).*RE.^(-1).*miu(i).*Uy(i)+(-2).*h1.^(-2).*RE.^(-1).*miu(i).*Vx(i),0,( ...
  -1).*EC.^(-1).*h1.^(-2).*PR.^(-1).*RE.^(-1).*Denx(i).*kar(i)+(-2).*EC.^(-1).*h1.^(-2).*PR.^(-1).*RE.^(-1).*kat(i).*Tx(i)+h1.^(-1).*Den(i).*et(i).*U(i) ...
  ];

%压力修正后的矩阵A
A2=[h1.^(-1).*U(i),h1.^(-1).*Den(i),0,0,0;h1.^(-1).*pr(i)+h1.^(-1).*U(i).^2+(-4/3).*h1.^(-2).*RE.^(-1).*miur(i).*Ux(i)+(-4/3).*h1.^(-2).*h11.*RE.^(-1).* ...
  miur(i).*V(i)+(2/3).*h1.^(-1).*RE.^(-1).*miur(i).*Vy(i),(-4/3).*h1.^(-2).*RE.^(-1).*Denx(i).*miur(i)+(-4/3).*h1.^(-2).*RE.^(-1).*miut(i).*Tx(i)+2.* ...
  h1.^(-1).*Den(i).*U(i),(-7/3).*h1.^(-2).*h11.*RE.^(-1).*miu(i)+(-1).*h1.^(-1).*RE.^(-1).*Deny(i).*miur(i)+(-1).*h1.^(-1).*RE.^(-1).*miut(i).*Ty(i),0, ...
  h1.^(-1).*pt(i)+(-4/3).*h1.^(-2).*RE.^(-1).*miut(i).*Ux(i)+(-4/3).*h1.^(-2).*h11.*RE.^(-1).*miut(i).*V(i)+(2/3).*h1.^(-1).*RE.^(-1).*miut(i).*Vy(i); ...
  h1.^(-2).*h11.*RE.^(-1).*miur(i).*U(i)+(-1).*h1.^(-1).*RE.^(-1).*miur(i).*Uy(i)+h1.^(-1).*U(i).*V(i)+(-1).*h1.^(-2).*RE.^(-1).*miur(i).*Vx(i),(7/3).* ...
  h1.^(-2).*h11.*RE.^(-1).*miu(i)+(2/3).*h1.^(-1).*RE.^(-1).*Deny(i).*miur(i)+(2/3).*h1.^(-1).*RE.^(-1).*miut(i).*Ty(i)+h1.^(-1).*Den(i).*V(i),(-1).* ...
  h1.^(-2).*RE.^(-1).*Denx(i).*miur(i)+(-1).*h1.^(-2).*RE.^(-1).*miut(i).*Tx(i)+h1.^(-1).*Den(i).*U(i),0,h1.^(-2).*h11.*RE.^(-1).*miut(i).*U(i)+(-1).* ...
  h1.^(-1).*RE.^(-1).*miut(i).*Uy(i)+(-1).*h1.^(-2).*RE.^(-1).*miut(i).*Vx(i);0,0,0,(-1).*h1.^(-2).*RE.^(-1).*Denx(i).*miur(i)+(-1).*h1.^(-2).*RE.^(-1) ...
  .*miut(i).*Tx(i)+h1.^(-1).*Den(i).*U(i),0;(-1).*EC.^(-1).*h1.^(-2).*PR.^(-1).*RE.^(-1).*kar(i).*Tx(i)+h1.^(-1).*e(i).*U(i)+h1.^(-1).*Den(i).*er(i).*U( ...
  i),h1.^(-1).*Den(i).*e(i)+h1.^(-1).*p(i)+(-8/3).*h1.^(-2).*RE.^(-1).*miu(i).*Ux(i)+(-8/3).*h1.^(-2).*h11.*RE.^(-1).*miu(i).*V(i)+(4/3).*h1.^(-1).* ...
  RE.^(-1).*miu(i).*Vy(i),2.*h1.^(-2).*h11.*RE.^(-1).*miu(i).*U(i)+(-2).*h1.^(-1).*RE.^(-1).*miu(i).*Uy(i)+(-2).*h1.^(-2).*RE.^(-1).*miu(i).*Vx(i),0,( ...
  -1).*EC.^(-1).*h1.^(-2).*PR.^(-1).*RE.^(-1).*Denx(i).*kar(i)+(-2).*EC.^(-1).*h1.^(-2).*PR.^(-1).*RE.^(-1).*kat(i).*Tx(i)+h1.^(-1).*Den(i).*et(i).*U(i) ...
  ];
A2(2,1)=A2(2,1)-h1.^(-1).*pr(i);
A2(2,5)=A2(2,5)-h1.^(-1).*pt(i);

B=[V(i),0,Den(i),0,0;h1.^(-1).*h11.*RE.^(-1).*miur(i).*U(i)+(-1).*RE.^(-1).*miur(i).*Uy(i)+U(i).*V(i)+(-1).*h1.^(-1).*RE.^(-1).*miur(i).*Vx(i),(-1).* ...
  h1.^(-1).*h11.*RE.^(-1).*miu(i)+(-1).*RE.^(-1).*Deny(i).*miur(i)+(-1).*RE.^(-1).*miut(i).*Ty(i)+Den(i).*V(i),(2/3).*h1.^(-1).*RE.^(-1).*Denx(i).*miur( ...
  i)+(2/3).*h1.^(-1).*RE.^(-1).*miut(i).*Tx(i)+Den(i).*U(i),0,h1.^(-1).*h11.*RE.^(-1).*miut(i).*U(i)+(-1).*RE.^(-1).*miut(i).*Uy(i)+(-1).*h1.^(-1).* ...
  RE.^(-1).*miut(i).*Vx(i);pr(i)+(2/3).*h1.^(-1).*RE.^(-1).*miur(i).*Ux(i)+(2/3).*h1.^(-1).*h11.*RE.^(-1).*miur(i).*V(i)+V(i).^2+(-4/3).*RE.^(-1).*miur( ...
  i).*Vy(i),(-1).*h1.^(-1).*RE.^(-1).*Denx(i).*miur(i)+(-1).*h1.^(-1).*RE.^(-1).*miut(i).*Tx(i),(-4/3).*h1.^(-1).*h11.*RE.^(-1).*miu(i)+(-4/3).*RE.^(-1) ...
  .*Deny(i).*miur(i)+(-4/3).*RE.^(-1).*miut(i).*Ty(i)+2.*Den(i).*V(i),0,pt(i)+(2/3).*h1.^(-1).*RE.^(-1).*miut(i).*Ux(i)+(2/3).*h1.^(-1).*h11.*RE.^(-1).* ...
  miut(i).*V(i)+(-4/3).*RE.^(-1).*miut(i).*Vy(i);0,0,0,(-1).*h1.^(-1).*h11.*RE.^(-1).*miu(i)+(-1).*RE.^(-1).*Deny(i).*miur(i)+(-1).*RE.^(-1).*miut(i).* ...
  Ty(i)+Den(i).*V(i),0;(-1).*EC.^(-1).*PR.^(-1).*RE.^(-1).*kar(i).*Ty(i)+e(i).*V(i)+Den(i).*er(i).*V(i),2.*h1.^(-1).*h11.*RE.^(-1).*miu(i).*U(i)+(-2).* ...
  RE.^(-1).*miu(i).*Uy(i)+(-2).*h1.^(-1).*RE.^(-1).*miu(i).*Vx(i),Den(i).*e(i)+p(i)+(4/3).*h1.^(-1).*RE.^(-1).*miu(i).*Ux(i)+(4/3).*h1.^(-1).*h11.*RE.^( ...
  -1).*miu(i).*V(i)+(-8/3).*RE.^(-1).*miu(i).*Vy(i),0,(-1).*EC.^(-1).*h1.^(-1).*h11.*PR.^(-1).*RE.^(-1).*ka(i)+(-1).*EC.^(-1).*PR.^(-1).*RE.^(-1).*Deny( ...
  i).*kar(i)+(-2).*EC.^(-1).*PR.^(-1).*RE.^(-1).*kat(i).*Ty(i)+Den(i).*et(i).*V(i)];

C=[0,0,0,Den(i),0;0,0,0,(2/3).*h1.^(-1).*RE.^(-1).*Denx(i).*miur(i)+(2/3).*h1.^(-1).*RE.^(-1).*miut(i).*Tx(i)+Den(i).*U(i),0;0,0,0,(2/3).*RE.^(-1).* ...
  Deny(i).*miur(i)+(2/3).*RE.^(-1).*miut(i).*Ty(i)+Den(i).*V(i),0;pr(i)+(2/3).*h1.^(-1).*RE.^(-1).*miur(i).*Ux(i)+(2/3).*h1.^(-1).*h11.*RE.^(-1).*miur( ...
  i).*V(i)+(2/3).*RE.^(-1).*miur(i).*Vy(i),(-1).*h1.^(-1).*RE.^(-1).*Denx(i).*miur(i)+(-1).*h1.^(-1).*RE.^(-1).*miut(i).*Tx(i),(-1/3).*h1.^(-1).*h11.* ...
  RE.^(-1).*miu(i)+(-1).*RE.^(-1).*Deny(i).*miur(i)+(-1).*RE.^(-1).*miut(i).*Ty(i),0,pt(i)+(2/3).*h1.^(-1).*RE.^(-1).*miut(i).*Ux(i)+(2/3).*h1.^(-1).* ...
  h11.*RE.^(-1).*miut(i).*V(i)+(2/3).*RE.^(-1).*miut(i).*Vy(i);0,0,0,Den(i).*e(i)+p(i)+(4/3).*h1.^(-1).*RE.^(-1).*miu(i).*Ux(i)+(4/3).*h1.^(-1).*h11.* ...
  RE.^(-1).*miu(i).*V(i)+(4/3).*RE.^(-1).*miu(i).*Vy(i),0];%

D=[h1.^(-1).*Ux(i)+h1.^(-1).*h11.*V(i)+Vy(i),h1.^(-1).*Denx(i),h1.^(-1).*h11.*Den(i)+Deny(i),0,0;h1.^(-1).*Denx(i).*prr(i)+h1.^(-1).*prt(i).*Tx(i)+h1.^( ...
  -2).*h11.^2.*RE.^(-1).*miur(i).*U(i)+h1.^(-1).*h11.*RE.^(-1).*Deny(i).*miurr(i).*U(i)+h1.^(-1).*h11.*RE.^(-1).*miurt(i).*Ty(i).*U(i)+(-4/3).*h1.^(-2) ...
  .*RE.^(-1).*Denx(i).*miurr(i).*Ux(i)+(-4/3).*h1.^(-2).*RE.^(-1).*miurt(i).*Tx(i).*Ux(i)+2.*h1.^(-1).*U(i).*Ux(i)+(-4/3).*h1.^(-2).*RE.^(-1).*miur(i).* ...
  Uxx(i)+(-1).*h1.^(-1).*h11.*RE.^(-1).*miur(i).*Uy(i)+(-1).*RE.^(-1).*Deny(i).*miurr(i).*Uy(i)+(-1).*RE.^(-1).*miurt(i).*Ty(i).*Uy(i)+(-1).*RE.^(-1).* ...
  miur(i).*Uyy(i)+(-4/3).*h1.^(-2).*h11.*RE.^(-1).*Denx(i).*miurr(i).*V(i)+(-4/3).*h1.^(-2).*h11.*RE.^(-1).*miurt(i).*Tx(i).*V(i)+2.*h1.^(-1).*h11.*U(i) ...
  .*V(i)+Uy(i).*V(i)+(-7/3).*h1.^(-2).*h11.*RE.^(-1).*miur(i).*Vx(i)+(-1).*h1.^(-1).*RE.^(-1).*Deny(i).*miurr(i).*Vx(i)+(-1).*h1.^(-1).*RE.^(-1).*miurt( ...
  i).*Ty(i).*Vx(i)+(-1/3).*h1.^(-1).*RE.^(-1).*miur(i).*Vxy(i)+(2/3).*h1.^(-1).*RE.^(-1).*Denx(i).*miurr(i).*Vy(i)+(2/3).*h1.^(-1).*RE.^(-1).*miurt(i).* ...
  Tx(i).*Vy(i)+U(i).*Vy(i),h1.^(-2).*h11.^2.*RE.^(-1).*miu(i)+h1.^(-1).*h11.*RE.^(-1).*Deny(i).*miur(i)+h1.^(-1).*h11.*RE.^(-1).*miut(i).*Ty(i)+2.*h1.^( ...
  -1).*Denx(i).*U(i)+2.*h1.^(-1).*Den(i).*Ux(i)+2.*h1.^(-1).*h11.*Den(i).*V(i)+Deny(i).*V(i)+Den(i).*Vy(i),(-4/3).*h1.^(-2).*h11.*RE.^(-1).*Denx(i).* ...
  miur(i)+(-4/3).*h1.^(-2).*h11.*RE.^(-1).*miut(i).*Tx(i)+2.*h1.^(-1).*h11.*Den(i).*U(i)+Deny(i).*U(i)+Den(i).*Uy(i),0,h1.^(-1).*Denx(i).*prt(i)+h1.^( ...
  -1).*ptt(i).*Tx(i)+h1.^(-1).*h11.*RE.^(-1).*Deny(i).*miurt(i).*U(i)+h1.^(-2).*h11.^2.*RE.^(-1).*miut(i).*U(i)+h1.^(-1).*h11.*RE.^(-1).*miutt(i).*Ty(i) ...
  .*U(i)+(-4/3).*h1.^(-2).*RE.^(-1).*Denx(i).*miurt(i).*Ux(i)+(-4/3).*h1.^(-2).*RE.^(-1).*miutt(i).*Tx(i).*Ux(i)+(-4/3).*h1.^(-2).*RE.^(-1).*miut(i).* ...
  Uxx(i)+(-1).*RE.^(-1).*Deny(i).*miurt(i).*Uy(i)+(-1).*h1.^(-1).*h11.*RE.^(-1).*miut(i).*Uy(i)+(-1).*RE.^(-1).*miutt(i).*Ty(i).*Uy(i)+(-1).*RE.^(-1).* ...
  miut(i).*Uyy(i)+(-4/3).*h1.^(-2).*h11.*RE.^(-1).*Denx(i).*miurt(i).*V(i)+(-4/3).*h1.^(-2).*h11.*RE.^(-1).*miutt(i).*Tx(i).*V(i)+(-1).*h1.^(-1).*RE.^( ...
  -1).*Deny(i).*miurt(i).*Vx(i)+(-7/3).*h1.^(-2).*h11.*RE.^(-1).*miut(i).*Vx(i)+(-1).*h1.^(-1).*RE.^(-1).*miutt(i).*Ty(i).*Vx(i)+(-1/3).*h1.^(-1).*RE.^( ...
  -1).*miut(i).*Vxy(i)+(2/3).*h1.^(-1).*RE.^(-1).*Denx(i).*miurt(i).*Vy(i)+(2/3).*h1.^(-1).*RE.^(-1).*miutt(i).*Tx(i).*Vy(i);Deny(i).*prr(i)+prt(i).*Ty( ...
  i)+h1.^(-2).*h11.*RE.^(-1).*Denx(i).*miurr(i).*U(i)+h1.^(-2).*h11.*RE.^(-1).*miurt(i).*Tx(i).*U(i)+(-1).*h1.^(-1).*h11.*U(i).^2+(7/3).*h1.^(-2).*h11.* ...
  RE.^(-1).*miur(i).*Ux(i)+(2/3).*h1.^(-1).*RE.^(-1).*Deny(i).*miurr(i).*Ux(i)+(2/3).*h1.^(-1).*RE.^(-1).*miurt(i).*Ty(i).*Ux(i)+(-1/3).*h1.^(-1).*RE.^( ...
  -1).*miur(i).*Uxy(i)+(-1).*h1.^(-1).*RE.^(-1).*Denx(i).*miurr(i).*Uy(i)+(-1).*h1.^(-1).*RE.^(-1).*miurt(i).*Tx(i).*Uy(i)+(4/3).*h1.^(-2).*h11.^2.* ...
  RE.^(-1).*miur(i).*V(i)+(2/3).*h1.^(-1).*h11.*RE.^(-1).*Deny(i).*miurr(i).*V(i)+(2/3).*h1.^(-1).*h11.*RE.^(-1).*miurt(i).*Ty(i).*V(i)+h1.^(-1).*Ux(i) ...
  .*V(i)+h1.^(-1).*h11.*V(i).^2+(-1).*h1.^(-2).*RE.^(-1).*Denx(i).*miurr(i).*Vx(i)+(-1).*h1.^(-2).*RE.^(-1).*miurt(i).*Tx(i).*Vx(i)+h1.^(-1).*U(i).*Vx( ...
  i)+(-1).*h1.^(-2).*RE.^(-1).*miur(i).*Vxx(i)+(-4/3).*h1.^(-1).*h11.*RE.^(-1).*miur(i).*Vy(i)+(-4/3).*RE.^(-1).*Deny(i).*miurr(i).*Vy(i)+(-4/3).*RE.^( ...
  -1).*miurt(i).*Ty(i).*Vy(i)+2.*V(i).*Vy(i)+(-4/3).*RE.^(-1).*miur(i).*Vyy(i),h1.^(-2).*h11.*RE.^(-1).*Denx(i).*miur(i)+h1.^(-2).*h11.*RE.^(-1).*miut( ...
  i).*Tx(i)+(-2).*h1.^(-1).*h11.*Den(i).*U(i)+h1.^(-1).*Denx(i).*V(i)+h1.^(-1).*Den(i).*Vx(i),(4/3).*h1.^(-2).*h11.^2.*RE.^(-1).*miu(i)+(2/3).*h1.^(-1) ...
  .*h11.*RE.^(-1).*Deny(i).*miur(i)+(2/3).*h1.^(-1).*h11.*RE.^(-1).*miut(i).*Ty(i)+h1.^(-1).*Denx(i).*U(i)+h1.^(-1).*Den(i).*Ux(i)+2.*h1.^(-1).*h11.* ...
  Den(i).*V(i)+2.*Deny(i).*V(i)+2.*Den(i).*Vy(i),0,Deny(i).*prt(i)+ptt(i).*Ty(i)+h1.^(-2).*h11.*RE.^(-1).*Denx(i).*miurt(i).*U(i)+h1.^(-2).*h11.*RE.^( ...
  -1).*miutt(i).*Tx(i).*U(i)+(2/3).*h1.^(-1).*RE.^(-1).*Deny(i).*miurt(i).*Ux(i)+(7/3).*h1.^(-2).*h11.*RE.^(-1).*miut(i).*Ux(i)+(2/3).*h1.^(-1).*RE.^( ...
  -1).*miutt(i).*Ty(i).*Ux(i)+(-1/3).*h1.^(-1).*RE.^(-1).*miut(i).*Uxy(i)+(-1).*h1.^(-1).*RE.^(-1).*Denx(i).*miurt(i).*Uy(i)+(-1).*h1.^(-1).*RE.^(-1).* ...
  miutt(i).*Tx(i).*Uy(i)+(2/3).*h1.^(-1).*h11.*RE.^(-1).*Deny(i).*miurt(i).*V(i)+(4/3).*h1.^(-2).*h11.^2.*RE.^(-1).*miut(i).*V(i)+(2/3).*h1.^(-1).*h11.* ...
  RE.^(-1).*miutt(i).*Ty(i).*V(i)+(-1).*h1.^(-2).*RE.^(-1).*Denx(i).*miurt(i).*Vx(i)+(-1).*h1.^(-2).*RE.^(-1).*miutt(i).*Tx(i).*Vx(i)+(-1).*h1.^(-2).* ...
  RE.^(-1).*miut(i).*Vxx(i)+(-4/3).*RE.^(-1).*Deny(i).*miurt(i).*Vy(i)+(-4/3).*h1.^(-1).*h11.*RE.^(-1).*miut(i).*Vy(i)+(-4/3).*RE.^(-1).*miutt(i).*Ty(i) ...
  .*Vy(i)+(-4/3).*RE.^(-1).*miut(i).*Vyy(i);0,0,0,h1.^(-1).*Denx(i).*U(i)+h1.^(-1).*Den(i).*Ux(i)+h1.^(-1).*h11.*Den(i).*V(i)+Deny(i).*V(i)+Den(i).*Vy( ...
  i),0;(-1).*EC.^(-1).*h1.^(-2).*PR.^(-1).*RE.^(-1).*Denx(i).*karr(i).*Tx(i)+(-1).*EC.^(-1).*h1.^(-2).*PR.^(-1).*RE.^(-1).*kart(i).*Tx(i).^2+(-1).*EC.^( ...
  -1).*h1.^(-2).*PR.^(-1).*RE.^(-1).*kar(i).*Txx(i)+(-1).*EC.^(-1).*h1.^(-1).*h11.*PR.^(-1).*RE.^(-1).*kar(i).*Ty(i)+(-1).*EC.^(-1).*PR.^(-1).*RE.^(-1) ...
  .*Deny(i).*karr(i).*Ty(i)+(-1).*EC.^(-1).*PR.^(-1).*RE.^(-1).*kart(i).*Ty(i).^2+(-1).*EC.^(-1).*PR.^(-1).*RE.^(-1).*kar(i).*Tyy(i)+2.*h1.^(-1).*Denx( ...
  i).*er(i).*U(i)+h1.^(-1).*Den(i).*Denx(i).*err(i).*U(i)+h1.^(-1).*Den(i).*ert(i).*Tx(i).*U(i)+h1.^(-1).*et(i).*Tx(i).*U(i)+(-1).*h1.^(-2).*h11.^2.* ...
  RE.^(-1).*miur(i).*U(i).^2+h1.^(-1).*e(i).*Ux(i)+h1.^(-1).*Den(i).*er(i).*Ux(i)+h1.^(-1).*pr(i).*Ux(i)+(-4/3).*h1.^(-2).*RE.^(-1).*miur(i).*Ux(i).^2+ ...
  2.*h1.^(-1).*h11.*RE.^(-1).*miur(i).*U(i).*Uy(i)+(-1).*RE.^(-1).*miur(i).*Uy(i).^2+h1.^(-1).*h11.*e(i).*V(i)+h1.^(-1).*h11.*Den(i).*er(i).*V(i)+2.* ...
  Deny(i).*er(i).*V(i)+Den(i).*Deny(i).*err(i).*V(i)+h1.^(-1).*h11.*pr(i).*V(i)+Den(i).*ert(i).*Ty(i).*V(i)+et(i).*Ty(i).*V(i)+(-8/3).*h1.^(-2).*h11.* ...
  RE.^(-1).*miur(i).*Ux(i).*V(i)+(-4/3).*h1.^(-2).*h11.^2.*RE.^(-1).*miur(i).*V(i).^2+2.*h1.^(-2).*h11.*RE.^(-1).*miur(i).*U(i).*Vx(i)+(-2).*h1.^(-1).* ...
  RE.^(-1).*miur(i).*Uy(i).*Vx(i)+(-1).*h1.^(-2).*RE.^(-1).*miur(i).*Vx(i).^2+e(i).*Vy(i)+Den(i).*er(i).*Vy(i)+pr(i).*Vy(i)+(4/3).*h1.^(-1).*RE.^(-1).* ...
  miur(i).*Ux(i).*Vy(i)+(4/3).*h1.^(-1).*h11.*RE.^(-1).*miur(i).*V(i).*Vy(i)+(-4/3).*RE.^(-1).*miur(i).*Vy(i).^2,h1.^(-1).*Denx(i).*e(i)+h1.^(-1).*Den( ...
  i).*Denx(i).*er(i)+h1.^(-1).*Den(i).*et(i).*Tx(i)+(-2).*h1.^(-2).*h11.^2.*RE.^(-1).*miu(i).*U(i)+2.*h1.^(-1).*h11.*RE.^(-1).*miu(i).*Uy(i)+2.*h1.^(-2) ...
  .*h11.*RE.^(-1).*miu(i).*Vx(i),h1.^(-1).*h11.*Den(i).*e(i)+Deny(i).*e(i)+Den(i).*Deny(i).*er(i)+h1.^(-1).*h11.*p(i)+Den(i).*et(i).*Ty(i)+(-8/3).*h1.^( ...
  -2).*h11.*RE.^(-1).*miu(i).*Ux(i)+(-8/3).*h1.^(-2).*h11.^2.*RE.^(-1).*miu(i).*V(i)+(4/3).*h1.^(-1).*h11.*RE.^(-1).*miu(i).*Vy(i),0,(-1).*EC.^(-1).* ...
  h1.^(-2).*PR.^(-1).*RE.^(-1).*Denx(i).*kart(i).*Tx(i)+(-1).*EC.^(-1).*h1.^(-2).*PR.^(-1).*RE.^(-1).*katt(i).*Tx(i).^2+(-1).*EC.^(-1).*h1.^(-2).*PR.^( ...
  -1).*RE.^(-1).*kat(i).*Txx(i)+(-1).*EC.^(-1).*PR.^(-1).*RE.^(-1).*Deny(i).*kart(i).*Ty(i)+(-1).*EC.^(-1).*h1.^(-1).*h11.*PR.^(-1).*RE.^(-1).*kat(i).* ...
  Ty(i)+(-1).*EC.^(-1).*PR.^(-1).*RE.^(-1).*katt(i).*Ty(i).^2+(-1).*EC.^(-1).*PR.^(-1).*RE.^(-1).*kat(i).*Tyy(i)+h1.^(-1).*Den(i).*Denx(i).*ert(i).*U(i)+ ...
  h1.^(-1).*Denx(i).*et(i).*U(i)+h1.^(-1).*Den(i).*ett(i).*Tx(i).*U(i)+(-1).*h1.^(-2).*h11.^2.*RE.^(-1).*miut(i).*U(i).^2+h1.^(-1).*Den(i).*et(i).*Ux( ...
  i)+h1.^(-1).*pt(i).*Ux(i)+(-4/3).*h1.^(-2).*RE.^(-1).*miut(i).*Ux(i).^2+2.*h1.^(-1).*h11.*RE.^(-1).*miut(i).*U(i).*Uy(i)+(-1).*RE.^(-1).*miut(i).*Uy( ...
  i).^2+Den(i).*Deny(i).*ert(i).*V(i)+h1.^(-1).*h11.*Den(i).*et(i).*V(i)+Deny(i).*et(i).*V(i)+h1.^(-1).*h11.*pt(i).*V(i)+Den(i).*ett(i).*Ty(i).*V(i)+( ...
  -8/3).*h1.^(-2).*h11.*RE.^(-1).*miut(i).*Ux(i).*V(i)+(-4/3).*h1.^(-2).*h11.^2.*RE.^(-1).*miut(i).*V(i).^2+2.*h1.^(-2).*h11.*RE.^(-1).*miut(i).*U(i).* ...
  Vx(i)+(-2).*h1.^(-1).*RE.^(-1).*miut(i).*Uy(i).*Vx(i)+(-1).*h1.^(-2).*RE.^(-1).*miut(i).*Vx(i).^2+Den(i).*et(i).*Vy(i)+pt(i).*Vy(i)+(4/3).*h1.^(-1).* ...
  RE.^(-1).*miut(i).*Ux(i).*Vy(i)+(4/3).*h1.^(-1).*h11.*RE.^(-1).*miut(i).*V(i).*Vy(i)+(-4/3).*RE.^(-1).*miut(i).*Vy(i).^2];

Hxx=-[0,0,0,0,0;0,(-4/3).*h1.^(-2).*RE.^(-1).*miu(i),0,0,0;0,0,(-1).*h1.^(-2).*RE.^(-1).*miu(i),0,0;0,0,0,(-1).*h1.^(-2).*RE.^(-1).*miu(i),0;0,0,0,0,(-1).* ...
  EC.^(-1).*h1.^(-2).*PR.^(-1).*RE.^(-1).*ka(i)];

Hyy=-[0,0,0,0,0;0,(-1).*RE.^(-1).*miu(i),0,0,0;0,0,(-4/3).*RE.^(-1).*miu(i),0,0;0,0,0,(-1).*RE.^(-1).*miu(i),0;0,0,0,0,(-1).*EC.^(-1).*PR.^(-1).*RE.^(-1).* ...
  ka(i)];

Hzz=-[0,0,0,0,0;0,(-1).*RE.^(-1).*miu(i),0,0,0;0,0,(-1).*RE.^(-1).*miu(i),0,0;0,0,0,(-4/3).*RE.^(-1).*miu(i),0;0,0,0,0,(-1).*EC.^(-1).*PR.^(-1).*RE.^(-1).* ...
  ka(i)];%

Hxy=-[0,0,0,0,0;0,0,(-1/3).*h1.^(-1).*RE.^(-1).*miu(i),0,0;0,(-1/3).*h1.^(-1).*RE.^(-1).*miu(i),0,0,0;0,0,0,0,0;0,0,0,0,0];

Hxz=-[0,0,0,0,0;0,0,0,(-1/3).*h1.^(-1).*RE.^(-1).*miu(i),0;0,0,0,0,0;0,(-1/3).*h1.^(-1).*RE.^(-1).*miu(i),0,0,0;0,0,0,0,0];%

Hyz=-[0,0,0,0,0;0,0,0,0,0;0,0,0,(-1/3).*RE.^(-1).*miu(i),0;0,0,(-1/3).*RE.^(-1).*miu(i),0,0;0,0,0,0,0];%

end