function [F0,Ft0,Fx0,Fy0,Fz0,Fxx0,Fyy0,Fzz0,Fxy0,Fxz0,Fyz0]=NPSE_distrub(xi,m,n,sz,st,MESH,NPSE)

Fai=NPSE.Fai;
alf=NPSE.alf;

%m_max=NPSE.m_max;
n_max=NPSE.n_max;
alfmn=alf(1+m,n_max+1+n,xi);
omega=NPSE.omega;
beta=NPSE.beta;
x=NPSE.X;
dx=NPSE.dx;
Ny=MESH.Ny;
Dy=MESH.D11;
Dyy=MESH.D22;

zi=sqrt(-1);

DD1=kron(Dy,eye(5,5));  DD2=kron(Dyy,eye(5,5)); 

JF=trapz(x(1:xi),alf(1+m,n_max+1+n,1:xi));

exp1=exp(zi*JF+zi*n*beta*sz-zi*m*omega*st);

F=Fai(:,m+1,n+n_max+1,xi);
Fx=(Fai(:,m+1,n+n_max+1,xi)-Fai(:,m+1,n+n_max+1,xi-1))/dx;
Fxx=0;
Fy=DD1*Fai(:,m+1,n+n_max+1,xi);%%%%%%%%%0.01s
Fyy=DD2*Fai(:,m+1,n+n_max+1,xi);%%%%%%%%%%%%0.01s
Fxy=(DD1*Fai(:,m+1,n+n_max+1,xi)-DD1*Fai(:,m+1,n+n_max+1,xi-1))/dx;%%%%%%%%%%%0.02s

alfx=0;

if m==0
F0=F.*exp1;
Ft0=-zi*m*omega*F.*exp1;
Fx0=(Fx+F*zi*alfmn).*exp1;
Fy0=Fy.*exp1;
Fz0=zi*n*beta*F.*exp1;
Fxx0=(Fxx+2*zi*Fx*alfmn+F*zi*alfx-F*alfmn^2).*exp1;
Fyy0=Fyy.*exp1;
Fzz0=(-n^2*beta^2*F).*exp1;
Fxy0=(Fxy+Fy*zi*alfmn).*exp1;
Fxz0=(Fx*zi*n*beta-F*alfmn*n*beta).*exp1;
Fyz0=(Fy*zi*n*beta).*exp1;
else
F0=2*real(F.*exp1);
Ft0=2*real(-zi*m*omega*F.*exp1);
Fx0=2*real((Fx+F*zi*alfmn).*exp1);
Fy0=2*real(Fy.*exp1);
Fz0=2*real(zi*n*beta*F.*exp1);
Fxx0=2*real((Fxx+2*zi*Fx*alfmn+F*zi*alfx-F*alfmn^2).*exp1);
Fyy0=2*real(Fyy.*exp1);
Fzz0=2*real((-n^2*beta^2*F).*exp1);
Fxy0=2*real((Fxy+Fy*zi*alfmn).*exp1);
Fxz0=2*real((Fx*zi*n*beta-F*alfmn*n*beta).*exp1);
Fyz0=2*real((Fy*zi*n*beta).*exp1);    
end

end