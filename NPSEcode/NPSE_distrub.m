function [distrub0]=NPSE_distrub(xi,m,n,sz,st,MESH,NPSE)

Fai=NPSE.Fai;
alf=NPSE.alf;

m_max=NPSE.m_max;
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

        distrub0.F=zeros(5*Ny,1);
        distrub0.Ft=zeros(5*Ny,1);
        distrub0.Fx=zeros(5*Ny,1);
        distrub0.Fy=zeros(5*Ny,1);
        distrub0.Fz=zeros(5*Ny,1);
        distrub0.Fxx=zeros(5*Ny,1);
        distrub0.Fyy=zeros(5*Ny,1);
        distrub0.Fzz=zeros(5*Ny,1);
        distrub0.Fxy=zeros(5*Ny,1);
        distrub0.Fxz=zeros(5*Ny,1);
        distrub0.Fyz=zeros(5*Ny,1);



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
distrub0.F=F.*exp1;
distrub0.Ft=-zi*m*omega*F.*exp1;
distrub0.Fx=(Fx+F*zi*alfmn).*exp1;
distrub0.Fy=Fy.*exp1;
distrub0.Fz=zi*n*beta*F.*exp1;
distrub0.Fxx=(Fxx+2*zi*Fx*alfmn+F*zi*alfx-F*alfmn^2).*exp1;
distrub0.Fyy=Fyy.*exp1;
distrub0.Fzz=(-n^2*beta^2*F).*exp1;
distrub0.Fxy=(Fxy+Fy*zi*alfmn).*exp1;
distrub0.Fxz=(Fx*zi*n*beta-F*alfmn*n*beta).*exp1;
distrub0.Fyz=(Fy*zi*n*beta).*exp1;
else
distrub0.F=2*real(F.*exp1);
distrub0.Ft=2*real(-zi*m*omega*F.*exp1);
distrub0.Fx=2*real((Fx+F*zi*alfmn).*exp1);
distrub0.Fy=2*real(Fy.*exp1);
distrub0.Fz=2*real(zi*n*beta*F.*exp1);
distrub0.Fxx=2*real((Fxx+2*zi*Fx*alfmn+F*zi*alfx-F*alfmn^2).*exp1);
distrub0.Fyy=2*real(Fyy.*exp1);
distrub0.Fzz=2*real((-n^2*beta^2*F).*exp1);
distrub0.Fxy=2*real((Fxy+Fy*zi*alfmn).*exp1);
distrub0.Fxz=2*real((Fx*zi*n*beta-F*alfmn*n*beta).*exp1);
distrub0.Fyz=2*real((Fy*zi*n*beta).*exp1);    
end

end