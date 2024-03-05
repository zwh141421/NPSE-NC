function [F0,Ft0,Fx0,Fy0,Fz0,Fxx0,Fyy0,Fzz0,Fxy0,Fxz0,Fyz0]=NPSE_disturb(xi,m,n,sz,st,X,alf,NPSE,f,fx,fxx,fy,fyy,fxy)

m_max = NPSE(1);
n_max = NPSE(2);
omega = NPSE(4);
beta  = NPSE(5);
x=X;
alfmn = alf(1,(xi-1)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n));
zz=sqrt(-1);


JF=trapz(x(1:xi),alf(1,m*(2*n_max+1)+(n_max+1+n):(1+m_max)*(2*n_max+1):(xi-1)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n)));

exp1=exp(zz*JF+zz*n*beta*sz-zz*m*omega*st);

alfx=0;

if m==0
F0=f (:,m*(2*n_max+1)+n+n_max+1).*exp1;
Ft0=-zz*m*omega*f (:,m*(2*n_max+1)+n+n_max+1).*exp1;
Fx0=(fx(:,m*(2*n_max+1)+n+n_max+1)+f(:,m*(2*n_max+1)+n+n_max+1)*zz*alfmn).*exp1;
Fy0=fy (:,m*(2*n_max+1)+n+n_max+1).*exp1;
Fz0=zz*n*beta*f (:,m*(2*n_max+1)+n+n_max+1).*exp1;
Fxx0=(fxx(:,m*(2*n_max+1)+n+n_max+1)+2*zz*fx(:,m*(2*n_max+1)+n+n_max+1)*alfmn+f(:,m*(2*n_max+1)+n+n_max+1)*zz*alfx-f(:,m*(2*n_max+1)+n+n_max+1)*alfmn^2).*exp1;
Fyy0=fyy (:,m*(2*n_max+1)+n+n_max+1).*exp1;
Fzz0=(-n^2*beta^2*f(:,m*(2*n_max+1)+n+n_max+1)).*exp1;
Fxy0=(fxy(:,m*(2*n_max+1)+n+n_max+1)+fy(:,m*(2*n_max+1)+n+n_max+1)*zz*alfmn).*exp1;
Fxz0=(fx(:,m*(2*n_max+1)+n+n_max+1)*zz*n*beta-f(:,m*(2*n_max+1)+n+n_max+1)*alfmn*n*beta).*exp1;
Fyz0=(fy(:,m*(2*n_max+1)+n+n_max+1)*zz*n*beta).*exp1;
else
F0=2*real(f (:,m*(2*n_max+1)+n+n_max+1).*exp1);
Ft0=2*real(-zz*m*omega*f (:,m*(2*n_max+1)+n+n_max+1).*exp1);
Fx0=2*real((fx(:,m*(2*n_max+1)+n+n_max+1)+f(:,m*(2*n_max+1)+n+n_max+1)*zz*alfmn).*exp1);
Fy0=2*real(fy (:,m*(2*n_max+1)+n+n_max+1).*exp1);
Fz0=2*real(zz*n*beta*f (:,m*(2*n_max+1)+n+n_max+1).*exp1);
Fxx0=2*real((fxx(:,m*(2*n_max+1)+n+n_max+1)+2*zz*fx(:,m*(2*n_max+1)+n+n_max+1)*alfmn+f(:,m*(2*n_max+1)+n+n_max+1)*zz*alfx-f(:,m*(2*n_max+1)+n+n_max+1)*alfmn^2).*exp1);
Fyy0=2*real(fyy (:,m*(2*n_max+1)+n+n_max+1).*exp1);
Fzz0=2*real((-n^2*beta^2*f(:,m*(2*n_max+1)+n+n_max+1)).*exp1);
Fxy0=2*real((fxy(:,m*(2*n_max+1)+n+n_max+1)+fy(:,m*(2*n_max+1)+n+n_max+1)*zz*alfmn).*exp1);
Fxz0=2*real((fx(:,m*(2*n_max+1)+n+n_max+1)*zz*n*beta-f(:,m*(2*n_max+1)+n+n_max+1)*alfmn*n*beta).*exp1);
Fyz0=2*real((fy(:,m*(2*n_max+1)+n+n_max+1)*zz*n*beta).*exp1);    
end

end
