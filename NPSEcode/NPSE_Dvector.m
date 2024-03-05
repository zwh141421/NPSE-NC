function [f,fx,fxx,fy,fyy,fxy]=NPSE_Dvector(xi,m,n,DD1,DD2,Fai,NPSE)

m_max = NPSE(1);
n_max = NPSE(2);

dx    = NPSE(8);



f=Fai(:,(xi-1)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n));
fx=(Fai(:,(xi-1)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n))-Fai(:,(xi-2)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n)))/dx;
fxx=0;
fy=DD1*Fai(:,(xi-1)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n));%%%%%%%%% 0.01s
fyy=DD2*Fai(:,(xi-1)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n));%%%%%%%%%%%% 0.01s
fxy=(DD1*Fai(:,(xi-1)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n))-DD1*Fai(:,(xi-2)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n)))/dx;%%%%%%%%%%% 0.02s


end
