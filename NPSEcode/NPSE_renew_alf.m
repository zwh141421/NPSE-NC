function[Residual]=NPSE_renew_alf(xi,m,n,Fai,y,Ny,NPSE)

m_max = NPSE(1);
n_max = NPSE(2);
dx    = NPSE(8);

uvw      = zeros(Ny,1);
uvw_all  = 0;
duvw_all = 0;

ux   = zeros(Ny,1);
vx   = zeros(Ny,1);
wx   = zeros(Ny,1);
duvw = zeros(Ny,1);
for i=1:Ny
    
    u  = Fai(5*(i-1)+2,(xi-1)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n));
    u1 = Fai(5*(i-1)+2,(xi-2)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n));

    v  = Fai(5*(i-1)+3,(xi-1)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n)); 
    v1 = Fai(5*(i-1)+3,(xi-2)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n));

    w  = Fai(5*(i-1)+4,(xi-1)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n)); 
    w1 = Fai(5*(i-1)+4,(xi-2)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n));
    
uvw(i) = abs(u)^2+abs(v)^2+abs(w)^2;
ux(i)  = (u-u1)/dx;
vx(i)  = (v-v1)/dx;
wx(i)  = (w-w1)/dx;

duvw(i)=conj(u)*ux(i)+conj(v)*vx(i)+conj(w)*wx(i);
if i>=2
    uvw_all  = uvw_all+(uvw(i)+uvw(i-1))*(y(i)-y(i-1))/2;
    duvw_all = duvw_all+(duvw(i)+duvw(i-1))*(y(i)-y(i-1))/2;   
end
    
end
Residual=duvw_all/uvw_all;


end