function[Residual]=NPSE_renew_alf(xi,m,n,Fai,MESH,NPSE)
Ny=MESH.Ny;
y=MESH.y;
n_max=NPSE.n_max;
dx=NPSE.dx;
uvw=zeros(Ny,1);
uvw_all=0;
duvw_all=0;

ux=zeros(Ny,1);
vx=zeros(Ny,1);
wx=zeros(Ny,1);
duvw=zeros(Ny,1);
for i=1:Ny
    
    u=Fai(5*(i-1)+2,1+m,n_max+n+1,xi); u1=Fai(5*(i-1)+2,1+m,n_max+n+1,xi-1);
    v=Fai(5*(i-1)+3,1+m,n_max+n+1,xi); v1=Fai(5*(i-1)+3,1+m,n_max+n+1,xi-1);
    w=Fai(5*(i-1)+4,1+m,n_max+n+1,xi); w1=Fai(5*(i-1)+4,1+m,n_max+n+1,xi-1);
    
uvw(i)=abs(u)^2+abs(v)^2+abs(w)^2;
ux(i)=(u-u1)/dx;
vx(i)=(v-v1)/dx;
wx(i)=(w-w1)/dx;

duvw(i)=conj(u)*ux(i)+conj(v)*vx(i)+conj(w)*wx(i);
if i>=2
    uvw_all=uvw_all+(uvw(i)+uvw(i-1))*(y(i)-y(i-1))/2;
    duvw_all=duvw_all+(duvw(i)+duvw(i-1))*(y(i)-y(i-1))/2;   
end
    
end
Residual=duvw_all/uvw_all;


end