function [ alf,Fai,modes ] = NPSE_initial(MESH,NPSE)

m=1;
n=0;
m_max=NPSE.m_max;
n_max=NPSE.n_max;
a1=NPSE.a1;
b1=NPSE.b1;
%计算初值
modes=zeros(1+m_max,2*n_max+1);           %参与计算的模态

Fai=zeros(5*MESH.Ny,1+m_max,2*n_max+1,MESH.Nx);
alf=zeros(1+m_max,2*n_max+1,MESH.Nx);

%for i=1:MESH.Ny
 % b1(5*(i-1)+4)=0;  
%end

Fai(:,1+m,n_max+1+n,1)=b1*NPSE.Amp;
alf(1+m,n_max+1+n,1)=a1;
modes(1+m,n_max+1+n)=1;
end