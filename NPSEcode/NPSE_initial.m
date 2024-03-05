function [ flow,alf,Fai,modes ] = NPSE_initial(Nx,Ny,flow0,a1,b1,NPSE)

m=1;
n=0;
m_max = NPSE(1);
n_max = NPSE(2);


modes=zeros(1+m_max,2*n_max+1);           %The modes involved in the calculation


flow=zeros(Nx*Ny,19);
flow(:,1)=reshape(flow0.U0,[Nx*Ny,1]);
flow(:,2)=reshape(flow0.Ux0,[Nx*Ny,1]);
flow(:,3)=reshape(flow0.Uy0,[Nx*Ny,1]);
flow(:,4)=reshape(flow0.Uxy0,[Nx*Ny,1]);
flow(:,5)=reshape(flow0.Uyy0,[Nx*Ny,1]);
flow(:,6)=reshape(flow0.V0,[Nx*Ny,1]);
flow(:,7)=reshape(flow0.Vx0,[Nx*Ny,1]);
flow(:,8)=reshape(flow0.Vy0,[Nx*Ny,1]);
flow(:,9)=reshape(flow0.Vyy0,[Nx*Ny,1]);
flow(:,10)=reshape(flow0.Vxy0,[Nx*Ny,1]);
flow(:,11)=reshape(flow0.T0,[Nx*Ny,1]);
flow(:,12)=reshape(flow0.Tx0,[Nx*Ny,1]);
flow(:,13)=reshape(flow0.Ty0,[Nx*Ny,1]);
flow(:,14)=reshape(flow0.Tyy0,[Nx*Ny,1]);
flow(:,15)=reshape(flow0.Den0,[Nx*Ny,1]);
flow(:,16)=reshape(flow0.Denx0,[Nx*Ny,1]);
flow(:,17)=reshape(flow0.Deny0,[Nx*Ny,1]);
flow(:,18)=reshape(flow0.Denyy0,[Nx*Ny,1]);
flow(1:Nx,19)=flow0.R;


Fai=zeros(5*Ny,Nx*(1+m_max)*(2*n_max+1));
alf=zeros(1,Nx*(1+m_max)*(2*n_max+1));


Fai(:,1*(1+m)*(n_max+1+n))=b1*NPSE(3);
alf(1,1*(1+m)*(n_max+1+n))=a1;
modes(1+m,n_max+1+n)=1;
end
