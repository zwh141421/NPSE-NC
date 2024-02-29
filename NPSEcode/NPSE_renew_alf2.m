function[Residual]=NPSE_renew_alf2(xi,m,n,Fai,MESH,NPSE)
Ny=MESH.Ny;
y=MESH.y;
n_max=NPSE.n_max;
dx=NPSE.dx;
Zi=sqrt(-1);

  uvw_all=0.0;
  duvw_all_r=0.0;
  duvw_all_i=0.0;
  
  
	for  n3=1:Ny

uvw(n3)=abs(Fai(Ny*5-3-(n3-1)*5,1+m,n_max+1+n,xi ))^2.0+abs(Fai(Ny*5-2-(n3-1)*5,1+m,n_max+1+n,xi ))^2.0+abs(Fai(Ny*5-1-(n3-1)*5,1+m,n_max+1+n,xi ))^2.0;
               
dux(n3)=(Fai(Ny*5-3-(n3-1)*5,1+m,n_max+1+n,xi)-Fai(Ny*5-3-(n3-1)*5,1+m,n_max+1+n,xi-1))/dx;
dvx(n3)=(Fai(Ny*5-2-(n3-1)*5,1+m,n_max+1+n,xi)-Fai(Ny*5-2-(n3-1)*5,1+m,n_max+1+n,xi-1))/dx;
dwx(n3)=(Fai(Ny*5-1-(n3-1)*5,1+m,n_max+1+n,xi)-Fai(Ny*5-1-(n3-1)*5,1+m,n_max+1+n,xi-1))/dx;


duvw(n3)=conj(Fai(Ny*5-3-(n3-1)*5,1+m,n_max+1+n,xi))*dux(n3)+conj(Fai(Ny*5-2-(n3-1)*5,1+m,n_max+1+n,xi))*dvx(n3)+conj(Fai(Ny*5-1-(n3-1)*5,1+m,n_max+1+n,xi))*dwx(n3);
      

      
      if (n3>=2)
          uvw_all=uvw_all+(uvw(n3)+uvw(n3-1))*(y(n3)-y(n3-1))*0.5;
          duvw_all_r=duvw_all_r+(real(duvw(n3))+real(duvw(n3-1)))*(y(n3)-y(n3-1))*0.5;
          duvw_all_i=duvw_all_i+(imag(duvw(n3))+imag(duvw(n3-1)))*(y(n3)-y(n3-1))*0.5;
      end 

	end  
  duvw_all=duvw_all_r+Zi*duvw_all_i;
	Residual=duvw_all/uvw_all;


end