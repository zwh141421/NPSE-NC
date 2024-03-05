 function [y,z,deltaz,dy,ddy] = NPSE_Grid(Ny,yi,ymax)
 N=Ny;

 y=zeros(N,1);
 dy=zeros(N,1);
 ddy=zeros(N,1);
 z=zeros(N,1);
 
    deltaz=2/(N-1);
          for i=1:1:N                
        z(i)=-1+(i-1)*deltaz;      
          end

            a=ymax*yi/(ymax-2*yi);     
            b=1+2*a/ymax;            

           for i=1:1:N
             y(i) = a*(1+z(i))/(b-z(i));      	%Primary grid point
             dy(i) = (a+a*b)/(y(i)+a)^2;        %dz/dy
             ddy(i) = -2*(a+a*b)/(y(i)+a)^3;	%d^2z/dy^2
          end
 

 end
