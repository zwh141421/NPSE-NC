 function [MESH] = NPSE_Grid(MESH)
 N=MESH.Ny;
 yi=MESH.yi;
 ymax=MESH.ymax;
 y=zeros(N,1);
 dy=zeros(N,1);
 ddy=zeros(N,1);
 z=zeros(N,1);
 
    deltaz=2/(N-1);
          for i=1:1:N                
        z(i)=-1+(i-1)*deltaz;      %在[-1,1]等分
          end

            a=ymax*yi/(ymax-2*yi);      %映射变换参数
            b=1+2*a/ymax;               %映射变换参数

           for i=1:1:N
             y(i) = a*(1+z(i))/(b-z(i));      	%按映射关系计算原网格点
             dy(i) = (a+a*b)/(y(i)+a)^2;       %映射坐标z关于y的一阶导数
             ddy(i) = -2*(a+a*b)/(y(i)+a)^3;	%映射坐标z关于y的二阶导数
          end
 
  MESH.y=y;   
  MESH.z=z;
  MESH.deltaz=deltaz;
  MESH.dy=dy;
  MESH.ddy=ddy;
 end
