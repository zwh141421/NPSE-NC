function [ MESH ] = NPSE_dimatrix(MESH)

N=MESH.Ny;
dy=MESH.dy;
ddy=MESH.ddy;
deltaz=MESH.deltaz;

%初始化离散矩阵
D1=zeros(N,N);
D2=zeros(N,N);

     for i=3:N-2
          D1(i,i-2:i+2)  = [1/12  -8/12  0  8/12  -1/12];                                     
          D2(i,i-2:i+2) = [-1/12  16/12 -30/12 16/12 -1/12];    
     end

       i=2;
           D1(i,i-1:i+3) = [-3/12  -10/12  18/12  -6/12  1/12];
           D2(i,i-1:i+4) = [10/12  -15/12  -4/12  14/12  -6/12  1/12];

       i=N-1;
          D1(i,i-3:i+1) = [-1/12  6/12  -18/12  10/12  3/12];
          D2(i,i-4:i+1) = [1/12  -6/12  14/12  -4/12  -15/12  10/12];

       i=1;
         D1(i,i:i+4) = [-25/12  48/12  -36/12  16/12  -3/12];

       i=N;
         D1(i,i-4:i) = [3/12  -16/12  36/12  -48/12  25/12];

        


  D11= D1.*dy/deltaz;
  D22= D2.*dy.^2/deltaz^2+D1.*ddy/deltaz; 
  MESH.D11=D11;
  MESH.D22=D22;

end