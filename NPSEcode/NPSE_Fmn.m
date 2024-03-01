function [ Fmn] = NPSE_Fmn(xi,MESH,flow,NPSE,modes)
Nz=2^5;
Nt=2^5;%z、t坐标上划分网格数
Ny=MESH.Ny;
m_max=NPSE.m_max;
n_max=NPSE.n_max;
omega=NPSE.omega;
beta=NPSE.beta;

nonlin=zeros(Ny,5,Nt,Nz);   %存储物理空间非线性项N 

if omega==0
    Tt=0;
else 
    Tt=2*pi/omega;
end

if beta==0
    Tz=0;
else
    Tz=2*pi/beta;
end

dz=Tz/Nz;
dt=Tt/Nt;

for zi=0:Nz-1
    for ti=0:Nt-1
        sz=dz*zi;
        st=dt*ti;
        
       %用以存储扰动项
        F=zeros(5*Ny,1);
        Ft=zeros(5*Ny,1);
        Fx=zeros(5*Ny,1);
        Fy=zeros(5*Ny,1);
        Fz=zeros(5*Ny,1);
        Fxx=zeros(5*Ny,1);
        Fyy=zeros(5*Ny,1);
        Fzz=zeros(5*Ny,1);
        Fxy=zeros(5*Ny,1);
        Fxz=zeros(5*Ny,1);
        Fyz=zeros(5*Ny,1);
        for m=0:m_max
            for n=-n_max:n_max
              if modes(1+m,n+n_max+1) 
                  [F0,Ft0,Fx0,Fy0,Fz0,Fxx0,Fyy0,Fzz0,Fxy0,Fxz0,Fyz0]=NPSE_distrub(xi,m,n,sz,st,MESH,NPSE);%计算一个m，n模态的扰动项
                  F=F+F0;                 %各模态扰动求和
                  Ft=Ft+Ft0;
                  Fx=Fx+Fx0;
                  Fy=Fy+Fy0;
                  Fz=Fz+Fz0;
                  Fxx=Fxx+Fxx0;
                  Fyy=Fyy+Fyy0;
                  Fzz=Fzz+Fzz0;
                  Fxy=Fxy+Fxy0;
                  Fxz=Fxz+Fxz0;
                  Fyz=Fyz+Fyz0;
              end
            end
        end
        
        [nonlin(:,:,ti+1,zi+1)]=NPSE_physical(flow,F,Ft,Fx,Fy,Fz,Fxx,Fyy,Fzz,Fxy,Fxz,Fyz,MESH,xi);                               %计算物理空间非线性项N
        %no=nonlin(:,:,1,1);
    end  
end

%FFT计算
FFT=zeros(5*Ny,Nt,Nz);
for i=1:Ny
    for j=1:5
        Fft=fft2(reshape(nonlin(i,j,:,:),[Nt,Nz]))/Nz/Nt;
        FFT(5*(i-1)+j,:,:)=reshape(Fft,1,Nt,Nz);
    end
end


%从FFT计算结果中取出需要的
Fmn=zeros(5*Ny,m_max+1,2*n_max+1);
 for n=-n_max:n_max
     for m=0:m_max
        if(n>0)
            N1=n+1;
        elseif(n==0)
            N1=1;
        else
            N1=Nz+n+1;
        end

        if(m>0)
            M1=Nt-(m-1);
        else
            M1=1;
        end
        
         Fmn(:,m+1,n+n_max+1)=FFT(:,M1,N1);
          
     end
 end
 
 %for i=1:32
   %  for j=1:32
     %        fprintf('MaxFmn( %d , %d) = %e\n',i,j,max(abs(FFT(:,i,j)))); 
    % end
% end

 Fmn(1:5,:,:)=0;
 Fmn(5*Ny-4:5*Ny,:,:)=0;


end