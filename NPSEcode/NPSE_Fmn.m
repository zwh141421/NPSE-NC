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
        
       %用以存储扰动项的结构体
        distrub.F=zeros(5*Ny,1);
        distrub.Ft=zeros(5*Ny,1);
        distrub.Fx=zeros(5*Ny,1);
        distrub.Fy=zeros(5*Ny,1);
        distrub.Fz=zeros(5*Ny,1);
        distrub.Fxx=zeros(5*Ny,1);
        distrub.Fyy=zeros(5*Ny,1);
        distrub.Fzz=zeros(5*Ny,1);
        distrub.Fxy=zeros(5*Ny,1);
        distrub.Fxz=zeros(5*Ny,1);
        distrub.Fyz=zeros(5*Ny,1);
        for m=0:m_max
            for n=-n_max:n_max
              if modes(1+m,n+n_max+1) 
                  [distrub0]=NPSE_distrub(xi,m,n,sz,st,MESH,NPSE);%计算一个m，n模态的扰动项
                  distrub.F=distrub.F+distrub0.F;                 %各模态扰动求和
                  distrub.Ft=distrub.Ft+distrub0.Ft;
                  distrub.Fx=distrub.Fx+distrub0.Fx;
                  distrub.Fy=distrub.Fy+distrub0.Fy;
                  distrub.Fz=distrub.Fz+distrub0.Fz;
                  distrub.Fxx=distrub.Fxx+distrub0.Fxx;
                  distrub.Fyy=distrub.Fyy+distrub0.Fyy;
                  distrub.Fzz=distrub.Fzz+distrub0.Fzz;
                  distrub.Fxy=distrub.Fxy+distrub0.Fxy;
                  distrub.Fxz=distrub.Fxz+distrub0.Fxz;
                  distrub.Fyz=distrub.Fyz+distrub0.Fyz;
              end
            end
        end
        
        [nonlin(:,:,ti+1,zi+1)]=NPSE_physical(flow,distrub,MESH,xi);                               %计算物理空间非线性项N
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