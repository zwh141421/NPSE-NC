function [ Fmn] = NPSE_Fmn(xi,Ny,DD1,DD2,flow,X,alf,Fai,NPSE,modes)
Nz=2^5;
Nt=2^5;

m_max = NPSE(1);
n_max = NPSE(2);
omega = NPSE(4);
beta  = NPSE(5);

nonlin=zeros(Ny*Nz*Nt,5);

f   = zeros(5*Ny,(2*n_max+1)*(m_max+1));
fx  = zeros(5*Ny,(2*n_max+1)*(m_max+1));
fxx = zeros(5*Ny,(2*n_max+1)*(m_max+1));
fy  = zeros(5*Ny,(2*n_max+1)*(m_max+1));
fyy = zeros(5*Ny,(2*n_max+1)*(m_max+1));
fxy = zeros(5*Ny,(2*n_max+1)*(m_max+1));

 for m=0:m_max
 for n=-n_max:n_max
    if modes(1+m,n+n_max+1) 
     [f(:,m*(2*n_max+1)+n+n_max+1),...
      fx(:,m*(2*n_max+1)+n+n_max+1),...
      fxx(:,m*(2*n_max+1)+n+n_max+1),...
      fy(:,m*(2*n_max+1)+n+n_max+1),...
      fyy(:,m*(2*n_max+1)+n+n_max+1),...
      fxy(:,m*(2*n_max+1)+n+n_max+1)]=NPSE_Dvector(xi,m,n,DD1,DD2,Fai,NPSE);
    end
 end
 end


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
                 [F0,Ft0,Fx0,Fy0,Fz0,Fxx0,Fyy0,Fzz0,Fxy0,Fxz0,Fyz0]=NPSE_disturb(xi,m,n,sz,st,X,alf,NPSE,f,fx,fxx,fy,fyy,fxy);% \:8ba1\:7b97\:4e00\:4e2am\:ff0cn\:6a21\:6001\:7684\:6270\:52a8\:9879
                  F=F+F0;                 
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
        
        [nonlin(zi*Nt+ti+1:Nt*Nz:Ny*Nt*Nz,:)]=NPSE_physical(flow,F,Ft,Fx,Fy,Fz,Fxx,Fyy,Fzz,Fxy,Fxz,Fyz,Ny,xi);
    end  
end

% FFT
FFT=zeros(5*Ny,Nt*Nz);
for i=1:Ny
    for j=1:5
        Fft=fft2 (reshape(nonlin(Nz*Nt*(i-1)+1:Nz*Nt*i,j),[Nt,Nz]))/Nz/Nt;
        FFT(5*(i-1)+j,:)=reshape(Fft,[Nt*Nz 1]);
    end
end


Fmn=zeros(5*Ny,(m_max+1)*(2*n_max+1));
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
        
         Fmn(:,m*(2*n_max+1)+(n+n_max+1))=FFT(:,(N1-1)*Nt+M1);
          
     end
 end
 

 Fmn(1:5,:,:)=0;
 Fmn(5*Ny-4:5*Ny,:,:)=0;


end
