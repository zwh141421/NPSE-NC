function[Residual,NPSE]=NPSE_PSE(xi,m,n,NPSE,MESH,flow,Fmn,Residual)

 zi=sqrt(-1);
 beta=NPSE.beta;
 omega=NPSE.omega;
 m_max=NPSE.m_max;
 n_max=NPSE.n_max;
 dx=NPSE.dx;
 x=NPSE.X;
  alf=NPSE.alf;
  Fai=NPSE.Fai;
  alfmn=alf(1+m,n_max+1+n,xi);
  alfx=0;
  
  if m==1&&n==0
      
  else
      alf(1+m,n_max+1+n,xi)=m/1*real(alf(1+1,n_max+1+0,xi));
  end
  JF=trapz(x(1:xi),alf(1+m,n_max+1+n,1:xi));


D11=MESH.D11;
D22=MESH.D22;
Ny=MESH.Ny;
 %参数矩阵
     for i=1:Ny
             [gamma(i*5-4:i*5,i*5-4:i*5),...
                 A(i*5-4:i*5,i*5-4:i*5),...
                 A2(i*5-4:i*5,i*5-4:i*5),...
                 B(i*5-4:i*5,i*5-4:i*5),...
                 C(i*5-4:i*5,i*5-4:i*5),...
                 D(i*5-4:i*5,i*5-4:i*5),...
                 Vxx(i*5-4:i*5,i*5-4:i*5),...
                 Vyy(i*5-4:i*5,i*5-4:i*5),...
                 Vzz(i*5-4:i*5,i*5-4:i*5),...
                 Vxy(i*5-4:i*5,i*5-4:i*5),...
                 Vxz(i*5-4:i*5,i*5-4:i*5),...
                 Vyz(i*5-4:i*5,i*5-4:i*5)] =NPSE_matrix_PSE(i,xi,MESH,flow);   %%%%%%1.6s
     end
 
     DD1=kron(D11,eye(5,5));  DD2=kron(D22,eye(5,5)); 
     
     A1=A2-2*zi*Vxx*alfmn-zi*n*beta*Vxz;
     B1=B-zi*alfmn*Vxy-zi*n*beta*Vyz;
     D1=D+zi*n*beta*C-zi*m*omega*gamma+zi*alfmn*A+n^2*beta^2*Vzz-zi*alfx*Vxx+alfmn^2*Vxx+alfmn*n*beta*Vxz;
     C1=Vyy;
     
     %求解φ
     %AAA=A1/dx;
     %BBB=C1*DD2-B1*DD1-D1+A1/dx;
     AAA=A1/dx+B1*DD1+D1-C1*DD2;
     BBB=A1/dx;
     
     %边界条件
     AAA(5,:)=AAA(1,:);
     AAA(1:4,:)=0; AAA(5*Ny-3:5*Ny,:)=0;
     BBB(5,:)=BBB(1,:);
     BBB(1:4,:)=0; BBB(5*Ny-3:5*Ny,:)=0;
     
    
     AAA(1,2)=1;%u=0
     AAA(2,3)=1;%v=0
     AAA(3,4)=1;%w=0
     AAA(4,5)=1;%T=0
     
     AAA(5*Ny-3,5*Ny-3)=1;%u=0
     AAA(5*Ny-1,5*Ny-1)=1;%w=0
     AAA(5*Ny,5*Ny)=1;    %T=0
     
     if m==0&&n==0
         AAA(5*Ny-2,5*Ny-2)=1.5;
         AAA(5*Ny-2,5*Ny-7)=-2;
         AAA(5*Ny-2,5*Ny-12)=0.5;
     else
         AAA(5*Ny-2,5*Ny-2)=1;%v=0
         
     end
     
     %a(:,1)=Fai(:,1+m,n_max+1+n,xi-1);
     %a(:,2)=Fai(:,1+m,n_max+1+n,xi);
     
     BBB2=(BBB*Fai(:,1+m,n_max+1+n,xi-1)+Fmn(:,1+m,n_max+n+1)/exp(zi*JF));
     Fai(:,1+m,n_max+1+n,xi)=AAA\BBB2;
     
     %b(:,1)=Fai(:,1+m,n_max+1+n,xi-1);
     %b(:,2)=Fai(:,1+m,n_max+1+n,xi);
     
     %更新α的值
     if m==1&&n==0
     [Residual]=NPSE_renew_alf(xi,m,n,Fai,MESH,NPSE);
     alf(1+m,n_max+1+n,xi)=alf(1+m,n_max+1+n,xi)-Residual*zi;
     end
    % Residual
     NPSE.Fai=Fai;
     NPSE.alf=alf;
end