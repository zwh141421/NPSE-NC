function[Residual,alf,Fai,NPSE]=NPSE_PSE(xi,m,n,X,alf,Fai,NPSE,y,Ny,DD1,DD2,Fmn,Residual,gamma,A,A2,B,C,D,Vxx,Vyy,Vzz,Vxy,Vxz,Vyz)

 m_max = NPSE(1);
 n_max = NPSE(2);
 omega = NPSE(4);
 beta  = NPSE(5);
 dx    = NPSE(8);

 x     = X;
 zz    = sqrt(-1);
 alfmn = alf(1,(xi-1)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n));
 alfx  = 0;
  
  if m==1&&n==0
      
  else
      alf(1,(xi-1)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n))=m/1*real(alf(1,(xi-1)*(1+m_max)*(2*n_max+1)+1*(2*n_max+1)+(n_max+1+0)));
  end

  JF=trapz(x(1:xi),alf(1,m*(2*n_max+1)+(n_max+1+n):(1+m_max)*(2*n_max+1):(xi-1)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n)));


     
     A1=A2-2*zz*Vxx*alfmn-zz*n*beta*Vxz;
     B1=B-zz*alfmn*Vxy-zz*n*beta*Vyz;
     D1=D+zz*n*beta*C-zz*m*omega*gamma+zz*alfmn*A+n^2*beta^2*Vzz-zz*alfx*Vxx+alfmn^2*Vxx+alfmn*n*beta*Vxz;
     C1=Vyy;
     
     
     % AAA=A1/dx;
     % BBB=C1*DD2-B1*DD1-D1+A1/dx;
     AAA=A1/dx+B1*DD1+D1-C1*DD2;
     BBB=A1/dx;
     
     % Boundary condition
     AAA(5,:)=AAA(1,:);
     AAA(1:4,:)=0; AAA(5*Ny-3:5*Ny,:)=0;
     BBB(5,:)=BBB(1,:);
     BBB(1:4,:)=0; BBB(5*Ny-3:5*Ny,:)=0;
     
    
     AAA(1,2)=1;% u=0
     AAA(2,3)=1;% v=0
     AAA(3,4)=1;% w=0
     AAA(4,5)=1;% T=0
     
     AAA(5*Ny-3,5*Ny-3)=1;% u=0
     AAA(5*Ny-1,5*Ny-1)=1;% w=0
     AAA(5*Ny,5*Ny)=1;    % T=0
     
     if m==0&&n==0
         AAA(5*Ny-2,5*Ny-2)=1.5;
         AAA(5*Ny-2,5*Ny-7)=-2;
         AAA(5*Ny-2,5*Ny-12)=0.5;
     else
         AAA(5*Ny-2,5*Ny-2)=1;% v=0
         
     end
     

     
     BBB2=(BBB*Fai(:,(xi-2)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n))+Fmn (:,m*(2*n_max+1)+(n+n_max+1))/exp(zz*JF));
     Fai(:,(xi-1)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n))=AAA\BBB2;
     

     
     % Iterative alpha
     if m==1&&n==0
     [Residual]=NPSE_renew_alf(xi,m,n,Fai,y,Ny,NPSE);
     alf(1,(xi-1)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n))=alf(1,(xi-1)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n))-Residual*zz;
     end

end
