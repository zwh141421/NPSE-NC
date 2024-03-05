function [Eval,Evec]=NPSE_eigenvalue(Ny,DD1,DD2,li,flow0,NPSE)
%parameter=NPSE_SetupParameter;
zz=sqrt(-1);

beta  = NPSE(5);
omega = NPSE(4);
 %alpha=parameter.alpha;
 

N=Ny;

     for i=1:N
             [gamma(i*5-4:i*5,i*5-4:i*5),...
                 A(i*5-4:i*5,i*5-4:i*5),...
                 B(i*5-4:i*5,i*5-4:i*5),...
                 C(i*5-4:i*5,i*5-4:i*5),...
                 D(i*5-4:i*5,i*5-4:i*5),...
                 Vxx(i*5-4:i*5,i*5-4:i*5),...
                 Vyy(i*5-4:i*5,i*5-4:i*5),...
                 Vzz(i*5-4:i*5,i*5-4:i*5),...
                 Vxy(i*5-4:i*5,i*5-4:i*5),...
                 Vxz(i*5-4:i*5,i*5-4:i*5),...
                 Vyz(i*5-4:i*5,i*5-4:i*5)] =NPSE_matrix_baseflow(i,Ny,flow0);   
     end
 

     
   
         if li
          A1=-zz*omega*gamma+zz*beta*C+D+beta^2*Vzz;
          B1=B-zz*beta*Vyz;
          C1=-Vyy;
          M1=-zz*A-beta*Vxz;
          N1=zz*Vxy;
          P1=-Vxx;
          
          AA1=A1+B1*DD1+C1*DD2;
          AA2=-M1-N1*DD1;
          BB2=P1;
          
           AAA=[ zeros(5*N,5*N)     eye(5*N);
                       AA1                      AA2                      ;];

           BBB=[eye(5*N)              zeros(5*N,5*N);
                     zeros(5*N,5*N)            BB2                ;];
                 
         else
             
        A1=-zz*omega*gamma+zz*beta*C+D+beta^2*Vzz;
          B1=B-zz*beta*Vyz;
          C1=-Vyy;
          M1=-zz*A-beta*Vxz;
          N1=zz*Vxy;
          
          AAA=A1+B1*DD1+C1*DD2;
          BBB=M1+N1*DD1;

         end
       
       %boundary condition
       if li
              %i=1
             AAA(5*N+5,:)=AAA(5*N+1,:); BBB(5*N+5,:)=BBB(5*N+1,:);
             AAA(5*N+1:5*N+4,:)=0;  BBB(5*N+1:5*N+4,:)=0;
             AAA(5*N+1,2)=1;
             AAA(5*N+2,3)=1;
             AAA(5*N+3,4)=1;
             AAA(5*N+4,5)=1;
             %AAA(5*N+4,5+5)=-1; 
             
             %i=N
             AAA(5*N-3:5*N,:)=0;   BBB(5*N-3:5*N,:)=0;
             AAA(5*N-3,5*N-3)=1;
             AAA(5*N-2,5*N-2)=1;
             AAA(5*N-1,5*N-1)=1;
             AAA(5*N,5*N)=1;
       else
           %i=1
             AAA(5,:)=AAA(1,:); BBB(5,:)=BBB(1,:);
             AAA(1:4,:)=0;  BBB(1:4,:)=0;
             AAA(1,2)=1;
             AAA(2,3)=1;
             AAA(3,4)=1;
             AAA(4,5)=1;
             %AAA(4,5+5)=-1;   
             
             %i=N
             AAA(5*N-3:5*N,:)=0;   BBB(5*N-3:5*N,:)=0;
             AAA(5*N-3,5*N-3)=1;
             AAA(5*N-2,5*N-2)=1;
             AAA(5*N-1,5*N-1)=1;
             AAA(5*N,5*N)=1;
       end
  
   [VE,D]=eig(AAA,BBB);
    e=diag(D);
    %[~,is]=sort(abs(real(e)));  
    %es1=e(is);                 
    %VE1=VE(:,is);          

    [~,is2]=sort(imag(e));      %Sort the imaginary part of e in ascending order, where 'is' is its directory number
    Eval=e(is2);                %The eigenvalues are ordered accordingly
    Evec=VE(:,is2);             %The eigenvectors are ordered accordingly
   
end