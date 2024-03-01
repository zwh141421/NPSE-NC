function [Eval,Evec]=NPSE_eigenvalue(MESH,flow,NPSE)
parameter=NPSE_SetupParameter;
zz=sqrt(-1);

 beta=NPSE.beta;
 omega=NPSE.omega;
 alpha=parameter.alpha;
 
D11=MESH.D11;
D22=MESH.D22;
N=MESH.Ny;

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
                 Vyz(i*5-4:i*5,i*5-4:i*5)] =NPSE_matrix_baseflow(i,MESH,flow);   
     end
 
     DD1=kron(D11,eye(5,5));  DD2=kron(D22,eye(5,5)); 
     
   
         if MESH.li
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
       
       %设定边界条件
       if MESH.li
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
    %[~,is]=sort(abs(real(e)));  %将e的实部按升序排列，is是其目录编号
    %es1=e(is);                         %将特征值按照实部排列
    %VE1=VE(:,is);                    %特征向量按照特征值的顺序

    [~,is2]=sort(imag(e));      %将e的虚部按升序排列，is是其目录编号
    Eval=e(is2);                     %特征值相应排序
    Evec=VE(:,is2);                 %特征向量相应排序
   
end