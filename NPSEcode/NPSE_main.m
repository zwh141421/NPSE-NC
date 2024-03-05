clc;
clear;
close all;
addpath('../initialization')

% LST

% Flow Parameter
parameter=NPSE_SetupParameter;%1.Pr;2.r;3.Ma;4.Te;5.Rg;6.Ec;7.Re0;8.F;9.omega;10.alpha;11.beta
 
% Mesh Parameter
Ny=301;                  
Nx=101;                 
ymax=100;                
yi=3;   

li=0;                 %1:alf^2   0:no alf^2

% NPSE Parameter
NPSE=zeros(8,1);
NPSE(1) = 5;                    % m mode
NPSE(2) = 0;                    % n mode
m=1;
n=0;        
NPSE(3) = 0.00176776695;          % Amp
NPSE(4) = m*parameter(9);       % omega
NPSE(5) = n*parameter(11);       % beta
NPSE(6) = parameter(7);            % X0
NPSE(7) = 2500;                  % Xmax
NPSE(8) = (NPSE(7)-NPSE(6))/(Nx-1); % dx
X       = linspace(NPSE(6),NPSE(7),Nx);
% 
[y,z,deltaz,dy,ddy] = NPSE_Grid(Ny,yi,ymax);               
[DD1,DD2]= NPSE_dimatrix(Ny,dy,ddy,deltaz);         % difference matrix
% Base Flow
[flow0]=NPSE_Baseflow(y,Nx,Ny,X,NPSE);

% Eigenvalue Calculation
[Eval,Evec]=NPSE_eigenvalue(Ny,DD1,DD2,li,flow0,NPSE);
scatter(real(Eval),imag(Eval),'filled','LineWidth',0.1);


e_is=382;%1426 %382

N=Ny;

a1=Eval(e_is);
b1=Evec(:,e_is);


u1=zeros(1,N);
for n3=1:1:N
      u1(n3)=abs(b1(5*(n3-1)+2));
end
b1=b1/max(u1);

b1=b1(1:5*N);

den1=zeros(1,N);
u1=zeros(1,N);
v1=zeros(1,N);
w1=zeros(1,N);
T1=zeros(1,N);
for n3=1:1:N
      den1(n3)=abs(b1(5*(n3-1)+1));
      u1(n3)=abs(b1(5*(n3-1)+2));
      v1(n3)=abs(b1(5*(n3-1)+3));
      w1(n3)=abs(b1(5*(n3-1)+4));
      T1(n3)=abs(b1(5*(n3-1)+5));
end

LW=3;
figure
hold on
plot(y,u1,'LineWidth',LW,'Color',[0 0 1])
plot(y,v1,'LineWidth',LW,'LineStyle','--','Color',[1 0 0])
plot(y,w1,'LineWidth',LW,'LineStyle',':','Color',[1 0 1])
plot(y,den1,'LineWidth',LW,'LineStyle','-.','Color',[0 0.498039215803146 0])
plot(y,T1,'LineWidth',LW,'LineStyle','--','Color',[0.600000023841858 0.200000002980232 0])
hold off
axis([0 20 0 1.05])
legend1=legend('u','v','w','\rho','T');

%%
profile on
% NPSE 

[ flow,alf,Fai,modes ] = NPSE_initial(Nx,Ny,flow0,a1,b1,NPSE);

 
 m_max=NPSE(1);
 n_max=NPSE(2);


 for xi=2:10
 %for xi=2:2
    fprintf('calculating   % d / %d \n',xi,Nx);
    
    for m=0:m_max
    for n=-n_max:n_max
       if modes(1+m,n+n_max+1) 
           alf(1,(xi-1)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n)) = alf(1,(xi-2)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n));
           Fai(:,(xi-1)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n)) = Fai(:,(xi-2)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n));
       end
    end
    end
    
  tic  
% Calculate Fmn
    [Fmn]=NPSE_Fmn(xi,Ny,DD1,DD2,flow,X,alf,Fai,NPSE,modes);
    toc
    for m=0:m_max
    for n=-n_max:n_max
        fprintf('MaxFmn( % d , % d) = % e \n',m,n,max(abs(Fmn(:,m*(2*n_max+1)+(n+n_max+1))))); 
    end
    end  
  


 dv=1e-12;%Setting deviation value
 Residual=1;% deviation value
 while abs(Residual)>dv
     % 
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
                 Vyz(i*5-4:i*5,i*5-4:i*5)] =NPSE_matrix_PSE(i,xi,Ny,flow); 
     end
    for m=0:m_max
    for n=-n_max:n_max
       if modes(1+m,n+n_max+1)  
        [Residual,alf,Fai,NPSE]=NPSE_PSE(xi,m,n,X,alf,Fai,NPSE,y,Ny,DD1,DD2,Fmn,Residual,gamma,A,A2,B,C,D,Vxx,Vyy,Vzz,Vxy,Vxz,Vyz);
       end
    end
    end  
    fprintf('Residual=  % .10e ,alf=  % .10f  % .10fi \n',abs(Residual),real(alf(1,(xi-1)*(1+m_max)*(2*n_max+1)+(2*n_max+1)+(n_max+1))),imag(alf(1,(xi-1)*(1+m_max)*(2*n_max+1)+(2*n_max+1)+(n_max+1))));
 end
 
 for m=0:m_max
     for n=-n_max:n_max
         fprintf('MaxAbsFai( % d , % d) = % .10e \n',m,n,max(abs(Fai(:,(xi-1)*(1+m_max)*(2*n_max+1)+m*(2*n_max+1)+(n_max+1+n)))))
     end
 end
 
 for m=0:m_max
     for n=-n_max:n_max
         if modes(1+m,n_max+1+n)==0
         if max(abs(Fmn(:,1+m,n_max+1+n)))>1e-16
             modes(1+m,n_max+1+n)=1;
             fprintf('Mode( % d , % d) Started \n',m,n);
         end
         end
     end
 end
 fprintf('\n')
    
 end
profile viewer
