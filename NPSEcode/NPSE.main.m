%LST前处理计算部分
global parameter
%气体基本参数
 parameter.Pr=0.72;              %普朗特数
 parameter. r=1.4;                  %气体常数
 parameter. Ma=0.01;                 %马赫数  
 parameter. Te=300;      
 parameter.Rg=1/parameter.r/parameter.Ma^2;   %无量纲普适气体常数
 parameter.Ec= parameter. Ma^2*( parameter. r-1);

 parameter.Re0=400;          %初始位置处，基于边界层厚度的雷诺数
 parameter.F=86*10^(-6);
 parameter.omega=parameter.F*parameter.Re0;
 parameter.alpha=0;
 parameter.beta=0;
 %parameter.K=-1/R;                                         %无量纲曲率
 
%计算用网格参数
MESH.Ny=301;                   %法向网格点数
MESH.Nx=101;                  %流向网格数
MESH.ymax=100;                 %网格点最大值
MESH.yi=3;                    %网格加密点范围参数
MESH.li=0;                 %1:考虑特征值平方项。   0：不考虑特征值平方项

%NPSE计算参数
NPSE.m_max=5;                    %m模态数
NPSE.n_max=0;                    %n模态数
m=1;
n=0;        
NPSE.Amp=0.00176776695;          %扰动的幅值/根号2
NPSE.omega=m*parameter.omega;
NPSE.beta=n*parameter.beta;
NPSE.X0=parameter.Re0;
NPSE.Xmax=2500;
NPSE.X=linspace(NPSE.X0,NPSE.Xmax,MESH.Nx);
NPSE.dx=(NPSE.Xmax-NPSE.X0)/(MESH.Nx-1);
%离散方法
MESH = NPSE_Grid(MESH);               %画网格
MESH = NPSE_dimatrix(MESH);         %微分矩阵
%基本流
[flow]=NPSE_Baseflow(MESH,NPSE);

%初始α和φ
[Eval,Evec]=NPSE_eigenvalue(MESH,flow,NPSE);
scatter(real(Eval),imag(Eval),'filled','LineWidth',0.1);%特征值

%绘图
e_is=382;%1426 %382

N=MESH.Ny;
y=MESH.y;
a1=Eval(e_is);
b1=Evec(:,e_is);


u1=zeros(1,N);
for n3=1:1:N
      u1(n3)=abs(b1(5*(n3-1)+2));
end
b1=b1/max(u1);

NPSE.a1=a1;
NPSE.b1=b1(1:5*N);

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

%NPSE计算部分
%x方向第一处位置处的初值
[ alf,Fai,modes ] = NPSE_initial(MESH,NPSE);
 NPSE.alf=alf;
 NPSE.Fai=Fai;
 
 m_max=NPSE.m_max;
 n_max=NPSE.n_max;

%load matlab2.mat
for xi=2:MESH.Nx
 %for xi=2:20
    fprintf('calculating   %d / %d\n',xi,MESH.Nx);
    
    for m=0:m_max
        for n=-n_max:n_max
           if modes(1+m,n+n_max+1) 
               NPSE.alf(1+m,n_max+1+n,xi)=NPSE.alf(1+m,n_max+1+n,xi-1);
               NPSE.Fai(:,1+m,NPSE.n_max+1+n,xi)= NPSE.Fai(:,1+m,n_max+1+n,xi-1);
           end
        end
    end
    
    
%计算Fmn
    [Fmn]=NPSE_Fmn(xi,MESH,flow,NPSE,modes);
    
    for m=0:m_max
        for n=-n_max:n_max
         fprintf('MaxFmn( %d , %d) = %e\n',m,n,max(abs(Fmn(:,1+m,n_max+n+1)))); 
        end
    end  
  

 %load matlab3.mat  
 dv=1e-12;%指定差值标准
 Residual=1;%实际残差
 while abs(Residual)>dv
    for m=0:m_max
        for n=-n_max:n_max
           if modes(1+m,n+n_max+1)  
           [Residual,NPSE]=NPSE_PSE(xi,m,n,NPSE,MESH,flow,Fmn,Residual);
           end
        end
    end  
    fprintf('Residual=  %.10e ,alf=  %.10f  %.10fi\n',abs(Residual),real(NPSE.alf(1+1,1,xi)),imag(NPSE.alf(1+1,1,xi)));
 end
 
 for m=0:m_max
     for n=-n_max:n_max
         fprintf('MaxAbsFai( %d , %d) = %.10e\n',m,n,max(abs(NPSE.Fai(:,1+m,n_max+1+n,xi))))
     end
 end
 
 for m=0:m_max
     for n=-n_max:n_max
         if modes(1+m,n_max+1+n)==0
         if max(abs(Fmn(:,1+m,n_max+1+n)))>1e-16
             modes(1+m,n_max+1+n)=1;
             fprintf('Mode( %d , %d) Started\n',m,n);
         end
         end
     end
 end
 fprintf('\n')
    
 end
