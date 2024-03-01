function dydt=NPSE_bvpeq(~,y)

parameter=NPSE_SetupParameter;
 Pr=parameter.Pr;              %普朗特数
 r=parameter.r;                  %气体常数
 Ma=parameter.Ma;                 %马赫数  
 Te=parameter.Te;              
a1=110.4/Te;
C1=sqrt(y(4))*(1+a1)/(y(4)+a1);
C11=(1+a1)*y(5)*(0.5*y(4)^(-0.5)*(y(4)+a1)-sqrt(y(4)))/((y(4)+a1)^2);%C1导数

dydt=[y(2);                                           %f的一阶导数
           y(3);                                           %f的二阶导数
           (-0.5*y(1)*y(3)-C11*y(3))/C1;    %f的三阶导数
           y(5);                                           %g的一阶导数
           (-Pr*(r-1)*Ma^2*C1*y(3)^2 -0.5*Pr*y(1)*y(5)-C11*y(5))/C1;%g的二阶导数
            y(4)];                                         %g
 %y=[f f' f'' g g' G]