function dydt=NPSE_bvpeq(~,y)

parameter=NPSE_SetupParameter;
 Pr = parameter(1);              
 r  = parameter(2);              
 Ma = parameter(3);               
 Te = parameter(4);              
a1=110.4/Te;
C1=sqrt(y(4))*(1+a1)/(y(4)+a1);
C11=(1+a1)*y(5)*(0.5*y(4)^(-0.5)*(y(4)+a1)-sqrt(y(4)))/((y(4)+a1)^2);%first derivative of C1

dydt=[y(2);                                           %first derivative of f
           y(3);                                           %second derivative of f
           (-0.5*y(1)*y(3)-C11*y(3))/C1;    %third derivative of f
           y(5);                                           %first derivative of g
           (-Pr*(r-1)*Ma^2*C1*y(3)^2 -0.5*Pr*y(1)*y(5)-C11*y(5))/C1;%second derivative of g
            y(4)];                                         %g
 %y=[f f' f'' g g' G]