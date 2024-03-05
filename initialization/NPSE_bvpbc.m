 function res = NPSE_bvpbc(y0,yinf)
    res = [y0(1);               %On the wall, f = 0
           y0(2);               %On the wall, f'= 0
           y0(5);               %On the wall, g'= 0(adiabatic condition)
           %y0(4)-3.98;         %Isothermal wall,g is the wall temperature  300K
           y0(6);               %On the wall, G = 0
           yinf(2) - 1;         %At infinity, f'= 1 
           yinf(4) - 1;];       %At infinity, g = 1
 end
