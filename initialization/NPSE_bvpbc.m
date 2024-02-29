 function res = NPSE_bvpbc(y0,yinf)
    res = [y0(1);               %在壁面，f为零
           y0(2);               %在壁面，f的一阶导数为零
           y0(5);               %在壁面，g的一阶导数为零 绝热壁面
           %y0(4)-3.98;          %等温壁面，g为壁温300K
           y0(6);                %在壁面，G为零
           yinf(2) - 1;         %在远场，f的一阶导数为1
           yinf(4) - 1;];       %在远场，g为1
 end