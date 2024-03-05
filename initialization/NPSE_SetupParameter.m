 function [parameter] = NPSE_SetupParameter
 parameter=zeros(11,1);
 parameter(1)  = 0.72;              %普朗特数
 parameter(2)  = 1.4;                  %气体常数
 parameter(3)  = 0.01;                 %马赫数  
 parameter(4)  = 300;      
 parameter(5)  = 1/parameter(2)/parameter(3)^2;   %无量纲普适气体常数
 parameter(6)  = parameter(3)^2*( parameter(2)-1);%Ec

 parameter(7)  = 400;          %初始位置处，基于边界层厚度的雷诺数
 parameter(8)  = 86*10^(-6);
 parameter(9)  = parameter(8)*parameter(7);
 parameter(10) = 0;
 parameter(11) = 0;
 %parameter.K=-1/R;                                         %无量纲曲率
 end