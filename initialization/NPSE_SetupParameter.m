 function [parameter] = NPSE_SetupParameter
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
 end