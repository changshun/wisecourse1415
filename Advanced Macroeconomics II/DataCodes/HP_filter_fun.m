function [xhp,xtr]= varargin(x,HP_LAMBDA);
% If x is a column vector of length LENGTH
% xtr=HP_mat\x; delivers the HP-trend and
% xhp=x-xtr; delivers the HP-filtered series
% This program computes HP_mat, given the length of
% some given column vector x
%disp('If x is a column vector of length LENGTH,');
%disp('xtr=HP_mat\x; delivers the HP-trend and');
%disp('xhp=x-xtr; delivers the HP-filtered series');
%disp('This program computes HP_mat, given the length of');
%disp('some given column vector x');
LENGTH = max(size(x));
% The following piece is due to Gerard A. Pfann
   HP_mat = [1+HP_LAMBDA, -2*HP_LAMBDA, HP_LAMBDA,              zeros(1,LENGTH-3);
             -2*HP_LAMBDA,1+5*HP_LAMBDA,-4*HP_LAMBDA,HP_LAMBDA, zeros(1,LENGTH-4);
                           zeros(LENGTH-4,LENGTH);
              zeros(1,LENGTH-4),HP_LAMBDA,-4*HP_LAMBDA,1+5*HP_LAMBDA,-2*HP_LAMBDA;                   zeros(1,LENGTH-3),          HP_LAMBDA,   -2*HP_LAMBDA, 1+HP_LAMBDA  ];
   for iiiii=3:LENGTH-2;
     HP_mat(iiiii,iiiii-2)=HP_LAMBDA;
     HP_mat(iiiii,iiiii-1)=-4*HP_LAMBDA;
     HP_mat(iiiii,iiiii)=1+6*HP_LAMBDA;
     HP_mat(iiiii,iiiii+1)=-4*HP_LAMBDA;
     HP_mat(iiiii,iiiii+2)=HP_LAMBDA;
   end;
xtr=HP_mat\x;
xhp=x-xtr; 
      
<center>页面执行时间：0 毫秒</center>