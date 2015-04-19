%Display power utility functions with different degrees of risk aversion.
clear all

C = linspace(0.8,1.2,100);
Eta0 = 0;
Eta2 = 2;
Eta5 = 5;
Eta8 = 8;

U0 = (C.^(1-Eta0)-1)/(1-Eta0);

U1 = log(C);

U2 = (C.^(1-Eta2)-1)/(1-Eta2);

U5 = (C.^(1-Eta5)-1)/(1-Eta5);

U8 = (C.^(1-Eta8)-1)/(1-Eta8);

figure
plot(C,U0,C,U1,C,U2,C,U5,C,U8);
legend('eta=0','eta=1','eta=2','eta=5','eta=8','Location','SouthEast');
xlim([0.8,1.2])<center>页面执行时间：15.625 毫秒</center>