% Author: Stephen Li
% Email: stephenli711@gmail.com
% Version: 2.0
% Date: July 31, 2018

% Instructions can be found in matlab-to-python v2 file.

% NOTE: Location where csv file and text file seem inconsistent. Not sure
% what the issue is.

%%
% Mldivide used on line
% Baseline for comparison against numpy.linalg.lstsq
% Answer is given by lineAns
% Results graphed in figure 1

rng default % for reproducibility
format long

% for line, y=a*x+b
clear all
close all

n = 50;
dx = 0.5;
scale = 2.0;
a = 1;
b = 0;

x = [0:dx:(n-1)*dx];
X(:,2) = x';
X(:,1) = 1;
B = [b; a];

Y = X*B+rand(n,1)*scale;

% fitting of line
Bsolve = X\Y;
Yfit = X*Bsolve;
figure, plot(x,Y, 'bo', x, Yfit, 'r-')
dlmwrite('csvlist.csv',X)
dlmwrite('csvlist.csv',Y,'-append')
lineAns = round(Bsolve,6)

fileID = fopen('resultsMatlab.txt','w');
fmt = '%.6f\r\n';
fprintf(fileID,fmt,lineAns);
fclose(fileID);

%%
% Mldivide used on curve
% Baseline for comparison against numpy.linalg.lstsq
% Answer is given by curveAns
% Results graphed in figure 2

% for 2d curve, y=a*x^2+b*x+c
clear all

n=50;
dx=0.5;
scale=20.0;
a=1;
b=2;
c=0;

x=[0:dx:(n-1)*dx];
X(:,2)=x';
X(:,3)=x'.^2;
X(:,1)=1;
B=[c; b; a];

Y=X*B+rand(n,1)*scale;

% fitting of line
Bsolve=X\Y;
Yfit=X*Bsolve;
dlmwrite('csvlist.csv',X,'-append')
dlmwrite('csvlist.csv',Y,'-append')

figure, plot(x,Y, 'bo', x, Yfit, 'r-')
curveAns = round(Bsolve,6)

fileID = fopen('resultsMatlab.txt','a');
fmt = '%.6f\r\n';
fprintf(fileID,fmt,curveAns);
fclose(fileID);

%%
% lsqlin function used
% Baseline for comparison against python's lsqlin.lsqlin function
% Answer is given by lsqlinAns
% 

% C = [0.9501    0.7620    0.6153    0.4057
%      0.2311    0.4564    0.7919    0.9354
%      0.6068    0.0185    0.9218    0.9169
%      0.4859    0.8214    0.7382    0.4102
%      0.8912    0.4447    0.1762    0.8936];
% A = [0.2027    0.2721    0.7467    0.4659
%      0.1987    0.1988    0.4450    0.4186
%      0.6037    0.0152    0.9318    0.8462];
% d = [0.0578, 0.3528, 0.8131, 0.0098, 0.1388]';
% b =[0.5251, 0.2026, 0.6721]';
C = randn(5,4);
dlmwrite('csvlist.csv',C,'-append')
A = randn(3,4);
dlmwrite('csvlist.csv',A,'-append')
d = randn(1,5)';
dlmwrite('csvlist.csv',d,'-append')
b = randn(1,3)';
dlmwrite('csvlist.csv',b,'-append')
% lb = -0.1*ones(4,1);
% ub = 2*ones(4,1);
% lsqlinAns = round(lsqlin(C, d, A, b, [], [], lb, ub),6)
lsqlinAns = round(lsqlin(C, d, A, b),6)

fileID = fopen('resultsMatlab.txt','a');
fmt = '%.6f\r\n';
fprintf(fileID,fmt,lsqlinAns);
fclose(fileID);

%%
% lsqnonneg function used
% Baseline for comparison against python's lsqlin.lsqnonneg function
% Answer is given by lsqlinAns
%

% C = [0.0372, 0.2869 
%      0.6861, 0.7071
%      0.6233, 0.6245 
%      0.6344, 0.6170];
% d = [0.8587, 0.1781, 0.0747, 0.8405]';
C = randn(4,2);
dlmwrite('csvlist.csv',C,'-append')
d = randn(1,4)';
dlmwrite('csvlist.csv',d,'-append')
lsqnonnegAns = round(lsqnonneg(C, d),6)

fileID = fopen('resultsMatlab.txt','a');
fmt = '%.6f\r\n';
fprintf(fileID,fmt,lsqnonnegAns);
fclose(fileID);

%%
% %lsqcurvefit
% xdata = ...
%   [0.0 1.0 2.0 3.0 4.0 5.0]
% ydata = ...
%   [0.1 0.9 2.2 2.8 3.9 5.1]
% fun = @(x,xdata)x(1)+x(2)*xdata+x(3)*xdata*xdata;
% x0 = [0.0, 0.0, 0.0];
% %x = lsqcurvefit(fun,x0,xdata,ydata)

% xdata = ...
%  [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
% ydata = ...
%  [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5];
% fun = @(x,xdata)x(1)*exp(x(2)*xdata);
% x0 = [100,-1];
% x = lsqcurvefit(fun,x0,xdata,ydata)
% xdata
% ydata
% fun(x0,xdata)

%%
% lsqnonlin function used
% Baseline for comparison against scipy.optimize.least_squares function
% Answer is given by lsqnonlinAns
% 

d = linspace(0,3);
dlmwrite('csvlist.csv',d,'-append')
y = exp(-1.3*d) + 0.05*randn(size(d));
dlmwrite('csvlist.csv',y,'-append')
fun = @(r)exp(-d*r)-y;
x0 = 4;
lsqnonlinAns = round(lsqnonlin(fun,x0),6)

fileID = fopen('resultsMatlab.txt','a');
fmt = '%.6f\r\n';
fprintf(fileID,fmt,lsqnonlinAns);
fclose(fileID);