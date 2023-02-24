function [alpha, b] = hyperfitellipse(V)
N = size(V,1);
b = max(abs(V(:)))/255;
x = V(:,1);
y = V(:,2);
Chi = [x.^2,y.^2,2*b*x,2*b*y,b^2*ones(N,1)];


X = mean(x);
Y = mean(y);
XY = mean(Chi(:,2))/2;
X2 = mean(Chi(:,1));
Y2 = mean(Chi(:,3));
W = [6*X2  X2+Y2 6*b*X 2*b*Y  b^2;...
     X2+Y2 6*Y2  2*b*X  6*b*Y b^2;...
     6*b*X  2*b*X 4*b^2  0     0 ;...
     2*b*Y  6*b*Y 0    4*b^2  0 ;...
     b^2    b^2   0     0     0];
 opt.issym = true;
 X = Chi'*Chi/N;
 [alpha,~] = eigs(W,X,1,'lm',opt);
end