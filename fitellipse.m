function w = fitellipse(x,y)
D = [x.*x, x.*y, y.*y, x, y,ones(size(x))];
S = D'*D;
G = zeros(6);
G(1,3) = 2; G(3,1) = 2; G(2,2) = -1;

[vec, val] = eig(S\G);
[~, idx] = find(val>0&~isinf(val));
w = vec(:,idx);
w = sqrt(1/(w'*S*w))*w;
end



