function [center, axis, theta] = calellipseparams(w)
a = w(1); 
b = w(2); 
c = w(3); 
d = w(4); 
e = w(5); 
f = w(6);

cx = (b*e-2*c*d)/(4*a*c-b^2);
cy = (b*d-2*a*e)/(4*a*c-b^2);
center = [cx,cy];
 
MA1 = sqrt(2*(a*cx^2+c*cy^2+b*cx*cy-f)/(a+c+sqrt((a-c)^2+b^2)));
MA2= sqrt(2*(a*cx^2+c*cy^2+b*cx*cy-f)/(a+c-sqrt((a-c)^2+b^2)));

axis = [max(MA1,MA2), min(MA1,MA2)];

if b==0
    if f*a>f*c
        theta = 0;
    else  
        theta = pi/2;
    end
else
    if f*a>f*c
        alpha = atan((a-c)/b);
        if alpha<0
            theta = 0.5*(-pi/2-alpha);
        else
            theta = 0.5*(pi/2-alpha);
        end
    else
        alpha = atan((a-c)/b);
        if alpha<0
            theta = pi/2+0.5*(-pi/2-alpha);
        else
            theta = pi/2+0.5*(pi/2-alpha);
        end
    end
end

theta = theta-pi/2;
end






