function ret = drawellipse(y0, x0, a, b, theta, Image, FilledColors)
% Set the elements of the Matrix Image which are in the interior of the
% ellipse E with the value 'color'. The ellipse E has center (y0, x0), the
% major axe = a, the minor axe = b, and teta is the angle macked by the
% major axe with the horizontal axe.
% ellipseMatrix(y0, x0, a, b, teta, Image, color)
%
% Function:  ellipseMatrix
% Version:   1.1
% Author:    Nicolae Cindea
% modified by Sahli Samir, Sept.4th 2009
% side conditions are not implemented
% teta = theta/180*pi;
teta = theta;
if size(FilledColors,2)==3
    col_ellipse = zeros(1,1,3);
    col_ellipse(1,1,:) = [FilledColors(1); FilledColors(2);FilledColors(3)];
else
    col_ellipse = FilledColors;
end
im = Image;
[ny, nx, D] = size(im);  % modifier pour [ny, nx, 3] 
imtemp = zeros(ny, nx,D);
% list = zeros(ny * nx, 2);  % leak of memory
list = zeros(1, 2);
toplist = 1;
c = sqrt(a * a - b * b);
x0 = round(x0);
y0 = round(y0);
list(toplist, 1) = y0;
list(toplist, 2) = x0;
im(y0, x0,:) = col_ellipse;
% correct size(im) => [:,:,3]
while (toplist > 0)
    y = list(toplist, 1);
    x = list(toplist, 2);
    toplist = toplist - 1;
    
    if local_isValid(y, x + 1, y0, x0, a, c, teta, imtemp, ny, nx, col_ellipse)
        toplist = toplist + 1;
        list(toplist, 1) = y;
        list(toplist, 2) = x + 1;
%         list(toplist, 1) = min(max(1, y), y0); % Sam, 220610
%         list(toplist, 2) = min(max(1, x+1), x0); % Sam, 220610
        im(list(toplist, 1), list(toplist, 2),:) = col_ellipse;
        imtemp(list(toplist, 1), list(toplist, 2),:) = col_ellipse;        
    end
    if local_isValid(y - 1, x, y0, x0, a, c, teta, imtemp, ny, nx, col_ellipse)
        toplist = toplist + 1;
        list(toplist, 1) = y - 1;
        list(toplist, 2) = x;
%         list(toplist, 1) = min(max(1, y-1), y0); % Sam, 220610
%         list(toplist, 2) = min(max(1, x), x0); % Sam, 220610
        im(list(toplist, 1), list(toplist, 2),:) = col_ellipse;
        imtemp(list(toplist, 1), list(toplist, 2),:) = col_ellipse;
    end
    if local_isValid(y, x - 1, y0, x0, a, c, teta, imtemp, ny, nx, col_ellipse)
        toplist = toplist + 1;
        list(toplist, 1) = y;
        list(toplist, 2) = x - 1;
%         list(toplist, 1) = min(max(1, y), y0); % Sam, 220610
%         list(toplist, 2) = min(max(1, x-1), x0); % Sam, 220610
        im(list(toplist, 1), list(toplist, 2),:) = col_ellipse;
        imtemp(list(toplist, 1), list(toplist, 2),:) = col_ellipse;        
    end
    
    % problème ICI... il faut checker apres...
    if local_isValid(y + 1, x, y0, x0, a, c, teta, imtemp, ny, nx, col_ellipse)==1
        toplist = toplist + 1;
        list(toplist, 1) = y + 1;
        list(toplist, 2) = x;
%         list(toplist, 1) = min(max(1, y+1), y0); % Sam, 220610
%         list(toplist, 2) = min(max(1, x), x0); % Sam, 220610
        im(list(toplist, 1), list(toplist, 2),:) = col_ellipse;
        imtemp(list(toplist, 1), list(toplist, 2),:) = col_ellipse;        
    end
  
end
ret = im;
%--------------------------------------------------------------------------
function is_val = local_isValid(y, x, y0, x0, a, c, teta, im, ny, nx, color)
if (y>ny || y<=0 || x>nx || x<=0)
    is_val = 0;
else
% above samir modification
    d1 = (x - x0 - c * cos(teta))^2 + (y - y0 - c * sin(teta))^2;
    d1 = sqrt(d1);
    d2 = (x - x0 + c * cos(teta))^2 + (y - y0 + c * sin(teta))^2;
    d2 = sqrt(d2);
    if (d1 + d2 <= 2*a) && (isequal( im(y, x,:), color)~=1) && (x>0) && (y>0) && ...
            (x <= nx) && (y <= ny)
        is_val = 1;
    else
        is_val = 0;
    end
end
% 
% Image = zeros(200,200,3);Image = draw_ellipse(100, 100, 200, 15, -10, Image, [1 0.5 0.25]); figure; imshow(Image);
% ??? Subscript indices must either be real positive integers or logicals.
% 
% Error in ==> draw_ellipse>local_isValid at 93
%     if (d1 + d2 <= 2*a) && (isequal( im(y, x,:), color)~=1) && (x>0) && (y>0) && ...
% 
% Error in ==> draw_ellipse at 60
%     if local_isValid(y, x - 1, y0, x0, a, c, teta, imtemp, ny, nx, col_ellipse)