function [C1,C2] = lineEllipse(a,b,O,A,B)
% Get points of intersection of line and Ellipse
% INPUT : a - major axis of Ellipse
%         b - minor axis of Ellipse
%         O - center of ellipse i,e (h,k)
%         A, P  - points of striaght line (x1,y1) and (x2,y2)
% OUTPUT : C - Point of intersection of line and ellipse

% Reference : http://www.ambrsoft.com/TrigoCalc/Circles2/Ellipse/EllipseLine.htm
% 
% Coded by :    Siva Srinivas Kolukula, PhD      
%               Tsunami and Storm surge Early WArning Group (TWS)
%               Indian National Centre for Ocean Information Services (INCOIS)
%               Hyderabad, INDIA
% E-mail   :    allwayzitzme@gmail.com                                        
% web-link :    https://sites.google.com/site/kolukulasivasrinivas/   
%% Ellipse center
h = O(1) ; k = O(2) ;
%% Line, Get slope and y-intercept of the line
AB = [A; B];
line = polyfit(AB(:,1),AB(:,2),1) ;
m = line(1) ;       % Slope of the line 
c = line(2) ;       % y-intercept of the line
%% Formula
eps = c-k ;delta = c+m*h ;
D = sqrt(a^2*m^2+b^2-delta^2-k^2+2*delta*k) ;
E = a^2*m^2+b^2 ;
x1 = (h*b^2-m*a^2*eps+a*b*D)/E ;x2 = (h*b^2-m*a^2*eps-a*b*D)/E ;
y1 = (b^2*delta+k*a^2*m^2+a*b*m*D)/E ;y2 = (b^2*delta+k*a^2*m^2-a*b*m*D)/E ;
C1 = [x1 y1] ; C2 = [x2 y2] ;
%% Check if intersection points exists
if ~isreal(C1)
    C1 = [NaN NaN] ;
    C2 = C1 ;
end
