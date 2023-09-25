function z = f_3piece(r0,x) 
% r is the Catersian coordinates of 2 points that define the 3 piece
% function (along with the origina & slope = 0 at origin)

r = r0;

%% middle parabola

y2 = r(1)*(x-r(2)).^2 + r(3);

%% calculating alpha_n and alpha_p
% left piece
y1 = r(5)*(x+r(4)-r(2)) + r(1)*(-r(4))^2 + r(3);

% right piece
y3 = -r(5)*(x-r(4)-r(2)) + r(1)*(r(4))^2 + r(3);

%% whole curve
z = y1 .* (x<-r(4)+r(2)) + y2 .* (x>=-r(4)+r(2) & x<r(4)+r(2)) + y3 .* (x>=r(4)+r(2));

end
