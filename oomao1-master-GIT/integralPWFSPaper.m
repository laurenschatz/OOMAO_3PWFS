y=-2*pi:0.1:2*pi;
f=@(x,y) (x.^2-y.^2)./(x-y);

g=integral(@(x)f(x,y), -1, 1);






 



