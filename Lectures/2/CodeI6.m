Matlab code 6: Limit cycle 

Filename: limitcycle.m

clear;
close;


a=0.1;
b=0.5;
options=[];

[t y]=ode23('cyclefunc',[0 50],[1 1],options,a,b);
plot(y(:,1),y(:,2)); 

Filename: cyclefunc.m

function dydt = f(t,y,flag,a,b)

dydt = [-y(1)+a*y(2)+y(1)*y(1)*y(2); b-a*y(2)-y(1)*y(1)*y(2)];