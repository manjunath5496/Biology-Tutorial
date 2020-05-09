Matlab code 2: Lambda genetic regulation kinetics

Filename: hasty.m

alpha=50;
gamma=20;
sigma1=1;
sigma2=5;
options=[];

[t1 y1]=ode23('hastyfunc',[0 10],[0],options,alpha,gamma,sigma1,sigma2);
[t2 y2]=ode23('hastyfunc',[0 10],[1],options,alpha,gamma,sigma1,sigma2);
plot(t1,y1(:,1),'b',t2,y2(:,1),'r');

Filename: hastyfunc.m 

function dydt = f(t,y,flag,alpha,gamma,sigma1,sigma2)
% [x] = y(1)


dydt = [alpha*y(1)^2/(1+(1+sigma1)*y(1)^2+sigma2*y(1)^4)-gamma*y(1)+1];