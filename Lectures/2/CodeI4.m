Matlab code 4: Spiro's model for chemotaxis

Filename: spiro.m

clear;
close;


To=8e-6;
Yo=20e-6;
Bo=1.7e-6;


options = odeset('RelTol',1e-9,'AbsTol',[1e-9 1e-9 1e-9 1e-9 1e-9]);
[t y]=ode23('spirofunc',[0 80],[4e-6 4e-6 0e-6 0.5e-6 10e-6],options);


Ptot=To-y(:,1)-y(:,2);
Bp=y(:,4);
Yp=y(:,5);
metlevel=1-(y(:,1)+y(:,3))/To;
phoslevel=1-(y(:,1)+y(:,2))/To;


subplot(2,2,1) plot(t,phoslevel,'bx');
axis([10 80 0 0.1]);
title('Phosphorylation Level');
subplot(2,2,2) plot(t,metlevel,'rx');
title('Methylation Level');
axis([10 80 0 1]);
subplot(2,2,3) plot(t,Bp/Bo,'gx');
axis([10 80 0 1]);
title('Bp/Btot');
subplot(2,2,4) plot(t,Yp/Yo,'yx');
axis([10 80 0 1]);
title('Yp/Ytot'); 

Filename: spirofunc.m


function dydt = f(t,y,flag)
% constants from Table 3 (Spiro et al.)
k1c=0.17; % 1/s
k3c=30*k1c; % 1/s
ratiok1bk1a=1.7e-6; % M
ratiok3ck3a=1.7e-6; % M
k_1=4e5; % 1/(Ms)
k_3=k_1; % 1/(Ms)|
k8=15; % 1/s
k9=3*k8; % 1/s
k11=0; % 1/s
k12=1.1*k8; % 1/s 
kb=8e5; % 1/(Ms) 
ky=3e7; % 1/(Ms) 
k_b=0.35; % 1/s 
k_y=5e5; % 1/(Ms) 
Kbind=1e6; % 1/M 


Yo=20e-6; % M 
Bo=1.7e-6; % M 
To=8e-6; % M 
Ro=0.3e-6; % M 
Zo=40e-6; % M 

% [T2]+[LT2] = y(1) 
% [T3]+[LT3] = y(2) 
% [T2p]+[LT2p] = y(3) 
% [Bp] = y(4) 
% [Yp] = y(5) 


cligand=1e-6; if t>20 cligand=1e-3; end; 
if t>50 cligand=1e-6; end; 
Vmaxunbound=k1c*Ro; % maximum turnover rate (MM kinetics) for unbound receptors
Vmaxbound=k3c*Ro; % maximum turnover rate (MM kinetics) for bound receptors 
KR=ratiok1bk1a; % Michaelis constant 
fb=Kbind*cligand/(1+Kbind*cligand); % fraction receptors bound to ligand
fu=1-fb; % fraction receptors not bound to ligand
kpt=ky*(Yo-y(5))+kb*(Bo-y(4));


ydot1=(-k8*fu-k11*fb)*y(1)+kpt*y(3)+(k_1*fu+k_3*fb)*y(2)*y(4)-Vmaxunbound*y(1)*fu/(KR+y(1)*fu)-Vmaxbound*y(1)*fb/(KR+y(1)*fb);
ydot2=(-k9*fu-k12*fb)*y(2)+kpt*(To-y(1)-y(2)-y(3))-(k_1*fu+k_3*fb)*y(2)*y(4)+Vmaxunbound*y(1)*fu/(KR+y(1)*fu)+Vmaxbound*y(1)*fb/(KR+y(1)*fb);
ydot3=(k8*fu+k11*fb)*y(1)-kpt*y(3)+(k_1*fu+k_3*fb)*(To-y(1)-y(2)-y(3))*y(4)-Vmaxunbound*y(3)*fu/(KR+y(3)*fu)+Vmaxbound*y(3)*fb/(KR+y(3)*fb);
ydot4=kb*(To-y(1)-y(2))*(Bo-y(4))-k_b*y(4);
ydot5=ky*(To-y(1)-y(2))*(Yo-y(5))-k_y*y(5)*Zo;


dydt=[ydot1; ydot2; ydot3; ydot4; ydot5];