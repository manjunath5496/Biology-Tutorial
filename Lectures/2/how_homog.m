clear;
close;

n=25; % number of discretized spatial points
L=2;   % spatial length
dx=L/n; % incremental change in x coordinate
dt=0.01;   % incremental change in time
T=50;   % total integration time
tplot=0.01; % plot frequency

% constants

D_D=0.28;
D_E=0.6;
s1=20;
s1p=0.028;
s2=0.0063;
s3=0.04;
s4=0.8;
s4p=0.027;
Dtot=1500;
Etot=85;

% initial conditions


Dini=Dtot*rand;
dini=Dtot-Dini;
Eini=Etot*rand;
eini=Etot-Eini;

for x=1:n,
    Dold(x)=Dini;
    dold(x)=dini;
    Eold(x)=Eini;
    eold(x)=eini;
end;


% plot initial conditions

xvec=linspace(0,1,n);
subplot(2,2,1);
plot(xvec,dold,'b');
ylabel('d');
subplot(2,2,2);
plot(xvec,Dold,'g');
ylabel('D');
subplot(2,2,3);
plot(xvec,eold,'r');
ylabel('e');
subplot(2,2,4);
plot(xvec,Eold,'y');
ylabel('E');

drawnow;
        
input('Hit a key to start simulation');
close;
figure;


% start of iterative time loop

for timestep=1:(T/dt);
    
    t=timestep*dt;
    
    % equation for internal points (from 2 to n-1)
    for x=2:n-1,
        dDdt(x)=-s1*Dold(x)/(1+s1p*eold(x))+s2*eold(x)*dold(x)+...
            D_D*(Dold(x-1)+Dold(x+1)-2*Dold(x))/dx^2;
        dddt(x)=s1*Dold(x)/(1+s1p*eold(x))-s2*eold(x)*dold(x);        
        dEdt(x)=s4*eold(x)/(1+s4p*Dold(x))-s3*Dold(x)*Eold(x)+...
            D_E*(Eold(x-1)+Eold(x+1)-2*Eold(x))/dx^2;
        dedt(x)=-s4*eold(x)/(1+s4p*Dold(x))+s3*Dold(x)*Eold(x);
    end;
        dDdt(1)=-s1*Dold(1)/(1+s1p*eold(1))+s2*eold(1)*dold(1)+...
            D_D*(Dold(2)-Dold(1))/dx^2;
        dddt(1)=s1*Dold(1)/(1+s1p*eold(1))-s2*eold(1)*dold(1);        
        dEdt(1)=s4*eold(1)/(1+s4p*Dold(1))-s3*Dold(1)*Eold(1)+...
            D_E*(Eold(2)-Eold(1))/dx^2;
        dedt(1)=-s4*eold(1)/(1+s4p*Dold(1))+s3*Dold(1)*Eold(1);
        
        dDdt(n)=-s1*Dold(n)/(1+s1p*eold(n))+s2*eold(n)*dold(n)+...
            D_D*(Dold(n-1)-Dold(n))/dx^2;
        dddt(n)=s1*Dold(n)/(1+s1p*eold(n))-s2*eold(n)*dold(n);        
        dEdt(n)=s4*eold(n)/(1+s4p*Dold(n))-s3*Dold(n)*Eold(n)+...
            D_E*(Eold(n-1)-Eold(n))/dx^2;
        dedt(n)=-s4*eold(n)/(1+s4p*Dold(n))+s3*Dold(n)*Eold(n);

    % update the new values.
    for x=1:n,
        dnew(x,timestep)=dold(x)+dddt(x)*dt;
        Dnew(x,timestep)=Dold(x)+dDdt(x)*dt;
        enew(x,timestep)=eold(x)+dedt(x)*dt;
        Enew(x,timestep)=Eold(x)+dEdt(x)*dt;
        dold(x)=dold(x)+dddt(x)*dt;
        Dold(x)=Dold(x)+dDdt(x)*dt;
        eold(x)=eold(x)+dedt(x)*dt;
        Eold(x)=Eold(x)+dEdt(x)*dt;
    end;
    
    if rem(t,tplot)==0

xvec=linspace(0,1,n);
subplot(2,2,1);
plot(xvec,dold,'b');
ylabel('d');
title(dold(1));
subplot(2,2,2);
plot(xvec,Dold,'g');
ylabel('D');
title(Dold(1));
subplot(2,2,3);
plot(xvec,eold,'r');
ylabel('e');
title(eold(1));
subplot(2,2,4);
plot(xvec,Eold,'y');
ylabel('E');
title(Eold(1));
drawnow;
end;   
end;
