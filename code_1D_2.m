clc;
clear;


a=20;
x0=-a; x1=a; h=1/20; xx=x0:h:x1; m=length(xx);

t0=0; dt=h; t1=0.9; tt=t0+dt:dt:t1 ; n=length(tt); 

x00=0; %initial point

f=@(x) 1;%drift term

p1=@(x,y) ((2*pi*dt).^(-1/2))*exp(-(x-y-f(y)*dt).^2/(2*dt));%General formula of probability density


P=zeros(n-1,m);
P1=zeros(m,1);
P2=zeros(m,1);
Q=zeros(m,m);




P1=p1(xx,x00)';
P(1,:)=P1;

for i=1:m
    Q(:,i)=p1(xx,xx(i))'; 
end

for k=2:n
P2=Q*P1*h; 
P1=P2;
P(k,:)=P2;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%exact solution 
% p=@(x,t) (pi*(1-exp(-2*t))).^(-1/2)*exp(-(x-exp(-t)*x0).^2/(1-exp(-2*t)));
p=@(x,t) (2*pi*t).^(-1/2)*exp(-(x-x00-t).^2/(2*t)); %drift term f==1;
for k=1:n
    for i=1:m
     T(k,i)=p(xx(i),tt(k));
    end
end
 l=n;
plot(xx,P(l,:),'.-')
hold on
plot(xx,T(l,:),'.-')

