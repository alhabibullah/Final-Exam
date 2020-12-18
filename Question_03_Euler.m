%% Eular method for Nx=10
clc;clear all
% Initial conditions
x0=0;
Lx=15;
Nx=10;
x=linspace(x0,Lx,Nx);
dx=x(2)-x(1);
t0=0;
t_final=150;
Nt=100;
t=linspace(t0,t_final,Nt);
dt=t(2)-t(1);
alpha=1;

%Boundary condition
T=zeros(Nx,Nt);
T_steady=x.^2.*exp(-x);
S_x=-(x.^2-4.*x+2).*exp(-x);
T(:,1)=0;
T(1,:)=0;
T(:,end)=T_steady*Lx;


figure(1)
plot(x,T_steady,'--gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
title('plot of Tsteady for Nx=10')
xlabel('x');
ylabel('Tsteady');

for i =1:length(t)-1
    for j=2:length(x)-1
        Txx=(T(j-1,i)-2*T(j,i)+T(j+1,i))/dx^2;
        T(j,i+1)=T(j,i)+dt*(alpha*Txx+S_x(j));
    end
end
figure(2)
for m=2:2:length(x)
    txt=['at x= ',num2str(x(m))];
    plot(t,T(m,:),'LineWidth',1.5,'DisplayName',txt);
    hold on
end
xlabel('Time');
ylabel('Temparature')
legend show 
title('Explicity Eular scheme for Nx=10| Number of timestep 100')

%% Eular Method for Nx=20
clc;clear all
% Initial conditions
x0=0;
Lx=15;
Nx=20;
x=linspace(x0,Lx,Nx);
dx=x(2)-x(1);
t0=0;
t_final=150;
Nt=400;
t=linspace(t0,t_final,Nt);
dt=t(2)-t(1);
alpha=1;

%Boundary condition
T=zeros(Nx,Nt);
T_steady=x.^2.*exp(-x);
S_x=-(x.^2-4.*x+2).*exp(-x);
T(:,1)=0;
T(1,:)=0;
T(:,end)=T_steady*Lx;


figure(3)
plot(x,T_steady,'--gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
title('plot of Tsteady for Nx=20')
xlabel('x');
ylabel('Tsteady');

for i =1:length(t)-1
    for j=2:length(x)-1
        Txx=(T(j-1,i)-2*T(j,i)+T(j+1,i))/dx^2;
        T(j,i+1)=T(j,i)+dt*(alpha*Txx+S_x(j));
    end
end
figure(4)
for m=2:2:length(x)
    txt=['at x= ',num2str(x(m))];
    plot(t,T(m,:),'LineWidth',1.5,'DisplayName',txt);
    hold on
end
xlabel('Time');
ylabel('Temparature')
legend show 
title('Explicity Eular scheme for Nx=20|Timestep=400')

