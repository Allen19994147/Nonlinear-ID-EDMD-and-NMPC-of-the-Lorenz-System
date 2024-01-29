%%% Nonlinear System Identification and Nonlienar Model Predictive Control
%%% Allen Lee
clc
clear all
%%% In this task, we use EDMDc to identify a Lorenz system and use NMPC to
%%% predict and control the behavior of the system using full state
%%% measurement.

global N n m x xr ur
global alpha  Rho beta Am Bm

alpha = 10; Rho = 28; beta = 8/3;

%%% System Parameters
n = 3; % # state
m = 2; % # input

x0 = [5 5 30]';      % IC
timestep = 0.001; RunTime = 10;
Sample_num = round(RunTime/timestep);
timeline = (0:Sample_num-1)*timestep;

X_dot_noisy = false;
X_noisy = false;
U_noisy = false;

ID_type = "state2state"; % approximate discrete system
% ID_type = "state2measurement"; % approximate continuous system

%%% Data Generation
X = zeros(length(x0),Sample_num); % record x
X_clean = zeros(size(X));
X_dot = zeros(length(x0),Sample_num); % record x
X_dot_clean = zeros(size(X_dot));
U = zeros(m,Sample_num); % record input
U_clean = zeros(size(U));
X_clean(:,1) = x0;
xm = x0;u = [0;0];
for i = 1:Sample_num
    u = -5+10.*rand(2,1);
    X_dot_clean(:,i) = Lorenz(X_clean(:,i),u);

    xm = X_clean(:,i) + timestep.*X_dot_clean(:,i);
    U_clean(:,i) = u;
    if(i==Sample_num)
        break
    end
    X_clean(:,i+1) = xm';
end

%%% Corrupted Measurement Data
if(X_dot_noisy)
    X_dot = X_dot_clean + wgn(size(X_dot_clean,1),size(X_dot_clean,2),-10); %0 dBW
else
    X_dot = X_dot_clean;
end
if(X_noisy)
    X = X_clean + wgn(size(X_dot_clean,1),size(X_dot_clean,2),-10);
else
    X = X_clean;
end
if(U_noisy)
    U = U_clean + wgn(size(U_clean,1),size(U_clean,2),-10);
else
    U = U_clean;
end
%% System Identification EDMDc

% Candidates Functions
Theta = Candidate_Library(X')';
Phi = Candidate_Library(U')';
switch(ID_type)
    case  "state2state"
        Theta_prime = Theta(:,2:end);
        Theta = Theta(:,1:end-1);
        ne = size(Theta,1); % to separate A,B later on
        % Theta_prime = A*Theta, A is a linear matrix
        Aext = Theta_prime*pinv([Theta;Phi(:,1:end-1)]);

    case "state2measurement"
        Theta_prime = Candidate_Library(X_dot')';
        ne = size(Theta,1); % to separate A,B later on
        % Theta_prime = A*Theta, A is a linear matrix
        Aext = Theta_prime*pinv([Theta;Phi]);
end


Am = Aext(:,1:ne);
Bm = Aext(:,ne+1:end);

%%% Data recovery
X_rcvA = zeros(size(X));
X_rcvA(:,1) = x0; % given IC
switch(ID_type)
    case "state2state"
        for i=1:length(Theta)
            dummy = Am*Candidate_Library(X_rcvA(:,i)')'...
                +Bm*Candidate_Library(U(:,i)')';
            X_rcvA(:,i+1) = dummy(1:n); % corres to x,y,z
        end
    case "state2measurement"
        for i=1:length(Theta)
            dummy = Am*Candidate_Library(X_rcvA(:,i)')'...
                +Bm*Candidate_Library(U(:,i)')';
            X_rcvA(:,i+1) = X_rcvA(:,i) + timestep.*dummy(1:n); % corres to x,y,z
            
        end
        X_rcvA = X_rcvA(:,1:end-1);
end
% X_rcvA = X_rcvA(:,1:end-1);
rmse(X',X_rcvA')


%% Data Comparison
for i=1:length(x0)
    figure
    hold on
    plot((0:Sample_num-1)*timestep,X(i,:));
    plot((0:Sample_num-1)*timestep,X_rcvA(i,:));
    legend('true','recovered')
    hold off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3D Plot
%%% patch
% plot3(X(:,1),X(:,2),X(:,3));
% plot3(X_rcvA(:,1),X_rcvA(:,2),X_rcvA(:,3));
% scatter3(X(:,1),X(:,2),X(:,3),s,c);
% scatter3(X_rcvA(:,1),X_rcvA(:,2),X_rcvA(:,3),s,c);
figure
hold on
colormap(hsv)
patch([X(1,:) nan],[X(2,:) nan],[X(3,:) nan],[timeline nan],'FaceColor','none','EdgeColor','interp')
patch([X_rcvA(1,:) nan],[X_rcvA(2,:) nan],[X_rcvA(3,:) nan],[timeline nan],'FaceColor','none','EdgeColor','interp')
colorbar
axis equal
xlabel('x(t)')
ylabel('y(t)')
zlabel('z(t)')
view(45,10)
hold off
%% plot3
figure
hold on
plot3(X(1,:),X(2,:),X(3,:),'r');
plot3(X_rcvA(1,:),X_rcvA(2,:),X_rcvA(3,:),'b');
axis equal
xlabel('x(t)')
ylabel('y(t)')
zlabel('z(t)')
view(45,10)
legend("true","recovered")
hold off
%% Animation of Lorentz System
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
% figure
% hold on
% colormap(hsv)
% colorbar
% axis equal
% xlabel('x(t)')
% ylabel('y(t)')
% zlabel('z(t)')
% view(45,10)
% xlim([-25,25])
% ylim([-25,25])
% zlim([0,50])
% for i=1:size(timeline,2)
%     scatter3(X(1,i),X(2,i),X(3,i),'r');
%     scatter3(X_rcvA(1,i),X_rcvA(2,i),X_rcvA(3,i),'b');
%     pause(0.00001)
% end
% hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Nonlinear Model Predictive Control
clc

N = 20; % Control/Prediction Horizen
xr = zeros(n,N+1); % State Reference
ur = zeros(m,N); % Manual

Run_Sample = 1000;
x0 = [-3.77 5.69 7.02]';
% x0 = [0 10 0]';
u0 = zeros(m,1);
Xall = zeros(n,Run_Sample+1); Uall = zeros(m,Run_Sample);
Xr = zeros(n,size(Xall,2)+N+1); Ur = zeros(m,size(Uall,2)+N);
x = x0; u = u0;

%%% Reference Tracking %%%
dummytimeline = (0:size(Xr,2)-1)*timestep;
Xr(1,:) = 10.*sin(1.*dummytimeline);
Xr(2,:) = 10.*cos(1.*dummytimeline); % References
for i=1:size(Ur,2)
    Ur(:,i) =  EquilibriumPoint(Xr(:,i));
end

%%
Xall = zeros(n,Run_Sample+1); Uall = zeros(m,Run_Sample);
Xall(:,1) = x0; Uall(:,1) = u0;
ur = Ur(:,1:1+N);
xr = Xr(:,1:1+N+1);

Xk = zeros(n*(N+1),1); % all future x from step k
Uk = zeros(m*N,1); % all future u from step k
Xk(1:n,1) = x0;
Uk(1:m,1) = u0;

zek = [Xk;Uk];
zek(1:n,1) = x0-xr(:,1); % Difference from the references
zek((N+1)*n+1:(N+1)*n+m) = u0-ur(:,1);

% Build QX,RU,H
xmin = [-10.1;-10.1;-100];
xmax = [10.1;10.1;100];
umin = [-100;-100].*1e2;
umax = [100;100].*1e2;

Fx = [eye(n);-eye(n)];
Fu = [eye(m);-eye(m)];

R = eye(m);
% Qx = 1.*eye(n);
Qx = 1e6.*diag([1 1 0]); % Don't care 3rd state
QX = zeros(n*(N+1),n*(N+1));
RU = zeros(m*N,m*N);
FX = zeros(2*n*(N+1),n*(N+1));
FU = zeros(2*m*N,m*N);
Gxe = zeros(2*n*(N+1),1);
Gue = zeros(2*m*N,1);
for i = 1:N
    QX((i-1)*n+1:i*n,(i-1)*n+1:i*n) = Qx;
    RU((i-1)*m+1:i*m,(i-1)*m+1:i*m) = R;

    FX((i-1)*2*n+1:i*2*n,(i-1)*n+1:i*n) = Fx;
    FU((i-1)*2*m+1:i*2*m,(i-1)*m+1:i*m) = Fu;
end

QX(N*n+1:end,N*n+1:end) = Qx;
FX(N*2*n+1:end,N*n+1:end) = Fx;
H = blkdiag(QX,RU.*0);

% Simulation with NMPC
ObjFunc = @(ze) (ze'*H*ze);
F = blkdiag(FX,FU);
Feq = []; Geq = []; lb = []; ub = []; % Linear constraints

nonlcon = @NonlinearConstraints;

options = optimoptions('fmincon','Display','none','Algorithm','sqp');
for i = 1:Run_Sample
    x = Xall(:,i);
    ur = Ur(:,i:i+N);
    xr = Xr(:,i:i+N+1);

    Gxe = zeros(2*n*(N+1),1);
    Gue = zeros(2*m*N,1);
    for j = 1:N
        Gxe((j-1)*2*n+1:j*2*n,1) = [xmax;-xmin] - Fx*xr(:,j);
        Gue((j-1)*2*m+1:j*2*m,1) =  [umax;-umin] - Fu*ur(:,j);
    end
    Gxe(N*2*n+1:end,1) = [xmax;-xmin] - Fx*xr(:,end);
    G = [Gxe;Gue];

    zek = fmincon(ObjFunc,zek,F,G,Feq,Geq,lb,ub,nonlcon,options);
    Uall(:,i) = zek(n*(N+1)+1:n*(N+1)+m,1) + ur(:,1);
    x_dot = Lorenz(Xall(:,i),Uall(:,i));
    if(i<=Run_Sample)
        Xall(:,i+1) =  Xall(:,i) +  x_dot.*timestep;
    end

end

%%
rmse(Xall(1:2,:)',Xr(1:2,1:size(Xall,2))')

%% Save and Figures
save("Lorentz_NMPC_Tracking.mat","Xall","Xr","Uall","Ur")

figure
hold on
plot(Xall(1,1:i))
plot(Xr(1,1:i))
hold off

figure
hold on
plot(Xall(2,1:i))
plot(Xr(2,1:i))
hold off

figure
hold on
plot([zeros(1,N) Xall(1,1:i)])
plot(Xr(1,1:i+N))
hold off

figure
hold on
plot([zeros(1,N) Xall(2,1:i)])
plot(Xr(2,1:i+N))
hold off
%% Visualization and Errors
timeframe = 0:timestep:(size(Xall,2)-1)*timestep;
figure
hold on
plot(timeframe,Xall(1,:),'r.')
plot(timeframe,Xr(1,1:size(Xall,2)),'b')
hold off
legend("Model","True","Noise")


figure
hold on
plot(timeframe,Xall(2,:),'r.')
plot(timeframe,Xr(2,1:size(Xall,2)),'b')
% plot(Xg(2,:),'g')
hold off
legend("Model","True","Noise")

figure
hold on
plot(timeframe,Xall(3,:),'r.')
plot(timeframe,Xr(3,1:size(Xall,2)),'b')
% plot(Xg(3,:),'g')
hold off
legend("Model","True","Noise")

figure
hold on
plot(Uall(1,:),'r.')
plot(Uall(2,:),'b.')
hold off
legend("u1","u2")



function x_dot = Lorenz(x,u)
global alpha Rho beta
x_dot = [alpha*(x(2)-x(1)) + 1/21.6*u(1)^3;
    Rho*x(1)-x(1)*x(3)-x(2) + u(2);
    x(1)*x(2)-beta*x(3)];
end


function [c,ceq] =NonlinearConstraints(z)
global N n m x xr ur Am Bm
c=[];
shfu = (N+1)*n;
for i=1:N
    x_dummy = z((i-1)*n+1:i*n,1)+xr(:,i);
    u_dummy = z(shfu+(i-1)*m+1:shfu+i*m)+ur(:,i);
    x_prime = Am*Candidate_Library(x_dummy')'+...
        Bm*Candidate_Library(u_dummy')';
    c1((i-1)*n+1:i*n,1)=z(i*n+1:(i+1)*n,1)+xr(:,i+1) - x_prime(1:n,1);

end
ceq = [z(1:n)-x+xr(:,i);c1];
end

function ur = EquilibriumPoint(xr)
global alpha Rho
ur = zeros(2,1);
ur(1) = nthroot(-21.6*alpha*(xr(2)-xr(1)),3);
ur(2) = -1*(Rho*xr(1)-xr(1)*xr(3)-xr(2));
end