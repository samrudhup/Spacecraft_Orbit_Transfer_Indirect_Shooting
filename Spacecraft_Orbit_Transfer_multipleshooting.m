clear all; close all; clc;
%OTMain
%Know conditions
comp_t  =tic;
T       = 0.1405;
Ve      = 1.8758344;
mu      = 1;
nx      = 5; %No of States
K       = 16; %No of intervals
%boundary conditions
r_0     = 1;
r_f     = 1.5;
theta_0 = 0;
vr_0    = 0;
vr_f    = 0;
vt_0    = sqrt(mu/r_0);
vt_f    = sqrt(mu/r_f);
m_0     = 1;
t_0     = 0;
t_f_guess = 3.5;

betaj=[];
%Initial Lambda guesses
lambda_r_0      = -2;
lambda_theta_0  = 0;
lambda_vr_0     = 2;
lambda_vt_0     = 2;
lambda_m_0      = 2;

P_0_guess = ones(nx*2, K-1);
lambda_0_guess = [lambda_r_0; lambda_theta_0; lambda_vr_0; lambda_vt_0; lambda_m_0; t_f_guess; P_0_guess(:)];
tau = linspace(-1,1,K+1);
options = optimset('Display','Iter','TolX',1e-8,'TolFun',1e-8);
z = fsolve(@OTError, lambda_0_guess, options, r_0, theta_0, vr_0, vt_0, m_0, r_f, vr_f, vt_f, T, mu, Ve, nx, K, tau, t_0);
%plot using integrated dynamics
[E,t,p] = OTError(z, r_0, theta_0, vr_0, vt_0, m_0, r_f, vr_f, vt_f, T, mu, Ve, nx, K, tau, t_0);

%Plots
figure (1)
plot(t,p(:,1),'r',t,p(:,2),'b',t,p(:,3),'m',t,p(:,4),'k', t, p(:,5),'g');
title('Spacecraft Orbit Transfer- State plots')
xlabel('time')
ylabel('States')
legend('r','theta','vr','vt','mass','Location','best');
grid on
hold on

%Recalculating beta for plot
for i=1:length(t)
    betaj(i)= atan2(p(i,8),p(i,9));
    betaj = unwrap(betaj);
   
end
figure (2)
plot(t,betaj);
title('Spacecraft Orbit Transfer- Input plot')
xlabel('time')
ylabel('Input')
legend('beta','Location','best');
grid on
hold on

figure (3)
polarplot(p(:,2),p(:,1));
hold on;

comp_tf=toc(comp_t);

%Orbit_transfer_myODE
function [P_dot]=OTODE(~,P,T,mu,Ve,t_0, t_f)
r       = P(1);
theta   = P(2);
vr      = P(3);
vt      = P(4);
m       = P(5);
lambda_r =P(6);
lambda_theta = P(7);
lambda_vr   = P(8);
lambda_vt   = P(9);
lambda_m    = P(10);


%solve for beta
beta    = atan2(lambda_vr,lambda_vt);

%first-order differential equations
r_dot   = vr;
theta_dot = vt/r;
vr_dot  = (T/m)*sin(beta) + ((vt^2)/r) - mu/r^2;
vt_dot  = (T/m)*cos(beta) - vr*vt/r;
m_dot   = -T/Ve;
lambda_r_dot =lambda_theta*vt/(r^2) + lambda_vr*(vt^2)/(r^2) - 2*lambda_vr*mu/(r^3) - lambda_vt*vr*vt/(r^2);
lambda_theta_dot = 0;
lambda_vr_dot = -lambda_r + lambda_vt*vt/r;
lambda_vt_dot = -lambda_theta/r - (2*vt*lambda_vr/r) + lambda_vt*vr/r;
lambda_m_dot = lambda_vt*T*cos(beta)/(m^2) + lambda_vr*T*sin(beta)/m^2;

P_dot   = [r_dot; theta_dot; vr_dot; vt_dot; m_dot; lambda_r_dot; lambda_theta_dot; lambda_vr_dot; lambda_vt_dot; lambda_m_dot];
P_dot   = (t_f-t_0)/2 * P_dot; %Scaling from tau to t

end

%OTError
function [E,t,P] = OTError(z_0, r_0, theta_0, vr_0, vt_0, m_0, r_f, vr_f, vt_f, T, mu, Ve, nx, K, tau, t_0)
t_f         = z_0(6); 
P_0_guess   = z_0(7:end);
P_0_guess   = reshape(P_0_guess, 2*nx, K-1);
options     = odeset('RelTol', 1e-8);
E=[];
t=[];
P=[];

for k=1:K
    if k==1
        %for the first interval since the inital boundary conditions are
        %known
        P_0 = [r_0; theta_0; vr_0; vt_0; m_0; z_0(1:5)];
    else
        P_0 = P_0_guess(:,k-1);
    end

    tspan =[tau(k), tau(k+1)];
    [tout,Pout] = ode113(@OTODE, tspan, P_0, options, T, mu, Ve, t_0, t_f);
    Pth = Pout(end,:).';
   
    if k<K
        E = [E;Pth-P_0_guess(:,k)];
    end
    t=[t;tout];
    P=[P;Pout];
end

r_tf        =Pth(1);
theta_tf    = Pth(2);
vr_tf       = Pth(3);
vt_tf       = Pth(4);
m_tf        = Pth(5);
lambda_r_tf = Pth(6);
lambda_theta_tf = Pth(7);
lambda_vr_tf = Pth(8);
lambda_vt_tf = Pth(9);
lambda_m_tf = Pth(10);

beta_tf     = atan2(lambda_vr_tf,lambda_vt_tf);
H_tf        = -1 + lambda_r_tf*vr_tf +lambda_theta_tf*vt_tf/r_tf +lambda_vr_tf*(sin(beta_tf)*(T/m_tf)+(vt_tf^2)/r_tf-mu/r_tf^2)+lambda_vt_tf*(cos(beta_tf)*(T/m_tf)-vr_tf*vt_tf/r_tf)-lambda_m_tf*(T/Ve); %beta uncertainty

E = [ E; r_tf-r_f; vr_tf-vr_f; vt_tf-vt_f; H_tf; lambda_theta_tf; lambda_m_tf-1];
end


