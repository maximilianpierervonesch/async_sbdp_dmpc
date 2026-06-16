%% Estimate proof constants for the distributed OCP
%
% Optimal control problem:
%
%   min_u  \sum_{i \in V} \int_0^T
%           x_i(t)'Q_i x_i(t) + u_i(t)'R_i u_i(t) dt
%
% subject to
%
%   \dot{x}_i = u_i + cpl * \sum_{j \in N_i} sin(x_j),
%
%   u_i \in [-2,2],    for all i \in V.
%
% This script numerically estimates the constants appearing in the
% stability and convergence proofs for several representative initial
% conditions.

clear; clc; close all

plot_on = false;
export = false;

%% Problem setup

T  = 3;      % MPC prediction horizon
Tf = 3;      % Closed-loop simulation horizon

N  = 100;    % Number of discretization intervals
dt = T/N;    % Sampling time

nx = 1;      % States per agent
nu = 1;      % Inputs per agent

cpl = 1;     % Coupling strength

%% Initial conditions used for constant estimation

x01 = [0.75; 0.5; 0.25];
x02 = -[0.75; 0.5; 0.25];
x03 = [0.5; 0; 0.5];
x04 = [0; 0.5; 0];
x05 = [0.1; 0.4; 0.8];

x0_list = {x01,x02,x03,x04,x05};

nScen = numel(x0_list);

%% Storage for estimated proof constants

% Lyapunov function bounds:
%   a1 ||x||^2 <= V(x) <= a2 ||x||^2
% and local decrease rate a3.
a1_scen = zeros(nScen,1);
a2_scen = zeros(nScen,1);
a3_scen = zeros(nScen,1);

% Maximum variation of the Lyapunov function along the trajectory.
L_V_scen      = zeros(nScen,1);
L_V_sqrt_scen = zeros(nScen,1);

% Local Lipschitz constants of the dynamics.
L_f_x_scen = zeros(nScen,1);
L_f_u_scen = zeros(nScen,1);

% Sensitivity of the optimal solution with respect to the initial state.
L_u_scen = zeros(nScen,1);
L_p_scen = zeros(nScen,1);

% Contraction and Lipschitz constants of the SBDP iteration.
C_SBDP_scen = zeros(nScen,1);
L_SBDP_scen = zeros(nScen,1);

% Nominal MPC contraction (cf. Reble 2012)
alpha_scen = zeros(nScen,1);

%% Estimate constants for each initial condition

for i = 1:nScen

    % Current initial condition
    x0 = x0_list{i};

    %% Verify controllability assumption and for given MPC sampling time and prediction horizon
    [alpha_mpc,dt_mpc] = check_controllability(x0,cpl,Tf,dt);
    alpha_scen(i) = alpha_mpc;

    %% MPC parameters

    Nm   = round(dt_mpc/dt);    % Integration steps during one MPC update interval
    Nsim = round(Tf/dt_mpc);    % Number of MPC updates in the closed-loop simulation

    %% Storage for closed-loop trajectories and proof-constant estimates
    Xcl = zeros(3,Nsim+1);
    Ucl = zeros(3,Nsim);
    Vcl = zeros(1,Nsim);
    Vcl_sqrt = zeros(1,Nsim);

    xnorm_squared = zeros(1,Nsim);
    xnorm = zeros(1,Nsim);
    delta_x_norm = zeros(1,Nsim);
    U_diff = zeros(1,Nsim);
    P_diff = zeros(1,Nsim);

    Xcl_SBDP = zeros(3,Nsim+1);
    Ucl_SBDP = zeros(3,Nsim);
    C_SBDP_cl = zeros(1,Nsim);
    L_SBDP_cl = zeros(1,Nsim);
    P_diff_SBDP = zeros(1,Nsim);
    x_norm_diff = zeros(1,Nsim);
    f_cl_x = zeros(1,Nsim);

    % initialize
    Xcl(:,1) = x0;
    Xcl_SBDP(:,1) = x0;

    %% Initial guess for SBDP in first MPC iteration
    u1_warm = zeros(nu,N);
    u2_warm = zeros(nu,N);
    u3_warm = zeros(nu,N);

    x1_warm = zeros(nx,N+1);
    x2_warm = zeros(nx,N+1);
    x3_warm = zeros(nx,N+1);

    l1_warm = zeros(nx,N+1);
    l2_warm = zeros(nx,N+1);
    l3_warm = zeros(nx,N+1);

    x_warm = [x1_warm;x2_warm;x3_warm];
    u_warm = [u1_warm;u2_warm;u3_warm];
    l_warm = [l1_warm;l2_warm;l3_warm];

    %% MPC loop
    for k = 1:Nsim

        %% Current state
        xk = Xcl(:,k);

        %% Solve finite horizon centralized OCP
        if k>1
            U_prev = U;
            P_prev = P;
        end

        [X1p,X2p,X3p,U1p,U2p,U3p,L1p,L2p,L3p,V] = ...
            solve_central_ocp(cpl,xk,N,dt);

        % Save primal and dual solution
        U = [U1p,U2p,U3p];
        P = [X1p,X2p,X3p,L1p,L2p,L3p];

        % store difference
        if k>1
            % Maximum pointwise change of the optimal input trajectory
            U_diff(k-1) = max(vecnorm(U' - U_prev','inf'));

            % Maximum pointwise change of the primal-dual solution trajectory
            P_diff(k-1) = max(vecnorm(P' - P_prev','inf'));
        end

        %% Apply first part of optimal input
        uk1 = U1p(1:Nm);
        uk2 = U2p(1:Nm);
        uk3 = U3p(1:Nm);

        %% Propagate system
        x = xk;

        for j = 1:Nm

            u = [uk1(j); uk2(j); uk3(j)];

            %% System dynamics
            xdot = dynamics(x,u,cpl);
            x = x + dt*xdot;
        end

        %% Store closed-loop state
        Xcl(:,k+1) = x;

        %% Store applied input
        Ucl(:,k) = [uk1(1); uk2(1); uk3(1)];

        %% store metric for analysis
        Vcl(k) = V;                                             % optimal value function
        Vcl_sqrt(k) = sqrt(V);                                  % square root of optimal value
        xnorm_squared(k) = norm(Xcl(:,k),'inf')^2;                  % squared state norm
        xnorm(k) = norm(Xcl(:,k),'inf');                            % norm of state
        delta_x_norm(k) = norm(x - Xcl(:,k),'inf');                 % difference between two states
        % Numerator used for estimating the state Lipschitz constant
        f_cl_x(k) = norm(dynamics(x,Ucl(:,k),cpl) - dynamics(Xcl(:,k),Ucl(:,k),cpl),'inf'); % difference of dynamics in x

        %% Solve finite-horizon OCP using SBDP

        xk = Xcl_SBDP(:,k);

        % Compute the exact centralized solution at the current state.
        % This solution is used as a reference for estimating the
        % SBDP contraction and Lipschitz constants.
        [X1ref,X2ref,X3ref,U1ref,U2ref,U3ref,L1ref,L2ref,L3ref,V] = ...
            solve_central_ocp(cpl,xk,N,dt);

        U_ref = [U1ref,U2ref,U3ref];
        P_ref = [X1ref,X2ref,X3ref,L1ref,L2ref,L3ref];

        % Solve the distributed OCP using warm-start information from
        % the previous MPC iteration.
        [X1_SBDP,X2_SBDP,X3_SBDP,...
            U1_SBDP,U2_SBDP,U3_SBDP,...
            L1_SBDP,L2_SBDP,L3_SBDP,...
            V_SBDP,C_SBDP_k,L_SBDP_k] = ...
            solve_distributed_ocp( ...
            cpl,xk,x_warm,u_warm,l_warm,...
            P_ref,U_ref,N,dt);

        %% Update warm start for the next MPC iteration

        x_warm = [X1_SBDP;X2_SBDP;X3_SBDP];
        u_warm = [U1_SBDP;U2_SBDP;U3_SBDP];
        l_warm = [L1_SBDP;L2_SBDP;L3_SBDP];

        P_SBDP = [x_warm;l_warm];

        %% Extract control inputs applied during the next MPC interval

        uk1 = U1_SBDP(1:Nm);
        uk2 = U2_SBDP(1:Nm);
        uk3 = U3_SBDP(1:Nm);

        %% Propagate the system using the SBDP control input

        x = Xcl_SBDP(:,k);

        for j = 1:Nm
            
            u = [uk1(j); uk2(j); uk3(j)];

            % System dynamics
            xdot = dynamics(x,u,cpl);

            % Forward Euler integration
            x = x + dt*xdot;
        end

        %% Store closed-loop state and applied input

        Xcl_SBDP(:,k+1) = x;

        Ucl_SBDP(:,k) = [uk1(1); uk2(1); uk3(1)];

        %% Store quantities used for proof-constant estimation

        C_SBDP_cl(k) = C_SBDP_k;     % SBDP contraction estimate
        L_SBDP_cl(k) = L_SBDP_k;     % Local Lipschitz estimate of SBDP

        % Closed-loop deviation between centralized MPC and SBDP-MPC
        x_norm_diff(k) = norm(Xcl_SBDP(:,k+1) - Xcl(:,k+1),'inf');

        % Maximum pointwise deviation of the primal-dual solution
        P_diff_SBDP(k) = max(vecnorm(P_SBDP - P','inf'));
    end

    %% Closed-loop time vectors
    t_x = 0:dt_mpc:Tf;              % state sampling grid
    t_u = 0:dt_mpc:Tf-dt_mpc;      % input sampling grid

    %% Plot MPC closed loop (optional)

    if plot_on
        figure

        subplot(4,1,1)
        plot(t_x,Xcl(1,:), ...
            t_x,Xcl(2,:), ...
            t_x,Xcl(3,:), ...
            'LineWidth',2)
        ylabel('States')
        legend('x_1','x_2','x_3')

        subplot(4,1,2)
        stairs(t_u,Ucl(1,:),'LineWidth',2)
        ylabel('u_1')

        subplot(4,1,3)
        stairs(t_u,Ucl(2,:),'LineWidth',2)
        ylabel('u_2')

        subplot(4,1,4)
        stairs(t_u,Ucl(3,:),'LineWidth',2)
        ylabel('u_3')
        xlabel('time')
    end

    %% =========================
    %  Constant estimation
    %% =========================

    delta_V        = Vcl(2:end) - Vcl(1:end-1);
    delta_V_sqrt   = Vcl_sqrt(2:end) - Vcl_sqrt(1:end-1);

    % Lyapunov bounds: V(x) ~ ||x||^2
    a1_scen(i) = min(Vcl ./ xnorm_squared);
    a2_scen(i) = max(Vcl ./ xnorm_squared);

    % Lyapunov decrease rate (discrete)
    a3_scen(i) = min(-delta_V ./ xnorm_squared(1:end-1));

    %% Lyapunov smoothness / variation constants

    L_V_scen(i)      = max(abs(delta_V)      ./ xnorm(1:end-1));
    L_V_sqrt_scen(i) = max(abs(delta_V_sqrt) ./ xnorm(1:end-1));

    %% Dynamics Lipschitz estimates

    L_f_x_scen(i) = max(f_cl_x ./ delta_x_norm);
    L_f_u_scen(i) = 1;   % known from model structure

    %% Sensitivity of optimal solution w.r.t. state

    L_u_scen(i) = max(U_diff(1:end-1) ./ delta_x_norm(1:end-1));
    L_p_scen(i) = max(P_diff(1:end-1) ./ delta_x_norm(1:end-1));

    %% SBDP contraction/Lipschitz estimates

    C_SBDP_scen(i) = max(C_SBDP_cl);
    L_SBDP_scen(i) = max(L_SBDP_cl);

end

%% ============================================================
%  Worst-case constants over all scenarios
%% ============================================================

a1 = min(a1_scen);
a2 = max(a2_scen);
a3 = min(a3_scen);

L_V        = max(L_V_scen);
L_V_sqrt   = max(L_V_sqrt_scen);

L_f_x = max(L_f_x_scen);
L_f_u = max(L_f_u_scen);

L_u = max(L_u_scen);
L_p = max(L_p_scen);

C_SBDP   = max(C_SBDP_scen);
L_SBDP   = max(L_SBDP_scen);

alpha_mpc = min(alpha_scen);

%% Gronwall-type bounds

b1_prime = dt_mpc * L_f_u * L_SBDP * exp(dt_mpc * L_f_x);
b2_prime = dt_mpc * (L_f_x + L_f_u * L_u) * exp(dt_mpc * L_f_x);
b3_prime = b1_prime;

%% Auxiliary bounds entering contraction estimate

b1 = L_V_sqrt * b1_prime;
b2 = L_p * sqrt(a2) * b2_prime;
b3 = 1 + L_p * b3_prime;

%% Estimated iteration complexity

a_tilde = 1 - sqrt(max(0,(1 - a3/a2)));

% safeguard against invalid log arguments
arg = a_tilde / (a_tilde*b3 + b1*b2);

if arg <= 0 || C_SBDP <= 0
    q = Inf;
else
    q = 1 + log(arg) / log(C_SBDP);
end

%% Export

if export

    %% Export results for reproducibility and post-processing

    results = struct();

    % Scenario-wise constants
    results.a1_scen = a1_scen;
    results.a2_scen = a2_scen;
    results.a3_scen = a3_scen;

    results.L_V_scen      = L_V_scen;
    results.L_V_sqrt_scen = L_V_sqrt_scen;

    results.L_f_x_scen = L_f_x_scen;
    results.L_f_u_scen = L_f_u_scen;

    results.L_u_scen = L_u_scen;
    results.L_p_scen = L_p_scen;

    results.C_SBDP_scen = C_SBDP_scen;
    results.L_SBDP_scen = L_SBDP_scen;

    % Worst-case constants
    results.a1 = a1;
    results.a2 = a2;
    results.a3 = a3;

    results.L_V      = L_V;
    results.L_V_sqrt = L_V_sqrt;

    results.L_f_x = L_f_x;
    results.L_f_u = L_f_u;

    results.L_u = L_u;
    results.L_p = L_p;

    results.C_SBDP = C_SBDP;
    results.L_SBDP = L_SBDP;

    results.alpha_mpc = alpha_mpc;

    % Derived bounds
    results.b1 = b1;
    results.b2 = b2;
    results.b3 = b3;

    results.a_tilde = a_tilde;
    results.q = q;

    results.dt = dt;
    results.dt_mpc = dt_mpc;
    results.T = T;
    results.Tf = Tf;

    % Save to file
    save('dmcp_sbdp_constants.mat','results');
end