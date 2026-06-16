function [X,U,L,V] = solve_local_ocp( ...
    id, ...
    uiq, ...
    xiq, ...
    xNeigh, ...
    lambdaNeigh, ...
    xk, ...
    cpl, ...
    N, dt)

%% Local optimal control problem (agent-wise SBDP subproblem)
%
% Each agent solves a reduced OCP using:
% - its own warm start
% - neighbor trajectories (Jacobi coupling)
% - fixed initial condition xk

nx = 1;
nu = 1;

nX = nx * (N+1);
nU = nu * N;

nz = nX + nU;

%% Initial guess (warm start)

z0 = zeros(nz,1);

% states warm start
z0(1:nX) = xiq;

% controls warm start
z0(nX+1:end) = uiq;

%% Bounds

lb = -inf(nz,1);
ub =  inf(nz,1);

% control constraints
lb(nX+1:end) = -2;
ub(nX+1:end) =  2;

% enforce initial condition
lb(1) = xk;
ub(1) = xk;

%% Solve NLP

options = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','off', ...
    'MaxFunctionEvaluations',1e6);

[z,V,~,~,lambda] = fmincon( ...
    @(z) local_objective( ...
        z, xiq, xNeigh, lambdaNeigh, cpl, N, dt), ...
    z0, ...
    [],[],[],[], ...
    lb,ub, ...
    @(z) local_constraints( ...
        z, xNeigh, cpl, N, dt), ...
    options);

%% Extract primal variables

X = z(1:N+1)';          % state trajectory
U = z(nX+1:end)';       % control trajectory

%% Extract dual variables

% NOTE:
% lambda.eqnonlin assumed to correspond to dynamics constraints
% ordered sequentially in time

L = [-lambda.lower(1); lambda.eqnonlin(:)]';