function [x1,x2,x3,u1,u2,u3,l1,l2,l3,V] = solve_central_ocp(cpl,x0,N,dt)

%% Centralized optimal control problem (M=3 agents)
%
% Decision variables:
%   z = [x(0:N); u(0:N-1)]
%
% dimensions
nx = 1;
nu = 1;
M  = 3;

%% Dimensions of optimization vector

nX = M * nx * (N+1);   % stacked states
nU = M * nu * N;       % stacked controls
nz = nX + nU;

z0 = zeros(nz,1);

%% Initial conditions

x10 = x0(1);
x20 = x0(2);
x30 = x0(3);

%% Bounds

lb = -inf(nz,1);
ub =  inf(nz,1);

% Control bounds
lb(nX+1:end) = -2;
ub(nX+1:end) =  2;

%% Fix initial conditions (state at time 0)

offset_x = N+1;
offset_u = N;

lb(1)              = x10;
ub(1)              = x10;

lb(offset_x+1)     = x20;
ub(offset_x+1)     = x20;

lb(2*offset_x+1)   = x30;
ub(2*offset_x+1)   = x30;

%% Solve NLP

options = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','iter', ...
    'MaxFunctionEvaluations',1e6);

[z,V,~,~,lambda] = fmincon( ...
    @(z) central_objective(z,N,dt), ...
    z0, ...
    [],[],[],[], ...
    lb,ub, ...
    @(z) central_constraints(z,cpl,N,dt), ...
    options);

%% =========================
% Extract primal solution
%% =========================

x1 = z(1:N+1);
x2 = z(offset_x+1 : offset_x+N+1);
x3 = z(2*offset_x+1 : 2*offset_x+N+1);

u1 = z(nX+1 : nX+N);
u2 = z(nX+offset_u +1 : nX+offset_u+N);
u3 = z(nX+2*offset_u +1 : nX+2*offset_u+N);

%% =========================
% Extract dual variables
%% =========================
%
% NOTE: assumes equality constraints are stacked as:
%   [dyn1(1); dyn2(1); dyn3(1); dyn1(2); ...]
%

l1 = [-lambda.lower(1);               lambda.eqnonlin(1:3:end)];
l2 = [-lambda.lower(offset_x+1);      lambda.eqnonlin(2:3:end)];
l3 = [-lambda.lower(2*offset_x+1);    lambda.eqnonlin(3:3:end)];