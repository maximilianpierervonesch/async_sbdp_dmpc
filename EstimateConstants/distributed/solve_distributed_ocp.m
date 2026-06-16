function [x1q,x2q,x3q,u1q,u2q,u3q,l1q,l2q,l3q,V,C_SBDP,L_SBDP] = ...
    solve_distributed_ocp(cpl,xk,x_warm,u_warm,l_warm,P_central,U_central,N,dt)

%% Distributed SBDP solver (Jacobi-type fixed-point iteration)
%
% Solves the distributed OCP using iterative local updates and warm starts.
% Compares iterates against centralized reference solution.
% Assumes that all three agents are coupled with eachother

%% Initialization (warm start)

u1q = u_warm(1,:);
u2q = u_warm(2,:);
u3q = u_warm(3,:);

x1q = x_warm(1,:);
x2q = x_warm(2,:);
x3q = x_warm(3,:);

l1q = l_warm(1,:);
l2q = l_warm(2,:);
l3q = l_warm(3,:);

xk = xk(:);

%% SBDP iteration parameters

maxSBDP = 5;

max_error = zeros(maxSBDP+1,1);
x_diff    = zeros(maxSBDP+1,1);
u_diff    = zeros(maxSBDP+1,1);
l_diff    = zeros(maxSBDP+1,1);

%% Central reference trajectories
% Structure assumption:
% P_central = [x1 x2 x3 l1 l2 l3]

p1_central = [P_central(:,1),P_central(:,4)]';
p2_central = [P_central(:,2),P_central(:,5)]';
p3_central = [P_central(:,3),P_central(:,6)]';

%% Initial error metrics (iteration 0)

error_p1 = max(vecnorm([x1q;l1q] - p1_central));
error_p2 = max(vecnorm([x2q;l2q] - p2_central));
error_p3 = max(vecnorm([x3q;l3q] - p3_central));

max_error(1) = max([error_p1,error_p2,error_p3]);

x_diff(1) = max(vecnorm([x1q;x2q;x3q] - P_central(:,1:3)'));
u_diff(1) = max(vecnorm([u1q;u2q;u3q] - U_central'));
l_diff(1) = max(vecnorm([l1q;l2q;l3q] - P_central(:,4:6)'));

%% SBDP fixed-point iterations

for q = 1:maxSBDP

    fprintf('SBDP Iteration %d\n', q);

    %% Neighbor information (Jacobi coupling)
    x1_neigh = [x2q; x3q];
    x2_neigh = [x1q; x3q];
    x3_neigh = [x1q; x2q];

    l1_neigh = [l2q; l3q];
    l2_neigh = [l1q; l3q];
    l3_neigh = [l1q; l2q];

    %% Local OCP solves
    [x1,u1,l1,V1] = solve_local_ocp(1,u1q,x1q,x1_neigh,l1_neigh,xk(1),cpl,N,dt);
    [x2,u2,l2,V2] = solve_local_ocp(2,u2q,x2q,x2_neigh,l2_neigh,xk(2),cpl,N,dt);
    [x3,u3,l3,V3] = solve_local_ocp(3,u3q,x3q,x3_neigh,l3_neigh,xk(3),cpl,N,dt);

    %% Error metrics w.r.t. centralized reference

    error_p1 = max(vecnorm([x1;l1] - p1_central));
    error_p2 = max(vecnorm([x2;l2] - p2_central));
    error_p3 = max(vecnorm([x3;l3] - p3_central));

    max_error(q+1) = max([error_p1,error_p2,error_p3]);

    x_diff(q+1) = max(vecnorm([x1;x2;x3] - P_central(:,1:3)','inf'));
    u_diff(q+1) = max(vecnorm([u1;u2;u3] - U_central','inf'));
    l_diff(q+1) = max(vecnorm([l1;l2;l3] - P_central(:,4:6)','inf'));

    %% Fixed-point update (currently no damping)

    alpha = 1;

    u1q = u1q + alpha*(u1 - u1q);
    u2q = u2q + alpha*(u2 - u2q);
    u3q = u3q + alpha*(u3 - u3q);

    x1q = x1q + alpha*(x1 - x1q);
    x2q = x2q + alpha*(x2 - x2q);
    x3q = x3q + alpha*(x3 - x3q);

    l1q = l1q + alpha*(l1 - l1q);
    l2q = l2q + alpha*(l2 - l2q);
    l3q = l3q + alpha*(l3 - l3q);

    %% Aggregate objective value
    V = V1 + V2 + V3;

end

%% ============================================================
% Contraction estimate from observed error decay
%% ============================================================

n_est = find(max_error(:) < 1e-4, 1);

if isempty(n_est)
    n_est = maxSBDP;
end

if n_est > 1
    C_SBDP = max(max_error(2:n_est) ./ max_error(1:n_est-1));
else
    C_SBDP = 0;
end

%% Sensitivity (global Lipschitz-style estimate)

L_SBDP = max( ...
    (x_diff + u_diff + l_diff) ./ max_error);