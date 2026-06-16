function J = local_objective(...
 z,...
 xiq,...
 xNeigh,...
 lambdaNeigh,...
 cpl,...
 N,dt)

% extract
xi=z(1:N+1);
ui=z((N+1)+1:end);

% Cost matrices 
R = 0.1;
Q = 1;

% number of neighbors 
num_neighbors = size(xNeigh,1);

J=0;

for k=1:N
    %% nominal stage cost
    l = Q*xi(k)^2 + R*ui(k)^2;

    %% Example Coulpling Jacobian
    g_ij = 0;
    for j=1:num_neighbors
        lambdaj = lambdaNeigh(j,k+1);  % first multiplier is for initial condition

         dfjdxi = cpl * cos(xiq(k));
        
         g_ij = g_ij - dfjdxi*lambdaj*(xi(k)-xiq(k));
    end 

    J = J + dt * (l + g_ij);
end

end