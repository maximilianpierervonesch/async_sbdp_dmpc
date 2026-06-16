function [c,ceq] = local_constraints(...
 z,...
 xNeigh,...
 cpl,...
 N,dt)

xi = z(1:N+1);
ui = z(N+2:end);

c = [];
ceq = zeros(N,1);

num_neighbors = size(xNeigh,1);

for k=1:N

    coupling = 0;

    for j=1:num_neighbors
        
        xj = xNeigh(j,k);
        coupling = coupling + cpl*sin(xj);

    end
    
    f = ui(k) + coupling;
    ceq(k) = xi(k+1) - (xi(k) + dt*f);

end

end