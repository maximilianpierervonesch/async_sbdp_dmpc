function [c,ceq] = central_constraints(z,cpl,N,dt)

c=[];

ceq      = zeros(3*N,1);
offset_x = (N+1);
offset_u = N;

x1 = z(1:N+1);
x2 = z(offset_x+1:offset_x+N+1);
x3 = z(2*offset_x+1:2*offset_x+N+1);

nX = 3*(N+1);

u1 = z(nX+1:nX+N);
u2 = z(nX+offset_u+1:nX+offset_u+N);
u3 = z(nX+2*offset_u+1:nX+2*offset_u+N);

for k=1:N

    %% Agent 1 continous time dynamics
    f1 = u1(k) + cpl * sin(x2(k)) + cpl * sin(x3(k));

    %% Agent 2 continous time dynamics
    f2 = u2(k) + cpl * sin(x1(k)) + cpl * sin(x3(k));

    % Agent 3 continous time dynamics
    f3 = u3(k) + cpl * sin(x1(k)) + cpl * sin(x2(k));


    ceq((k-1)*3 + 1:k*3) = ...
        [x1(k+1)-(x1(k)+dt*f1);
        x2(k+1)-(x2(k)+dt*f2);
        x3(k+1)-(x3(k)+dt*f3);
        ];
end

end