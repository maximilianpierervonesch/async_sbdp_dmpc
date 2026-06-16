function J = integral_cost_function(z,N)

Q = eye(1);
R = 0.1;

offset_x= (N+1);
offset_u= N;

x1=z(1:N+1);
x2=z(offset_x+1:offset_x+N+1);
x3=z(2*offset_x+1:2*offset_x+N+1);

nX=3*(N+1);

u1=z(nX+1:nX+N);
u2=z(nX+offset_u+1:nX+offset_u+N);
u3=z(nX+2*offset_u+1:nX+2*offset_u+N);

J=zeros(N,1);

for k=1:N

    L1=x1(k)'*Q*x1(k) + R*u1(k)^2;
    L2=x2(k)'*Q*x2(k) + R*u2(k)^2;
    L3=x3(k)'*Q*x3(k) + R*u3(k)^2;

    J(k)=L1+L2+L3;

end

end