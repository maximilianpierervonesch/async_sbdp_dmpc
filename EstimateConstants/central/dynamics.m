function x_dot = dynamics(x,u,cpl)

     x_dot = [   
            u(1) + cpl * sin(x(2)) + cpl * sin(x(3));
            u(2) + cpl * sin(x(1)) + cpl * sin(x(3));
            u(3) + cpl * sin(x(1)) + cpl * sin(x(2));];
end 