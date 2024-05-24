function yprpr = hessianfunc(x)
    syms x1 x2;
    F=x1.^2+ 2*x2.^2-0.3*cos(3*pi*x1+4*pi*x2)+0.3;
    H=hessian(F);
    yprpr=double(subs(H,[x1,x2],[x(1),x(2)]));
end
