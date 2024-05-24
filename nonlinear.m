clear all
close all
clc

X1=-100:0.1:100;
X2=-100:0.1:100;
[x1,x2]=meshgrid(X1,X2);
F=x1.^2+ 2*x2.^2-0.3*cos(3*pi*x1+4*pi*x2)+0.3;
realFMin=min(min(F))
mesh(x1,x2,F)

figure
contourf(x1,x2,F)
hold on

% %% Newton-Raphson
 fprintf('Newton-Raphson Algorithm\n');
  x=rand(2,1);
% x=[1;1];
 epsilon=10^(-4);
 
 tic
 fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n',x(1),x(2),func(x))
 plot(x(1),x(2),'r.')
 x_next=x-inv(hessianfunc(x))*gradfunc(x);
 fprintf('k=2, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
 plot(x_next(1),x_next(2),'r*')
 k=3;
 %while(abs(func(x_next)-func(x))>epsilon)
 while(norm(gradfunc(x_next))>epsilon)
 x=x_next;
 x_next=x-inv(hessianfunc(x))*gradfunc(x);
 fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',k,x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
 plot(x_next(1),x_next(2),'r*')
 k=k+1;
 end
toc
title('Newton-Raphson Algorithm')
 set(gca,'fontsize',35)
 set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

% %% Hestenes-Stiefel Algorithm
 figure
 contourf(x1,x2,F)
 hold on
% 

 fprintf('Hestenes-Stiefel Algorithm\n');
 x=rand(2,1);
% x=[1;1];
epsilon=10^(-4);

tic
 fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n',x(1),x(2),func(x))
 plot(x(1),x(2),'r.')
g=gradfunc(x);
d=-g;

alpha=0:0.01:1;
funcalpha=zeros(length(alpha),1);
for i=1:length(alpha)
    funcalpha(i)=func(x+alpha(i)*d);
end    
[val,ind]=min(funcalpha);
alpha=alpha(ind);
x_next=x+alpha*d;
g_next=gradfunc(x_next);
beta=(g_next'*(g_next-g))/(d'*(g_next-g));
d_next=-g_next+beta*d;
fprintf('k=2, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
plot(x_next(1),x_next(2),'r*')
k=3;
 while(norm(gradfunc(x_next))>epsilon)
     x=x_next;
     g=g_next;
     d=d_next;
      alpha=0:0.01:1;
      funcalpha=zeros(length(alpha),1);
        for i=1:length(alpha)
             funcalpha(i)=func(x+alpha(i)*d);
        end    
              [val,ind]=min(funcalpha);
              alpha=alpha(ind);

  x_next=x+alpha*d;
  g_next=gradfunc(x_next);
  beta=(g_next'*(g_next-g))/(d'*(g_next-g));
  d_next=-g_next+beta*d;

  fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',k,x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
  plot(x_next(1),x_next(2),'r*')
  k=k+1;
 end

 toc
 title('Hestenes-Stiefel Algorithm')
 set(gca,'fontsize',35)
 set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

 % %% Polak-Ribi`ere Algorithm
 figure
 contourf(x1,x2,F)
 hold on
% 

 fprintf('Polak-Ribi`ere\n');
 x=rand(2,1);
% x=[1;1];
epsilon=10^(-4);

tic
 fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n',x(1),x(2),func(x))
 plot(x(1),x(2),'r.')
g=gradfunc(x);
d=-g;

alpha=0:0.01:1;
funcalpha=zeros(length(alpha),1);
for i=1:length(alpha)
    funcalpha(i)=func(x+alpha(i)*d);
end    
[val,ind]=min(funcalpha);
alpha=alpha(ind);
x_next=x+alpha*d;
g_next=gradfunc(x_next);
  beta=(g_next'*(g_next-g))/(g'*g);
d_next=-g_next+beta*d;
fprintf('k=2, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
plot(x_next(1),x_next(2),'r*')
k=3;
 while(norm(gradfunc(x_next))>epsilon)
     x=x_next;
     g=g_next;
     d=d_next;
      alpha=0:0.01:1;
      funcalpha=zeros(length(alpha),1);
        for i=1:length(alpha)
             funcalpha(i)=func(x+alpha(i)*d);
        end    
              [val,ind]=min(funcalpha);
              alpha=alpha(ind);

  x_next=x+alpha*d;
  g_next=gradfunc(x_next);
  beta=(g_next'*(g_next-g))/(g'*g);
  d_next=-g_next+beta*d;

  fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',k,x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
  plot(x_next(1),x_next(2),'r*')
  k=k+1;
 end

 toc
 title('Polak-Ribi`ere')
 set(gca,'fontsize',35)
 set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

  % %% Fletcher-Reeves Algorithm
 figure
 contourf(x1,x2,F)
 hold on
% 

 fprintf('Fletcher-Reeves\n');
 x=rand(2,1);
% x=[1;1];
epsilon=10^(-4);

tic
 fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n',x(1),x(2),func(x))
 plot(x(1),x(2),'r.')
g=gradfunc(x);
d=-g;

alpha=0:0.01:1;
funcalpha=zeros(length(alpha),1);
for i=1:length(alpha)
    funcalpha(i)=func(x+alpha(i)*d);
end    
[val,ind]=min(funcalpha);
alpha=alpha(ind);
x_next=x+alpha*d;
g_next=gradfunc(x_next);
beta=(g_next'*g_next)/(g'*g);
d_next=-g_next+beta*d;
fprintf('k=2, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
plot(x_next(1),x_next(2),'r*')
k=3;
 while(norm(gradfunc(x_next))>epsilon)
     x=x_next;
     g=g_next;
     d=d_next;
      alpha=0:0.01:1;
      funcalpha=zeros(length(alpha),1);
        for i=1:length(alpha)
             funcalpha(i)=func(x+alpha(i)*d);
        end    
              [val,ind]=min(funcalpha);
              alpha=alpha(ind);

  x_next=x+alpha*d;
  g_next=gradfunc(x_next);
  beta=(g_next'*g_next)/(g'*g);
  d_next=-g_next+beta*d;

  fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',k,x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
  plot(x_next(1),x_next(2),'r*')
  k=k+1;
 end

 toc
 title('Fletcher-Reeves')
 set(gca,'fontsize',35)
 set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);