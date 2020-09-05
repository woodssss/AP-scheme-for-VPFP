function PN_1D_evo(epsi,Nx,L,lambda,tau,T,IC)
%%%%%%%%%%%%% Proximal Quasi newton method for scaled Fokker planck %%%%%%
% This solver aims to solve
% \patial_t f = 1/\epsi FP(f).
% Stiffness comes from 1/\epsi
% tau is the outer time step
% Nx is the number of uniform grid point
% T is final time
% IC stand for initial condition, pick from 1-4
% Lambda is the stepsize of each inner iteration.
% Author: Wuzhe Xu
% Date: 09/05/2020
% Code of paper VARIATIONAL ASYMPTOTIC PRESERVING SCHEME FORVLASOV-POISSON-FOKKER-PLANCK SYSTEM
%%%%%%%%%%%%%
warning off
close all;
dx=2*L/Nx;
x=(-L+dx/2:dx:L)';
switch IC
    case 1
        f0 = exp(-(x-1).^2);
        disp('One bump to the left')
    case 2
        f0=exp(-(x+2).^2);
        disp('One bump to the right')
    case 3
        f0 = 2*exp(-5*(x-1.5).^2/1.2)+0.5*exp(-2*(x+1.5).^2/1.5);
        disp('Two bumps')
    case 4
        f0=2*exp(-4*(x-2.5).^2/1.2)+2*exp(-4*(x+2.5).^2/1.5);
        disp('Two bumps')
    otherwise
        error('Case not available, choose from {1,2,3,4}')
end
Mas=sum(f0)*dx;
M=exp(-(x).^2/2)*Mas/(2*pi)^0.5;

%% prepare A
A_rho=eye(Nx);
A_m=zeros(Nx,Nx);
h=1/2/dx;
for i = 2:Nx-1
    A_m(i,i+1)=h;
    A_m(i,i-1)=-h;
end
A_m(end,end-1)=-h;
A_m(1,2)=h;
A_m(1,1)=h;
A_m(end,end)=-h;
A=[A_rho A_m];

%% Initiate u0
f=f0;
num=1;
maxiter=200;
b=f;
t=0;
tol=1e-7;
while t<T
    u=[f;zeros(Nx,1)];
    tic;
    for k = 1:maxiter
        H=hess_diag(u,epsi,tau,Nx);
        HI=inverse(H);
        u_next=proxG_H(H,u-lambda*HI*grad(u,epsi,M,Nx,tau),A,b);
        [f,m]=decomp(u,Nx);
        if min(f)<0
            disp('f lose positivity, check parameter')
            return
        end
        error_E = abs(energyy(u_next,tau,dx,M,epsi,Nx)-energyy(u,tau,dx,M,epsi,Nx))/abs(energyy(u,tau,dx,M,epsi,Nx));
        error_u = abs(norm(u_next)-norm(u))/abs(norm(u));
        
        if (error_E <=tol && error_u<=tol)
            u=u_next;
            break
        end
        u=u_next;
        
    end
    toc;
    t=t+tau
    [f,m]=decomp(u,Nx);
    num=num+1;
    b=f;
end
plot(x,f,'r-*',x,M,'g',x,f0,'b')
legend('f(T,x)', 'Exact Equilibrium', 'f(0,x)')
drawnow
end

function [rho,m] = decomp(u,Nx)
rho=u(1:Nx);
m = u(Nx+1:end);
end

%% New energy with right entropy
function e=energyy(u,tau,dx,M,epsi,Nx)
[f,m]=decomp(u,Nx);
e1=epsi*m.^2./f+2*tau*f.*log(f./M);
e=sum(e1)*dx;
end

%% Proximal of F
function p=proxG_H(H,u,A, b)
HI=inverse(H);
lambda=(A*HI*A')^(-1)*(b-A*u);
p=u+HI*A'*lambda;
end


function HI=inverse(H)
[n,m]=size(H);
HI=0*H;
for i = 1:n
    if H(i,i)>0
        HI(i,i)=1/H(i,i);
    end
end
end

function g=grad(u,epsi,M,Nx,tau)
[f,m]=decomp(u,Nx);
gE1=-epsi*m.^2./f.^2;
gE2=2*tau*(log(f./M)+ones(Nx,1));
g1=gE1+gE2;
g2=2*epsi*m./f;
g1(f==0)=0; 
g2(f==0)=0; 
g=[g1;g2];
end


function H=hess_diag(u,epsi,tau,Nx)
[f,m]=decomp(u,Nx);
hh1=(2*epsi*m.^2./f.^3+2*tau./f);
hh4=(2*epsi./f);
hh1(hh1==0)=0;
hh4(hh4==0)=0;
H1=diag(hh1);
H4=diag(hh4);
H2=0*diag(-2*epsi*m./f.^2);
H=[H1 H2;H2 H4];
end


