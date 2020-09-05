function PN_VPFP_1D_mix_muscl(Nx,Nv,Lx,Lv,tau,T,opt1)
%%%%%%%%%%%%% Proximal Quasi newton method for mix regime VPFP %%%%%%
% This solver aims to solve
% \patial_t f + v \nabla_x f = 1/\epsi ( \nabla_x \phi \nabla_v f + FP(f)).
% -\Delta_x \phi = \rho -h 
% tau is the preset outer time step
% Nx Nv are the number of uniform grid point
% T is final time
% Author: Wuzhe Xu
% Date: 09/05/2020
% Code of paper VARIATIONAL ASYMPTOTIC PRESERVING SCHEME FORVLASOV-POISSON-FOKKER-PLANCK SYSTEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
epsi0=1e-3;
dv=2*Lv/Nv;
v=(-Lv+dv/2:dv:Lv)';
dx = 2*Lx/Nx;
x=(-Lx+dx/2:dx:Lx)';
tau=min(dx/15,tau);
dt = tau;
h_x = 1.6711/2.5321*exp(cos(pi*x));
rho_x0 = 1/6*(2+sin(pi*(x)));
epsix=x;
if nargin<7
    opt1=0;
    disp('Defualt setting: do not show result of each step')
end
for i=1:length(x)
    if x(i)<=0.3
        epsix(i)=epsi0+0.5*(tanh(5-10*x(i))+tanh(5+10*x(i)));
    else
        epsix(i)=epsi0;
    end
end
phi_x0 = getdp(rho_x0-h_x,Nx);
f0 = zeros(Nv,Nx);
for i = 1 :Nx
    for j = 1:Nv
        f0(j,i) = rho_x0(i)*exp(-(v(j)+phi_x0(i))^2/2);
    end
end
t=0;
f=f0;
M = zeros(Nv,Nx);
rho=zeros(1,Nx);
%% prepare A
A_rho=eye(Nv);
A_m=zeros(Nv,Nv);
h=1/2/dv;
for i = 2:Nv-1
    A_m(i,i+1)=h;
    A_m(i,i-1)=-h;
end
A_m(end,end-1)=-h;
A_m(1,2)=h;
A_m(1,1)=h;
A_m(end,end)=-h;
A=[A_rho A_m];
tfinal=T;
count=0;
rho_new=zeros(1,Nx);
while t < tfinal
    %% Advection by MUSCL
    f=MUSCL(v,f,Nx,dx,dt);
    %% get grad_theta_x and M
    for i = 1:Nx
        rho(i) = sum(f(:,i))*dv;
    end
    grad_theta_x = getdp(rho'-h_x,Nx);
    for i = 1:Nx
        for j = 1:Nv
            M(j,i) = exp(-(v(j)+grad_theta_x(i))^2/2);
        end
    end
    for i = 1:Nx
        M(:,i) = rho(i)*M(:,i);
    end
    
    %% Collision
    parfor i = 1:Nx
        f(:,i) = PN(f(:,i),M(:,i),rho(i),epsix(i),Nv,A,tau,dv);
    end
    if opt1 ==1
        for i = 1:Nx
            rho_new(i)=sum(f(:,i))*dv;
        end
        plot(x,rho_new)
        drawnow
    end
    t = t+dt;
    count=count+1;
    
    if rem(t,0.1)<1e-5
        for i = 1:Nx
            rho_new(i)=sum(f(:,i))*dv;
        end
        plot(x,rho_new)
        drawnow
        filename=['PN_vpfp_t', num2str_decimal(t)];
        save(filename)
        disp(['t=',num2str(t)])
    end
end

end


function f_out = MUSCL(v,f,Nx,dx,dt)
f=f';
f1=f;
for k = 3 : Nx-2
    f1(k,:) = f(k,:) -dt/dx*max(v',0).*(f(k,:)-f(k-1,:)) - dt/dx*min(v',0).*(f(k+1,:)-f(k,:))...
        -dt/dx*(Fluc(f(k,:), f(k+1,:), f(k-1,:), f(k+2,:), v', dx, dt)...
        -Fluc(f(k-1,:), f(k,:), f(k-2,:), f(k+1,:), v', dx, dt));
end
f1(1,:) = f(1,:) -dt/dx*max(v',0).*(f(1,:)-f(Nx,:)) - dt/dx*min(v',0).*(f(2,:)-f(1,:))...
    -dt/dx*(Fluc(f(1,:), f(2,:), f(Nx,:), f(3,:), v', dx, dt)...
    -Fluc(f(Nx,:), f(1,:), f(Nx-1,:), f(2,:), v', dx, dt));
f1(2,:) = f(2,:) -dt/dx*max(v',0).*(f(2,:)-f(1,:)) - dt/dx*min(v',0).*(f(3,:)-f(2,:))...
    -dt/dx*(Fluc(f(2,:), f(3,:), f(1,:), f(4,:), v', dx, dt)...
    -Fluc(f(1,:), f(2,:), f(Nx,:), f(3,:), v', dx, dt));
f1(Nx-1,:) = f(Nx-1,:) -dt/dx*max(v',0).*(f(Nx-1,:)-f(Nx-2,:)) - dt/dx*min(v',0).*(f(Nx,:)-f(Nx-1,:))...
    -dt/dx*(Fluc(f(Nx-1,:), f(Nx,:), f(Nx-2,:), f(1,:), v', dx, dt)...
    -Fluc(f(Nx-2,:), f(Nx-1,:), f(Nx-3,:), f(Nx,:), v', dx, dt));
f1(Nx,:) = f(Nx,:) -dt/dx*max(v',0).*(f(Nx,:)-f(Nx-1,:)) - dt/dx*min(v',0).*(f(1,:)-f(Nx,:))...
    -dt/dx*(Fluc(f(Nx,:), f(1,:), f(Nx-1,:), f(2,:), v', dx, dt)...
    -Fluc(f(Nx-1,:), f(Nx,:), f(Nx-2,:), f(1,:), v', dx, dt));
f_out = f1';
end



function e=energyy(u,M,epsi,tau,dv,Nv)
[f,m]=decomp(u,Nv);
e1=epsi*m.^2./f+2*tau*f.*log(f./M);
e=sum(e1)*dv;
end


function [rho,m] = decomp(u,Nv)
rho=u(1:Nv);
m = u(Nv+1:end);
end

function dphi=getdp(r,Nx)
rho_hat =  fft(r);
phi_hat = [1 1:Nx/2 -Nx/2+1:-1]'.^(-2).* rho_hat;
phi_hat(1) = 1;  %this value is assigned to 0 mode of phi_hat
dphi = real(ifft(phi_hat.*[0 1:Nx/2 -Nx/2+1:-1]'*1i))/pi;
end


function H=hess_diag(u,epsi,tau,Nv)
[f,m]=decomp(u,Nv);
f=max(f,1e-10);
H1=diag(2*epsi*m.^2./f.^3+2*tau./f);
%H1=min(H1,1e10);
H2=0*diag(-2*epsi*m./f.^2);
H4=diag(2*epsi./f);
H=[H1 H2;H2 H4];
end

function p=proxG(H,u,b,A)
lambda=(A*H^(-1)*A')^(-1)*(b-A*u);
p=u+H^(-1)*A'*lambda;
end


function fv=PN(f_in,M_in,mass,epsi,Nv,A,tau,dv)
f=f_in/mass;
M=M_in/mass;
u=[f;zeros(Nv,1)];
b=f;
maxiter=100;
tol=1e-6;
for k = 1:maxiter
    H=hess_diag(u,epsi,tau,Nv);
    vv=proxG(H,u-H^(-1)*grad(u,M,epsi,Nv,tau),b,A)-u;
    tk=1;
    while abs(imag((energyy(u+tk*vv,M,epsi,tau,dv,Nv))))>1e-15  || energyy(u+tk*vv,M,epsi,tau,dv,Nv)>energyy(u,M,epsi,tau,dv,Nv)+0.1*tk*grad(u,M,epsi,Nv,tau)'*vv  %|| norm(A*(u+tk*v)-b)>1e-6
        tk=0.5*tk;
    end
    u_next=u+tk*vv;
    if (sum(abs(u_next-u))/sum(abs(u))<tol) && (sum(abs(energyy(u_next,M,epsi,tau,dv,Nv)-energyy(u,M,epsi,tau,dv,Nv)))/sum(abs(energyy(u,M,epsi,tau,dv,Nv)))<tol)
        u=u_next;
        break
    end
    
    u=u_next;
    
    
end
[f_new,m]=decomp(u,Nv);
fv=mass*f_new;
end


function g=grad(u,M,epsi,Nv,tau)
[f,m]=decomp(u,Nv);
gE1=-epsi*m.^2./f.^2;
gE2=2*tau*(log(f./M)+ones(Nv,1));
g1=gE1+gE2;
g2=2*epsi*m./f;
g=[g1;g2];
end


function [F] = Fluc(fl, fr, fll, frr, v, dx, dt)
%fl,fr,fll,frr are sorted as fll, fl, fr, frr, like i-2, i-1, i, i+1, and
%the flucation we want is at i-1/2
% fl = [1 2 3];
% fr = [2 3 4];
% fll = [6 5 7];
% frr = [8 4 1];
% v = [2 -5 3];
% dx =1;
% dt =1;
theta = max(sign(v),0).*(fl-fll)./(fr-fl) - min(sign(v),0).*(frr-fr)./(fr-fl);
%vtheta=0*theta;
%for i = 1:length(theta)
%    vtheta(i)=max(0,min(2*theta(i),(theta(i)+1)/2,2));
%end
vtheta=max(0,min(2*theta,min((theta+1)/2,2)));
%F = 1/2*abs(v).*(1-dt/dx*abs(v)).*(fr-fl).*VanLeer(theta);
F = 1/2*abs(v).*(1-dt/dx*abs(v)).*(fr-fl).*vtheta;
end

function name=num2str_decimal(a)
s=num2str(a);
c='';
for i = 1:length(s)
    if s(i)=='0'
        c(i)='z';
    elseif s(i)=='.'
        c(i)='p';
    elseif s(i)=='-'
        c(i)='n';
    else
        c(i)=s(i);
    end
end
name=c;
end


