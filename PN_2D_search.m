function PN_2D_search(epsi,L,Nx,Ny,tau,T,IC)
%%%%%%%%%%%%% Proximal Quasi newton method for 2D scaled Fokker planck %%%%%%
% This solver aims to solve
% \patial_t f = 1/\epsi FP(f).
% Stiffness comes from 1/\epsi
% tau is the outer time step
% Nx Ny are the number of uniform grid point
% T is final time
% IC stand for initial condition, pick from 1-4
%%%%%%%%%%%%%
warning off
close all
%% Initial
dx=2*L/Nx;
dy=2*L/Ny;
x=(-L+dx/2:dx:L);
y=(-L+dy/2:dy:L);
N=Nx*Ny;
[X,Y]=meshgrid(x,y);
switch IC
    case 1
        f0 = exp(-((X-1).^2+(Y-1).^2));
        disp('One bump to the upper left')
    case 2
        f0=exp(-((X+1).^2+(Y+2).^2))+0.5*exp(-((X-1).^2+(Y-1).^2));
        disp('Two asymmetric bumps')
    case 3
        f0=exp(-((X+2).^2+(Y+2).^2))+exp(-((X-2).^2+(Y-2).^2));
        disp('Two symmetric bumps')
    case 4
        f=zeros(Nx,Ny);
        for i=1:Nx
            for j = 1:Ny
                f(i,j)=1.5*(1+((((x(i)-2)^2+(y(j)-2)^2))^0.5-2)^2)^(-10);
            end
        end
        for i=1:Nx
            for j = 1:Ny
                f(i,j)=f(i,j)+2*(1+((((x(i)+2)^2+(y(j)+2)^2))^0.5-2)^2)^(-10);
            end
        end
        f0=f;
        disp('Two torus')
    otherwise
        error('Case not available, choose from {1,2,3,4}')
end
Mas=sum(sum(f0))*dx*dy;
M = 1/2/pi*exp(-((X).^2+(Y).^2)/2);
M = M/(sum(sum(M))*dx*dy)*Mas;
%% prepare matrix A
A_rho=eye(N);
R_x=eye(Nx);
L_x=-eye(Nx);
structure_R=zeros(Nx,Nx);
structure_L=zeros(Nx,Nx);
for i = 2:Nx-1
    structure_R(i,i+1) = 1;
    structure_L(i,i-1) = 1;
end
structure_R(1,2) = 1;
structure_R(1,1) = 1;
structure_L(end,end-1) = 1;
structure_L(end,end) = 1;
A_mx=(kron(structure_L,L_x)+kron(structure_R,R_x))./2./dx;
DY=zeros(Ny,Ny);
for i = 2:Ny-1
    DY(i,i+1)=1;
    DY(i,i-1)=-1;
end
DY(1,2)=1;
DY(1,1)=1;
DY(end,end-1)=-1;
DY(end,end)=-1;
A_my=kron(eye(Ny),DY)./2./dy;
A=[A_rho A_mx A_my];
t=0;
%% stack matrix as vector
f_vec=m2v(f0);
m0=zeros(2*N,1);
u0=[f_vec;m0];
b=f_vec;
%% set up parameter, stepsize
tol=1e-6;
maxiter=200;
KK=180;
u=u0;
u_com(:,1)=u;
f_com(:,1)=f_vec;
u=u0;
number=0;
while t<T
    for k=1:maxiter   
        %% proximal quasi newton
        H=hess_2D_diag(u,epsi,tau,N);
        v=proxG(H,u-H^(-1)*grad_2D(u,N,epsi,tau,M),A,b)-u;
        tk=1;
        while abs(imag((energyy(u+tk*v,dx,dy,tau,M,epsi,N))))>1e-20  || energyy(u+tk*v,dx,dy,tau,M,epsi,N)>energyy(u,dx,dy,tau,M,epsi,N)+0.01*tk*grad_2D(u,N,epsi,tau,M)'*v  %|| norm(A*(u+tk*v)-b)>1e-6
            tk=0.5*tk;
        end
        u_next=u+tk*v;
        error_E = abs(energyy(u_next,dx,dy,tau,M,epsi,N)-energyy(u,dx,dy,tau,M,epsi,N))/abs(energyy(u,dx,dy,tau,M,epsi,N));
        error_u = abs(norm(u_next)-norm(u))/abs(norm(u));   
        if (error_E <=tol && error_u<=tol)
            u=u_next;
            break
        end
        u=u_next;
        [f,mx,my] = decomp(u,N);
        %mxx=v2m(mx,Nx,Ny);
        %myy=v2m(my,Nx,Ny);
        %Mf=v2m(f);
        %mesh(X,Y,mxx./Mf)
        %figure(2)
        %mesh(X,Y,myy./Mf)
        %drawnow
        if min(f)<0
            break
        end  
        if k<KK
            u_com(:,k+1) = u;
            f_com(:,k+1) = f;
        end
        %figure(2)
        %mesh(X,Y,Mf)
        %drawnow
        
    end
    [f,mx,my] = decomp(u,N);
    Mf=v2m(f,Nx,Ny);
    mesh(X,Y,Mf)
    drawnow
    t=t+tau
    number=number+1;
    b=f;
    %filename=['PN_2D_search_sm25_2n1_ex_torus_',num2str(number)];
    %save(filename)
end

end

function [f,mx,my]=decomp(u,N)
f=u(1:N);
mx=u(N+1:2*N);
my=u(2*N+1:end);
end


function v=m2v(F)
v=reshape(F',1,[])';        
end

function m=v2m(f,Nx,Ny)
m=reshape(f,Nx,Ny);
end

function E=energyy(u,dx,dy,tau,M,epsi,N)
[f,mx,my]=decomp(u,N);
En1=(mx.^2+my.^2)./f;
E1=epsi*sum(sum(En1))*dx*dy;
E2=sum(f.*(log(f./m2v(M))))*dx*dy;
E=E1+2*tau*E2;
end


function gE=grad_2D(u,N,epsi,tau,M)
gE1=zeros(3*N,1);
[f,mx,my]=decomp(u,N);
gE1(1:N)=-epsi*(mx.^2+my.^2)./f.^2+2*tau*(log(f./m2v(M)) + ones(N,1));
gE1(N+1:2*N) = 2*epsi*(mx)./f;
gE1(2*N+1:end) = 2*epsi*(my)./f;
gE=gE1;
end


function Hd=hess_2D_diag(u,epsi,tau,N)
[f,mx,my]=decomp(u,N);
%f=max(f,5e-4);
H1=diag(2*epsi*(mx.^2+my.^2)./f.^3+2*tau./f);
H2=diag(-2*epsi*mx./f.^2)*0;
H3=diag(-2*epsi*my./f.^2)*0;
H4=diag(2*epsi./f);
[n1,n2]=size(H4);
Hd=[H1 H2 H3; H2 H4 zeros(n1,n1); H3 zeros(n1,n1) H4];
end



function p=proxG(H,u,A, b)
lambda=(A*H^(-1)*A')^(-1)*(b-A*u);
p=u+H^(-1)*A'*lambda;
end

