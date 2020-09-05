function PN_3D_search_evo(epsi,L,Nx,Ny,Nz,tau,T)
%%%%%%%%%%%%% Proximal Quasi newton method for 3D scaled Fokker planck %%%%%%
% This solver aims to solve
% \patial_t f = 1/\epsi FP(f).
% Stiffness comes from 1/\epsi
% tau is the outer time step
% Nx Ny Nz are the number of uniform grid point
% T is final time
% Author: Wuzhe Xu
% Date: 09/05/2020
% Code of paper VARIATIONAL ASYMPTOTIC PRESERVING SCHEME FORVLASOV-POISSON-FOKKER-PLANCK SYSTEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off
close all
%% Initial
dx=2*L/Nx;
dy=2*L/Ny;
dz=2*L/Nz;
x=(-L+dx/2:dx:L);
y=(-L+dy/2:dy:L);
z=(-L+dz/2:dz:L);
N=Nx*Ny*Nz;
[X,Y]=meshgrid(x,y);
[YY,ZZ]=meshgrid(y,z);
% Define IC
f0_xy=1/2/pi*(exp(-((X-1).^2+(Y+1).^2))+exp(-((X+1).^2+(Y-1).^2)));
f0_xy=f0_xy/(sum(sum(f0_xy))*dx*dy);
for i = 1:Nz
    F(:,:,i)=exp(-z(i)^2/2)*f0_xy/(2*pi)^0.5;
end
Mxy = 1/2/pi*exp(-((X).^2+(Y).^2)/2);
M=[];
for i = 1:Nz
    M(:,:,i)=Mxy*exp(-z(i)^2/2)/(2*pi)^0.5;
end
error=0;
for i = 1:Nz
    error=error+norm(F(:,:,i)-M(:,:,i));
end
f0=m2v(F);
%% prepare matrix A
A_rho=eye(N);
DY=zeros(Ny,Ny);
for i = 2:Ny-1
    DY(i,i+1)=1;
    DY(i,i-1)=-1;
end
DY(1,2)=1;
DY(1,1)=1;
DY(end,end-1)=-1;
DY(end,end)=-1;
A_my=kron(eye(Nx*Nz),DY)./2./dy;
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
A_mx_temp=(kron(structure_L,L_x)+kron(structure_R,R_x))./2./dx;
A_mx=kron(eye(Nz),A_mx_temp);
structure=zeros(Nz,Nz);
for i=2:Nz-1
    structure(i,i+1)=1;
    structure(i,i-1)=-1;
end
structure(1,1)=1;
structure(1,2)=1;
structure(end,end)=-1;
structure(end,end-1)=-1;
A_mz=kron(structure,eye(Nx*Ny))./2./dz;
A=[A_rho A_mx A_my A_mz];

%% stack matrix as vector
f_vec=f0;
b=f_vec;
%% set up parameter
maxiter=2e2;
f=f_vec;
tol=8e-6;
t=0;
while t<T
    index=5;
    u=[f;zeros(3*N,1)];
    Fyz=zeros(Ny,Nz);
    Myz=zeros(Ny,Nz);
    for i = 1:Nz
        Fyz(:,i)=F(:,index,i);
        Myz(:,i)=M(:,index,i);
    end
    for k=1:maxiter 
        %% PN method
        HI=hess_3D_diag_inverse(u,epsi,tau,N);
        tk=1;
        v=proxG(HI,u-HI*grad_3D(u,N,epsi,tau,M),A,b)-u;
        G = grad_3D(u,N,epsi,tau,M);
        while abs(imag((energyy(u+tk*v,dx,dy,dz,tau,M,epsi,N))))>1e-30  || energyy(u+tk*v,dx,dy,dz,tau,M,epsi,N)>energyy(u,dx,dy,dz,tau,M,epsi,N)+0.05*tk*G'*v  || norm(A*(u+tk*v)-b)>1e-6
            tk=0.5*tk;
        end
        u_next=u+tk*v;
        
        [f,mx,my,mz] = decomp(u_next,N);
        FF=v2m(f,Nx,Ny,Nz);
        for i = 1:Nz
            error=error+norm(FF(:,:,i)-M(:,:,i));
        end
        if min(f)<0
            min(f)
            break
        end
        error_E = abs(energyy(u_next,dx,dy,dz,tau,M,epsi,N)-energyy(u,dx,dy,dz,tau,M,epsi,N))/abs(energyy(u,dx,dy,dz,tau,M,epsi,N));
        error_u = abs(norm(u_next)-norm(u))/abs(norm(u));
        
        if (error_E <=tol && error_u<=tol)
            u=u_next;
            break
        end
        u=u_next;
    end
    [f,mx,my,mz] = decomp(u,N);
    b=f;
    FF=v2m(f,Nx,Ny,Nz);
    figure(1)
    mesh(X,Y,FF(:,:,index))
    title(['f slice of z=', num2str(z(index)) ])
    drawnow
    
    Fyz=zeros(Ny,Nz);
    for i = 1:Nz
        Fyz(:,i)=FF(:,index,i);
    end
    figure(2)
    mesh(X,Y,Fyz)
    title(['f slice of x=', num2str(x(index)) ])
    drawnow
    t=t+tau
    %filename=['PN_3D_hess_new_2n1_',num2str(number)];
    %save(filename)
end 
end



function gE=grad_3D(u,N,epsi,tau,M)
gE1=zeros(4*N,1);
[f,mx,my,mz]=decomp(u,N);
f=max(f,1e-10);
gE1(1:N)=-epsi*(mx.^2+my.^2+mz.^2)./f.^2+2*tau*(log(f./m2v(M)) + ones(N,1));
gE1(N+1:2*N) = 2*epsi*(mx)./f;
gE1(2*N+1:3*N) = 2*epsi*(my)./f;
gE1(3*N+1:end) = 2*epsi*(mz)./f;
gE=gE1;
end

function [f,mx,my,mz]=decomp(u,N)
f=u(1:N);
mx=u(N+1:2*N);
my=u(2*N+1:3*N);
mz=u(3*N+1:end);
end

function v=m2v(F)
v=reshape(F,1,[])';        
end

function m=v2m(f,Nx,Ny,Nz)
m=reshape(f,Nx,Ny,Nz);
end

function p=proxG(HI,u,A, b)
lambda=(A*HI*A')^(-1)*(b-A*u);
p=u+HI*A'*lambda;
end

function Hd=hess_3D_diag_inverse(u,epsi,tau,N)
[f,mx,my,mz]=decomp(u,N);
H1=diag(1./(2*epsi*(mx.^2+my.^2+mz.^2)./f.^3+2*tau./f));
H2=diag(-2*epsi*mx./f.^2)*0;
H0=0*H2;
H3=diag(-2*epsi*my./f.^2)*0;
H4=diag(-2*epsi*mz./f.^2)*0;
H5=diag(1./(2*epsi./f));
Hd=[H1 H2 H3 H4; H0 H5 H0 H0; H0 H0 H5 H0; H0 H0 H0 H5];
end


function E=energyy(u,dx,dy,dz,tau,M,epsi,N)
[f,mx,my,mz]=decomp(u,N);
En1=(mx.^2+my.^2+mz.^2)./f;
E1=epsi*sum(sum(En1))*dx*dy*dz;
E2=sum(f.*(log(f./m2v(M))))*dx*dy*dz;
E=E1+2*tau*E2;
end
