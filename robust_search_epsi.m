function robust_search_epsi(Nx,L,tau,IC)
%%%%%%%%%%%%% Robustness to epsi %%%%%%
% Use line search
% tau is the outer time step
% Nx is the number of uniform grid point
% T is final time
% v \in [-L,L]
% Use IC 3
% Author: Wuzhe Xu
% Date: 09/05/2020
% Code for Figure 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error_u=[];
error_f=[];
KKK=50;
NN=6;
legend_name = cell(1,NN);
for ii = 1:NN
    temp_eps = 10^(1-ii);
    
    PN_1D_search_end(temp_eps,Nx,L,tau,tau,IC);
    filename = ['PN_search_tau_', num2str_decimal(tau), '_T_', num2str_decimal(tau), '_epsi_',num2str_decimal(temp_eps),'_end'];
    load(filename)
    u_star=u;
    f_star=f;
    for j = 1:KKK
        error_u(ii,j) = sum(abs(u_com(:,j)-u_star))/sum(abs(u_star));
        error_f(ii,j) = sum(abs(f_com(:,j)-f_star))/sum(abs(f_star));
    end
    if temp_eps==1
        legend_name{ii}=['\epsilon =', num2str(1)];
    else
        temp_s = regexprep(cellstr(num2str(temp_eps.', '%.0e')), '(?<=e[-+])0*', '');
        legend_name{ii}=['\epsilon =', temp_s{1}];
    end
end
xk=linspace(1,KKK,KKK);
figure(1)
semilogy(xk,error_u(1,:),'k-p',xk,error_u(2,:),'r-o',xk,error_u(3,:),'b.-',xk,error_u(4,:),'c-*',xk,error_u(5,:),'m-^',xk,error_u(6,:),'g--','Linewidth',2)
title('1D convergence rate')
%legend('\tau=0.0125','\tau=0.025','\tau=0.05','\tau=0.1','\tau=0.2','\tau=0.4')
legend(legend_name)
xlabel('Iterations','Fontsize',25)
ylabel('error_k','Fontsize',25)
ylim([1e-22,2])
set(gca,'FontSize',30)
set(gcf,'position',[1,1,1440,900])
figure(2)
semilogy(xk,error_f(1,:),'k-p',xk,error_f(2,:),'r-o',xk,error_f(3,:),'b.-',xk,error_f(4,:),'c-*',xk,error_f(5,:),'m-^',xk,error_f(6,:),'g--','Linewidth',2)
title('1D convergence rate')
%legend('\tau=0.0125','\tau=0.025','\tau=0.05','\tau=0.1','\tau=0.2','\tau=0.4')
legend(legend_name)
xlabel('Iterations','Fontsize',25)
ylabel('error_k','Fontsize',25)
ylim([1e-22,2])
set(gca,'FontSize',30)
set(gcf,'position',[1,1,1440,900])

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
