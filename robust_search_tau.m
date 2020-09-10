function robust_search_tau(Nx,L,IC)
%%%%%%%%%%%%% Robustness to tau %%%%%%
% Use line search
% tau is the outer time step
% Nx is the number of uniform grid point
% T is final time
% v \in [-L,L]
% Use IC 3
% Author: Wuzhe Xu
% Date: 09/05/2020
% Code for Figure 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error_u_epsi_1=[];
error_u_epsi_n5=[];
KKK=50;
NN=6;
legend_name = cell(1,NN);
for ii = 1:NN
    temp_tau = 0.4/2^(ii-1);
    PN_1D_search_end(1,Nx,L,temp_tau,temp_tau,IC);
    filename = ['PN_search_tau_', num2str_decimal(temp_tau), '_T_', num2str_decimal(temp_tau), '_epsi_',num2str_decimal(1),'_end'];
    load(filename)
    u_star=u;
    for j = 1:KKK
        error_u_epsi_1(ii,j) = sum(abs(u_com(:,j)-u_star))/sum(abs(u_star));
    end
    legend_name{ii}=['\tau =', num2str(temp_tau)];
end
for ii = 1:NN
    temp_tau = 0.4/2^(ii-1);
    PN_1D_search_end(1e-5,Nx,L,temp_tau,temp_tau,IC);
    filename = ['PN_search_tau_', num2str_decimal(temp_tau), '_T_', num2str_decimal(temp_tau), '_epsi_',num2str_decimal(1e-5),'_end'];
    load(filename)
    u_star=u;
    for j = 1:KKK
        error_u_epsi_n5(ii,j) = sum(abs(u_com(:,j)-u_star))/sum(abs(u_star));
    end
    legend_name{ii}=['\tau =', num2str(temp_tau)];
end
xk=linspace(1,KKK,KKK);
figure(1)
semilogy(xk,error_u_epsi_1(1,:),'k-p',xk,error_u_epsi_1(2,:),'r-o',xk,error_u_epsi_1(3,:),'b.-',xk,error_u_epsi_1(4,:),'c-*',xk,error_u_epsi_1(5,:),'m-^',xk,error_u_epsi_1(6,:),'g--','Linewidth',2)
title('1D convergence rate with \epsilon =1')
%legend('\tau=0.0125','\tau=0.025','\tau=0.05','\tau=0.1','\tau=0.2','\tau=0.4')
legend(legend_name)
xlabel('Iterations','Fontsize',25)
ylabel('error_k','Fontsize',25)
ylim([1e-22,2])
set(gca,'FontSize',30)
set(gcf,'position',[1,1,1440,900])
figure(2)
semilogy(xk,error_u_epsi_n5(1,:),'k-p',xk,error_u_epsi_n5(2,:),'r-o',xk,error_u_epsi_n5(3,:),'b.-',xk,error_u_epsi_n5(4,:),'c-*',xk,error_u_epsi_n5(5,:),'m-^',xk,error_u_epsi_n5(6,:),'g--','Linewidth',2)
title('1D convergence rate with \epsilon =10^{-5}')
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
