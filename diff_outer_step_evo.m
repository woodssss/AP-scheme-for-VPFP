function diff_outer_step_evo(Nv,L,tau,T)
%%%%%%%%%%%%% test convergence at different time steps %%%%%%
% Use line search
% tau is the outer time step
% Nx is the number of uniform grid point
% T is final time
% v \in [-L,L]
% Use IC 3
% Author: Wuzhe Xu
% Date: 09/05/2020
% Code for Figure 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NN=6;
for ii=1:NN
    temp_epsi = 10^(1-ii);
    if ii<=3
        lambda=0.5;
    else
        lambda=0.4;
    end
    PN_1D_evo(temp_epsi,Nv,L,lambda,tau,T,3)
    filename = ['PN_tau_', num2str_decimal(tau), '_epsi_',num2str_decimal(temp_epsi),'_evo'];
    load(filename)
    num_iter_epsi(ii,:) = numer_iter;
    if temp_epsi==1
        legend_name{ii}=['\epsilon =', num2str(1)];
    else
        temp_s = regexprep(cellstr(num2str(temp_epsi.', '%.0e')), '(?<=e[-+])0*', '');
        legend_name{ii}=['\epsilon =', temp_s{1}];
    end
end
ll=size(num_iter_epsi,2);
xk=(1:1:ll);
figure(1)
%semilogy(xk,num_iter_epsi(1,:),'ro',xk,num_iter_epsi(2,:),'b^',xk,num_iter_epsi(3,:),'gx',xk,num_iter_epsi(4,:),'k-p',xk,num_iter_epsi(5,:),'m-^',xk,num_iter_epsi(6,:),'g--','Linewidth',2)
plot(xk,num_iter_epsi(1,:),'r-o',xk,num_iter_epsi(2,:),'b-^',xk,num_iter_epsi(3,:),'g-x',xk,num_iter_epsi(4,:),'k-p',xk,num_iter_epsi(5,:),'m-^',xk,num_iter_epsi(6,:),'g--','Linewidth',2)
title('Convergence behavior at different t')
%legend('\tau=0.0125','\tau=0.025','\tau=0.05','\tau=0.1','\tau=0.2','\tau=0.4')
legend(legend_name)
xk_step = round(ll/10);
label_seq = round((1:xk_step:ll)*tau,3);
tick_seq = (1:xk_step:ll);
xticks(tick_seq);
xlabel('Time','Fontsize',25)
ylabel('Iteration till convergence','Fontsize',25)
%ylim([1e-22,2])
set(gca, 'xticklabel', label_seq);
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