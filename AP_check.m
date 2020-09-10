function AP_check(Nx,Nv,Lx,Lv,tau,T)
%%%%%%%%%%%%% Check AP for 1D case %%%%%%
% This solver aims to solve
% \patial_t f + v \nabla_x f = 1/\epsi ( \nabla_x \phi \nabla_v f + FP(f)).
% -\Delta_x \phi = \rho -h 
% tau is the preset outer time step
% Nx Nv are the number of uniform grid point
% T is final time
% Author: Wuzhe Xu
% Date: 09/05/2020
% Code for figure 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NN=4;
legend_name = cell(1,NN);
for ii = 1:NN
    temp_epsi = 10^(-1-ii);
    PN_search_VPFP_1D(temp_epsi,Nx,Nv,Lx,Lv,tau,T);
    filename=['PN_vpfp_epsi', num2str_decimal(temp_epsi),'_t_'   ,  num2str_decimal(T)];
    load(filename)
    E2M_vec(ii,:) = E2M;
    temp_s = regexprep(cellstr(num2str(temp_epsi.', '%.0e')), '(?<=e[-+])0*', '');
    legend_name{ii}=['\epsilon =', temp_s{1}];
end
ll=length(E2M_vec)
xk=(1:1:ll);
semilogy(xk,E2M_vec(1,:),'ro',xk,E2M_vec(2,:),'b^',xk,E2M_vec(3,:),'gx',xk,E2M_vec(4,:),'mp','Linewidth',2)
title('1D VPFP: asymptotic behavior')
legend(legend_name)
xlabel('Time','Fontsize',25)
ylabel('log||f^n-M^n||_1','Fontsize',20)
xk_step = round(ll/10);
label_seq = round((1:xk_step:ll-1)*tau,3);
tick_seq = (1:xk_step:ll-1);
xticks(tick_seq);
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