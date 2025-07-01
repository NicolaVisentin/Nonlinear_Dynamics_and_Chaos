%% F critic condition

clear
close all
clc

% define parameters

om0=10;
al3=1;
al4=0;
al5=0;
G=2.5e-4;

% compute alpha2 coefficient such that sigma_effective=0

syms al2

s_eff = al5 - 11/12*(al2^2/om0^3)^2 + (53/20*al3)*(al2/om0^2)^2 - (14*al2*al4)/(5*om0^2) + (3/80)*(al3/om0)^2;
al2=solve(s_eff==0,al2);
for ii=1:length(al2)
    if isreal(al2(ii))
        al2=al2(ii);
        break
    end
end

s_eff=eval(s_eff);
al2=double(al2);
s_eff=double(s_eff);
g_eff = al3 - 10/9*(al2/om0)^2;

% compute F_cr and the force vector

F_cr=sqrt( (256*om0^2*G^3) / (9*sqrt(3)*abs(g_eff)) );

F_vect=linspace(0.1*F_cr,10*F_cr,10);
F_cr_indx=find(F_vect>F_cr,1);
F_vect(F_cr_indx)=F_cr;

% compute and plot the backbone curve

om_res=@(a_res) om0+(3*g_eff)/(8*om0)*a_res.^2+(5*s_eff)/(16*om0)*a_res.^4;

a_max=F_vect(end)/(2*om0*G);
a_backb=linspace(0,1.1*a_max,100);
om_backb=om_res(a_backb);

figure
h=plot(om_backb,a_backb,'k--','LineWidth',1.5);
hold on
grid on
xlabel('\omega')
ylabel('a_{ss}')

% compute and plot response curves for different values of the force F

F_vect=linspace(F_cr/10,F_cr*10,10);
F_vect(find(F_vect>F_cr,1))=F_cr;
om=linspace(9.995,10.005,50);
colors=jet(length(F_vect));

indx_sn=zeros(1,length(F_vect));
for ii=1:length(F_vect)
    count=0;

    % compute parameters for the the ii-th response curve
    F=F_vect(ii);
    plot_color=colors(ii,:);
    
    % solve numerically the steady-state response equation for each omega
    % point (x axis)
    syms a_ss    
    for jj=1:length(om)
        eqn=a_ss^2-((F/(2*om0))^2*(G^2+(om0-om(jj)+(3*g_eff)/(8*om0)*a_ss^2+(5*s_eff)/(16*om0)*a_ss^4)^2)^-1);
        sol=vpasolve(eqn==0,a_ss); % gives many roots (real and complex)
        for kk=1:length(sol)
            if isreal(sol(kk))
                plot(om(jj),abs(sol(kk)),'.','Color',plot_color)    % we want to plot only the real ones...
                count=count+1;                                      % ...and check if there is more than one real solution for some omega (saddle-node)
            end
        end
    end
    if count/2~=length(om)  % actually there are 2 equal real solutions
        indx_sn(ii)=1;
    end

end
legend(h,'backbone curve','Location','best')
title('Response curves for different values of \Gamma')
axis([9.995 10.005 0 a_max])
hold off

% print info

fprintf('alpha2=%f \nsigma_eff=%f \nF_cr=%f \n\nF_cr corresponds to the curve n. %d \nSaddle-node is found from curve n. %d\n\n', al2,s_eff,F_cr,F_cr_indx,find(indx_sn==1,1))
