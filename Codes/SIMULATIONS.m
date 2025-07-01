clear
close all
clc

%% Backbone curves

clear
close all
clc
tic

% define parameters

om0=10;
al2_vect=linspace(9.488,9.64,8);
al3=1;
al4=0;
al5=0;

% iterate different values of alpha_2 parameter

aZD_points_sym=cell(1,length(al2_vect));
aZD_points=zeros(1,length(al2_vect));
omegaZD_points=zeros(1,length(al2_vect));
aZD_points_paper=zeros(1,length(al2_vect));
omegaZD_points_paper=zeros(1,length(al2_vect));

figure
for ii=1:length(al2_vect)
    
    al2=al2_vect(ii);

    % compute and plot the corresponding backbone curve
    g_eff = al3 - 10/9*(al2/om0)^2;
    s_eff = al5 - 11/12*(al2^2/om0^3)^2 + (53/20*al3)*(al2/om0^2)^2 - (14*al2*al4)/(5*om0^2) + (3/80)*(al3/om0)^2;
    om_res=@(a_res) om0+(3*g_eff)/(8*om0)*a_res.^2+(5*s_eff)/(16*om0)*a_res.^4;
    
    fplot(om_res,[0,2])
    hold on
    grid on
    xlabel('a_{res}')
    ylabel('\omega_{res}')
    title('Backbone curves for different values of \alpha_2')
    
    % compute the corresponding zero dispersion point and save it in a
    % vector
    syms a
    f=om0+(3*g_eff)/(8*om0)*a^2+(5*s_eff)/(16*om0)*a^4;
    dfda=diff(f,1,a);
    aZD_points_sym{ii}=solve(dfda==0,a);
    aZD_points(ii)=aZD_points_sym{ii}(3,1);
    omegaZD_points(ii)=om_res(aZD_points(ii));

    % compute the corresponding zero dispersion point (formula of the
    % paper) and save it in a vector
    aZD_points_paper(ii)=sqrt( (-3*g_eff)/(5*s_eff) );
    omegaZD_points_paper(ii)=om0-(9*g_eff^2)/(80*om0*s_eff);

end

% plot zero-dispersion points locus

zd_plot=plot(aZD_points,omegaZD_points,'k--','LineWidth',1.5);
ylim([9.999,10.005])
legend(zd_plot,'zero-dispersions points','Location','best')

% plot zero-dispersion points locus (formula in the paper)

plot(aZD_points_paper,omegaZD_points_paper,'b-.','LineWidth',1);

hold off

toc
%% Response curves
clear
close all
clc
tic

% define parameters

ep=1;
om0=10;
al2=9.64;
al3=1;
al4=0;
al5=0;
F=1e-2;

% compute and plot the backbone curve

g_eff = al3 - 10/9*(al2/om0)^2;
s_eff = al5 - 11/12*(al2^2/om0^3)^2 + (53/20*al3)*(al2/om0^2)^2 - (14*al2*al4)/(5*om0^2) + (3/80)*(al3/om0)^2;
om_res=@(a_res) om0+(3*g_eff)/(8*om0)*a_res.^2+(5*s_eff)/(16*om0)*a_res.^4;

a_backb=linspace(0,2,100);
om_backb=om_res(a_backb);

figure
h=plot(om_backb,a_backb,'k--','LineWidth',1.5);
hold on
grid on
title('Response curves for different values of \Gamma')
xlabel('\omega')
ylabel('a_{ss}')

% compute and plot response curves for different values of Gamma damping
% coefficient

Gamma=@(a_res) F/(2*om0*a_res);
a_res_vect=linspace(0.2,2,10);
om=linspace(9.995,10.005,1000);
colors=jet(length(a_res_vect));

for ii=1:length(a_res_vect)

    % compute parameters for the the ii-th response curve
    a_res=a_res_vect(ii);
    G=F/(2*om0*a_res);
    plot_color=colors(ii,:);
    
    % solve numerically the steady-state response equation for each omega
    % point (x axis)
    syms a_ss    
    for jj=1:length(om)
        eqn=a_ss^2-((F/(2*om0))^2*(G^2+(om0-om(jj)+(3*g_eff)/(8*om0)*a_ss^2+(5*s_eff)/(16*om0)*a_ss^4)^2)^-1);
        sol=vpasolve(eqn==0,a_ss); % gives many roots (we only want to plot the non-complex ones)
        for kk=1:length(sol)
            if isreal(sol(kk))
                plot(om(jj),abs(sol(kk)),'.','Color',plot_color)
            end
        end
    end

end
legend(h,'backbone curve','Location','best')
axis([9.995 10.005 0 2.1])
hold off

toc
%% Numerical simulation of the system
clear
close all
clc
tic

% system data

om0=10;
al2=9.64;
al3=1;
al4=0;
al5=0;
ep=1;
G=2.5e-4;
F=0.01;

% simulation parameters

tspan=[0 30000];
y0=[0.1 0]';

% solution with ode45

[t,y]=ode45(@equaz_moto,tspan,y0);
x=y(:,1);
xp=y(:,2);

% plot the steady-state solution

    % extract steady-state solution
    
    regime_start_index=floor(0.9*length(t));
    t_regime=t(regime_start_index:end);
    x_regime=x(regime_start_index:end);
    xp_regime=xp(regime_start_index:end);
    
    % time domain

    figure
    plot(t,x)
    grid on
    xlabel('t [s]')
    ylabel('x')
    title('Full system solution')

    figure
    plot(t_regime,x_regime)
    grid on
    xlabel('t [s]')
    ylabel('x')
    title('Steady-state solution')

    % spectrum
    
    X=fft(x_regime);
    Ts=mean(diff(t_regime));
    fsamp=1/Ts;
    N=length(t_regime); 
    f0=om0/(2*pi);

    P2=abs(X/N);
    P1=P2(1:N/2+1);
    P1(2:end-1)=2*P1(2:end-1);
    P1_dB=20*log10(P1);
    f=fsamp*(0:(N/2))/N;

    figure
    plot(f/f0,P1_dB,'k')
    grid on
    xlabel('\omega/\omega_0')
    ylabel('Spectrum amplitude [dB]')
    title('Steady state solution - frequency content')
    xlim([0 7])
    ylim([-200 20])

    % comparison with the theoretical harmonic amplitudes
    
    c02=-al2/(2*om0^2);
    c04=-((19*al2^3-45*om0^2*al2*al3+27*om0^4*al4) / (72*om0^6));
    c22=al2 / (6*om0^2);
    c24=-((14*al2^3+45*om0^2*al2*al3-48*om0^4*al4) / (288*om0^6));
    c33=(2*al2^2+3*om0^2*al3)/(96*om0^4);
    c35=-((580*al2^4+3240*om0^2*al2^2*al3-405*om0^4*al3^2+648*om0^4*al2*al4-2700*om0^6*al5) / (69120*om0^8));
    c44=(10*al2^3+45*om0^2*al2*al3+36*om0^4*al4) / (4320*om0^6);
    c55=(100*al2^4+900*om0^2*al3*al2^2+1584*om0^4*al2*al4+405*om0^4*al3^2+1080*al5*om0^6) / (414720*om0^8);
    
    a0=@(a) c02*(ep*a).^2+c04*(ep*a).^4;
    a2=@(a) c22*(ep*a).^2+c24*(ep*a).^4;
    a3=@(a) c33*(ep*a).^3+c35*(ep*a).^5;
    a4=@(a) c44*(ep*a).^4;
    a5=@(a) c55*(ep*a).^5;
    
    pks=findpeaks(P1_dB);
    a_res=max(pks);
    a_res=db2mag(a_res);

    a0_res=a0(a_res);
    a2_res=a2(a_res);
    a3_res=a3(a_res);
    a4_res=a4(a_res);
    a5_res=a5(a_res);

    hold on
    bottom=-500;
    p0=plot([0 0],[bottom mag2db(abs(a0_res))],'r-',LineWidth=3);
    p0.Color(4)=0.2;
    p1=plot([1 1],[bottom mag2db(abs(a_res))],'g-',LineWidth=3);
    p1.Color(4)=0.2;
    p2=plot([2 2],[bottom mag2db(abs(a2_res))],'b-',LineWidth=3);
    p2.Color(4)=0.2;
    p3=plot([3 3],[bottom mag2db(abs(a3_res))],'y-',LineWidth=3);
    p3.Color(4)=0.2;
    p4=plot([4 4],[bottom mag2db(abs(a4_res))],'c-',LineWidth=3);
    p4.Color(4)=0.2;
    p5=plot([5 5],[bottom mag2db(abs(a5_res))],'m-',LineWidth=3);
    p5.Color(4)=0.2;
    hold off

    a_axis=linspace(0,2,100);
    figure
    plot(a_axis,a0(a_axis)*1e1,'r','LineWidth',1.5);
    hold on
    plot(a_axis,a2(a_axis)*1e2,'b','LineWidth',1.5);
    plot(a_axis,a3(a_axis)*1e3,'y','LineWidth',1.5);
    plot(a_axis,a4(a_axis)*1e4,'c','LineWidth',1.5);
    plot(a_axis,a5(a_axis)*1e5,'m','LineWidth',1.5);
    xline(a_res,'k')
    grid on
    xlabel('a')
    ylabel('a_k')
    title('Amplitudes of overtones and DC')
    legend('a_0\cdot10','a_2\cdot10^2','a_3\cdot10^3','a_4\cdot10^4','a_5\cdot10^5','Location','northwest')
    hold off

    a_res_th=F/(2*om0*G);
    fprintf('Amplitude of the main harmonic computed with the hybrid method: a_{res}=%f\n',a_res_th)

    % phase space (for unforced case)

    % icond={y0};
    % Xlim=[-4 4];
    % Ylim=[-4 4];
    % figure
    % PhasePlane(@equaz_moto,tspan,icond,'Xlim',Xlim,'Ylim',Ylim,'hx',0.5,'hy',0.5,'scale',1,'PlotSingularPoints',true);

toc

%% FUNCTION DEFINITION: equation of motion

function dydt=equaz_moto(t,y)
    if nargin==1
           y=t;
    end
    
    G=2.5e-4;
    g1=0;
    g2=0;
    g3=0;
    g4=0;
    om0=10;
    al2=9.64;
    al3=1;
    al4=0;
    al5=0;
    h1=0;
    h2=0;
    F=0.01;
    Om=3e-3+om0;

    g=@(y) 1+g1*y(1)+g2*y(1).^2+g3*y(1).^3+g4*y(1).^4 + 0*y(2);
    f=@(y) al2*y(1).^2+al3*y(1).^3+al4*y(1).^4+al5*y(1).^5 + 0*y(2);
    h=@(y) 1+h1*y(1)+h2*y(1).^2 + 0*y(2);

    dydt=[y(2); -2*G.*g(y).*y(2)-om0^2*y(1)-f(y)+h(y)*F*cos(Om*t)];
end
