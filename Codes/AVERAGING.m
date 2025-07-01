clear
close all
clc

%% Symbolic variables definition
tic

syms ep Om t Phi Phip om0 
syms a a0 a2 a3 a4 a5 ap
syms al2 al3 al4 al5
syms x xp xpp
syms c02 c04 c22 c24 c33 c35 c44 c55
syms Phi Phip
syms G g1 g2 g3 g4
syms h1 h2 F

%% Definition of ak coefficients

c02=-al2/(2*om0^2);
c04=-((19*al2^3-45*om0^2*al2*al3+27*om0^4*al4) / (72*om0^6));
c22=al2 / (6*om0^2);
c24=-((14*al2^3+45*om0^2*al2*al3-48*om0^4*al4) / (288*om0^6));
c33=(2*al2^2+3*om0^2*al3)/(96*om0^4);
c35=-((580*al2^4+3240*om0^2*al2^2*al3-405*om0^4*al3^2+648*om0^4*al2*al4-2700*om0^6*al5) / (69120*om0^8));
c44=(10*al2^3+45*om0^2*al2*al3+36*om0^4*al4) / (4320*om0^6);
c55=(100*al2^4+900*om0^2*al3*al2^2+1584*om0^4*al2*al4+405*om0^4*al3^2+1080*al5*om0^6) / (414720*om0^8);

a0=c02*(ep*a)^2+c04*(ep*a)^4;
a0=expand(a0);
a2=c22*(ep*a)^2+c24*(ep*a)^4;
a2=expand(a2);
a3=c33*(ep*a)^3+c35*(ep*a)^5;
a3=expand(a3);
a4=c44*(ep*a)^4;
a4=expand(a4);
a5=c55*(ep*a)^5;
a5=expand(a5);

%% Definition of x, xp, xpp

x = ep*a*cos(Om*t+Phi) + a0 + a2*cos(2*(Om*t+Phi)) + a3*cos(3*(Om*t+Phi)) + a4*cos(4*(Om*t+Phi)) + a5*cos(5*(Om*t+Phi));
x=expand(x);

xp = - ep*Om*a*sin(Om*t+Phi) - 2*Om*a2*sin(2*(Om*t+Phi)) - 3*Om*a3*sin(3*(Om*t+Phi)) - 4*Om*a4*sin(4*(Om*t+Phi)) - 5*Om*a5*sin(5*(Om*t+Phi));
xp=expand(xp);

xpp = - ep^5*Om*ap*sin(Om*t+Phi) - ep^5*Om*(Om/ep^4+Phip)*a*cos(Om*t+Phi) - 4*Om^2*a2*cos(2*(Om*t+Phi)) - 9*Om^2*a3*cos(3*(Om*t+Phi)) ...
     - 16*Om^2*a4*cos(4*(Om*t+Phi)) - 25*Om^2*a5*cos(5*(Om*t+Phi));
xpp=expand(xpp);

%% Equation of motion

expr=[xpp, 2*G*xp, om0^2*x, al2*x^2, al3*x^3, al4*x^4, al5*x^5, -F*cos(Om*t)];

expanded_expr=sym(zeros(1,length(expr)));
parfor ii=1:length(expr)
    expanded_expr(ii)=expand(expr(ii));
end

espressione_espansa=sum(expanded_expr);

%% Remove order 6 or higher in epsilon 

comp=children(espressione_espansa);
logic_6=zeros(1,size(comp,2));
parfor ii=1:size(comp,2)
    logic_6(ii)=(has(comp{ii},ep^6) || has(comp{ii},ep^7) || has(comp{ii},ep^8) || has(comp{ii},ep^9) || has(comp{ii},ep^10) || has(comp{ii},ep^11) || has(comp{ii},ep^12) || has(comp{ii},ep^13) || has(comp{ii},ep^14) || has(comp{ii},ep^15) || has(comp{ii},ep^16) || has(comp{ii},ep^17) || has(comp{ii},ep^18) || has(comp{ii},ep^19) || has(comp{ii},ep^20) || has(comp{ii},ep^21) || has(comp{ii},ep^22) || has(comp{ii},ep^23) || has(comp{ii},ep^24) || has(comp{ii},ep^25));
end
idx=find(~logic_6);

comp_av=sym(0);
for ii=1:length(idx)
    comp_av=comp_av+comp{idx(ii)};
end

%% Solve for ap and Phip, together with the constaint equation

espr_new=comp_av;
const_eqn=ap*cos(Om*t+Phi)-a*Phip*sin(Om*t+Phi);

[ap_sol,Phip_sol]=solve([espr_new==0, const_eqn==0],[ap,Phip]);

ap_sol=expand(ap_sol);
Phip_sol=expand(Phip_sol);

%% Average over a period

T=2*pi/Om;

comp_ap=children(ap_sol);
comp_Phip=children(Phip_sol);

comp_ap_av=sym(0);
for ii=1:size(comp_ap,2)
    comp_ap_av=comp_ap_av+int(comp_ap{ii},t,0,T);
end

comp_Phip_av=sym(0);
for ii=1:size(comp_Phip,2)
    comp_Phip_av=comp_Phip_av+int(comp_Phip{ii},t,0,T);
end

comp_ap_av=(comp_ap_av)/T;
comp_Phip_av=(comp_Phip_av)/T;

comp_ap_av=simplify(comp_ap_av);
comp_Phip_av=simplify(comp_Phip_av);

comp_ap_av=expand(comp_ap_av);
comp_Phip_av=expand(comp_Phip_av);

%% Evaluate for epsilon=1

ep=1;

comp_ap_av=eval(comp_ap_av)
comp_Phip_av=eval(comp_Phip_av)

toc
%% Saving workspace
tic

timestamp=datetime('now','TimeZone','local','Format','d_MMM_y__HH_mm_ss');
timestamp=string(timestamp);
save(['workspace_averaging__', num2str(timestamp), '.mat'])

toc