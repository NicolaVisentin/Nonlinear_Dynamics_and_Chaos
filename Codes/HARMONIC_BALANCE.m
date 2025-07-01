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

%% Definition of ak coefficients

a0=c02*(ep*a)^2+c04*(ep*a)^4;
a2=c22*(ep*a)^2+c24*(ep*a)^4;
a3=c33*(ep*a)^3+c35*(ep*a)^5;
a4=c44*(ep*a)^4;
a5=c55*(ep*a)^5;

%% Definition of x, xp, xpp

x = ep*a*cos(Om*t+Phi) + a0 + a2*cos(2*(Om*t+Phi)) + a3*cos(3*(Om*t+Phi)) + a4*cos(4*(Om*t+Phi)) + a5*cos(5*(Om*t+Phi));
x=expand(x);

xp = - ep*Om*a*sin(Om*t+Phi) - 2*Om*a2*sin(2*(Om*t+Phi)) - 3*Om*a3*sin(3*(Om*t+Phi)) - 4*Om*a4*sin(4*(Om*t+Phi)) - 5*Om*a5*sin(5*(Om*t+Phi));
xp=expand(xp);

xpp = - ep^5*Om*ap*sin(Om*t+Phi) - ep^5*Om*(Om/ep^4+Phip)*a*cos(Om*t+Phi) - 4*Om^2*a2*cos(2*(Om*t+Phi)) - 9*Om^2*a3*cos(3*(Om*t+Phi)) ...
     - 16*Om^2*a4*cos(4*(Om*t+Phi)) - 25*Om^2*a5*cos(5*(Om*t+Phi));
xpp=expand(xpp);

%% Equation of motion

expr=[xpp, om0^2*x, al2*x^2, al3*x^3, al4*x^4, al5*x^5];
expanded_expr=sym(zeros(1,length(expr)));

parfor ii = 1:length(expr)
    expanded_expr(ii)=expand(expr(ii));
end

espressione_espansa = sum(expanded_expr);
espressione_semplificata=simplify(espressione_espansa,'Steps',1);

%% Harmonic extraction

% k=2,3,4,5 cases

comp=children(espressione_semplificata);
k=[2 3 4 5];

parfor jj=1:length(k) 

    logic_cos = false(1, size(comp, 2));
    for ii=1:size(comp,2)
        logic_cos(ii)=has(comp{ii}, cos(k(jj)*Phi + k(jj)*Om*t));
    end
    idx=find(logic_cos);

    ki_current=sym(0);
    for ii=1:length(idx)
        ki_current=ki_current+comp{idx(ii)}/cos(k(jj)*Phi + k(jj)*Om*t);
    end
    ki{jj}=sym(0);
    ki{jj}=simplify(ki_current,'Steps',1);

end

% k=0 case

logic_cos = false(1, size(comp, 2));
logic_sin = false(1, size(comp, 2));
for ii=1:size(comp,2)
    logic_cos(ii)=has(comp{ii}, 'cos');
    logic_sin(ii)=has(comp{ii}, 'sin');
end
logic=logic_sin+logic_cos;
idx=find(~logic);

ki_current=sym(0);
for ii=1:length(idx)
    ki_current=ki_current+comp{idx(ii)};
end
ki{5}=sym(0);
ki{5}=simplify(ki_current,'Steps',1);

%% Epsilon powers extraction

% eliminate all high powers of epsilon

for jj=1:size(ki,2)
    comp=children(expand(ki{jj}));
    logic_6=zeros(1,size(comp,2));
    parfor ii=1:size(comp,2)
        logic_6(ii)=(has(comp{ii},ep^6) || has(comp{ii},ep^7) || has(comp{ii},ep^8) || has(comp{ii},ep^9) || has(comp{ii},ep^10) || has(comp{ii},ep^11) || has(comp{ii},ep^12) || has(comp{ii},ep^13) || has(comp{ii},ep^14) || has(comp{ii},ep^15) || has(comp{ii},ep^16) || has(comp{ii},ep^17) || has(comp{ii},ep^18) || has(comp{ii},ep^19) || has(comp{ii},ep^20) || has(comp{ii},ep^21) || has(comp{ii},ep^22) || has(comp{ii},ep^23) || has(comp{ii},ep^24) || has(comp{ii},ep^25));
    end
    idx=find(~logic_6);
    
    comp_av=sym(0);
    for ii=1:length(idx)
        comp_av=comp_av+comp{idx(ii)};
    end
    ki{jj}=expand(comp_av);
end

% extract corresponding epsilon powers

c=[2 3 4 5];
for ii=1:size(ki,2)
    comp_eps=children(ki{ii});

    for jj=1:length(c)
        epj=ep^c(jj);

        logic_eps = false(1, size(comp_eps, 2));
        parfor kk=1:size(comp_eps,2)
            logic_eps(kk)=( has(comp_eps{kk}, epj) && ~has(comp_eps{kk} , epj^2));
        end
        idx=find(logic_eps);

        c_current=sym(0);
        for ww=1:length(idx)
            c_current=c_current+comp_eps{idx(ww)}/ep^c(jj);
        end
        ci{ii,jj}=sym(0);
        ci{ii,jj}=simplify(c_current,'Steps',1);
    end
end

for ii=1:5
    for jj=1:4
        ci{ii,jj}=simplify(ci{ii,jj});
    end
end

%% Solution of the resulting algebraic system

% solve for the coefficients

equations=[ci{1,1}==0,ci{1,2}==0,ci{1,3}==0,ci{1,4}==0,ci{2,1}==0,ci{2,2}==0,ci{2,3}==0,ci{2,4}==0,ci{3,1}==0,ci{3,2}==0,ci{3,3}==0,ci{3,4}==0,ci{4,1}==0,ci{4,2}==0,ci{4,3}==0,ci{4,4}==0,ci{5,1}==0,ci{5,2}==0,ci{5,3}==0,ci{5,4}==0];
C=solve(equations, [c02 c04 c22 c24 c33 c35 c44 c55],'IgnoreAnalyticConstraints',true)

toc
%% Saving workspace
tic

timestamp=datetime('now','TimeZone','local','Format','d_MMM_y__HH_mm_ss');
timestamp=string(timestamp);
save(['workspace_harmonic_balance__', num2str(timestamp), '.mat'])

toc