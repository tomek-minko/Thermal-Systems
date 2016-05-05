close all; clear all;
import casadi.*
%% Constants

Pc_polyfit  = 70;
Pe_MT   = 24;
T_sh1 = 10;
nis  = 0.8;

Pc_Hi=50;
Pcnd = 49;
m    = 0.067;
c_pr=1.2;
M_MT=1000;
h_gco=223;
h_eo=450;
UA_MT=0.6;

M_rec =100.;  
c_pw  =4.18;% [kJ/kgK]
M_wt  =1000.;
UA=0.3;
%####### Building Zone #################################################################################################
a     = 100.; 
b     = 20.;
c     = 5.;
m_a   = 1.25; %# air mass 
M_zone=a*b*c*m_a; %#[kg]
C_pair  = 1.000;  %#[kJ/kgK]
UA_ind=a*c*0.000051+a*c*0.000051+2*b*c*0.000051+a*b*0.000072 + a*b*0.00003;

tau=0.1;
UA=0.3; 
M_air_hc=50.; M_air_ca=50.;
    
%% Simulation data
T = 2*60*60;
nmin=15;
Ts=nmin*60; 
N=T/Ts;
Nhours=14;
Hphours=24;
Nh=   Nhours*60/nmin;  % SIMULATION STEPS
Hp = Hphours*60/nmin;                     % Prediction Horizon
Hu = Hphours*60/nmin;                      % Control Horizon

load Data.mat

load N_p_nonoise.mat
y=y_nonoise;
%% Electricity price from year 2015
price=importdata('Market_Data.xlsx','%f %f %f %f %f %f %f %f %f'); %kr/MWh
[k,l]=size(price.textdata(3:end,3));
Day_of_the_ear=10;
pr=[];
for i=1:Nhours+Hphours
    ff=str2num(price.textdata{Day_of_the_ear*24+i,3});
    [a,b]=size(ff);
    if b>1
        h_pr(1,i)=ff(1)*1000+ff(2);
    else
        h_pr(1,i)=ff; %eur/MW
    end
    pr=[pr repmat(h_pr(1,i),1,60)];
end
% pr_est=pr/1000; % eur/kWh
pr_est=pr+4.1;
p_1st=181;
np_1st=8*60;
T_1st=8*60;
% pr_est=pr_est-min(pr_est)+0.1;
%%
%% Disturbances 
Tout_var = 8.; % variation (you can change it according to your preferences)
stopparam=24.; % simulation time (please choose it such that we get full days 24*nr )
tday=24.;      % 24 hours 
nr_days=stopparam/tday ;% number of days
per=[]; 
for i =1:60*tday*nr_days-1
        per=vertcat([per,i+1]) ;        % period with one minute time step (one minut is the sampling time for my plant)
end
        T_outdoor =  Tout_var*sin(2*nr_days/2*pi*per/(60*tday*nr_days)); % Temperature vector
        Tout_avr = 2;                 % average value (you can change it according to your preferences)
        T_outdoor = abs(T_outdoor)+Tout_avr;
        t_amb = T_outdoor(1:N);
        t_ind = 22*ones(1,N);
%% Occupancy vector 
        n_p=3*ones(1,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% v3=[4.2334,3.84279];

%% Declare model variables
%% States
T_mt  = MX.sym('T_mt');
T_co  = MX.sym('T_co');
T_ci  = MX.sym('T_ci');
T_wt  = MX.sym('T_wt');
h_ho  = MX.sym('h_ho');
T_ind  = MX.sym('T_ind');
T_sa  = MX.sym('T_sa');
T_ca  = MX.sym('T_ca');

%% Controls
m_c   = MX.sym('m_c');
m_r   = MX.sym('m_r');
m_gca = MX.sym('m_gca');
P_hi  = MX.sym('P_hi');
OD_1  = MX.sym('OD_1');
UAc   = MX.sym('UAc');
UAh   = MX.sym('UAh');
Q_d   = MX.sym('Q_d');
h_gci  = MX.sym('h_gci');
T_gci  = MX.sym('T_gci');
T_gca = MX.sym('T_gca');
Q_h   = MX.sym('Q_h');
Q_c   = MX.sym('Q_c');
m_sa  = MX.sym('m_sa');
Q_d2  = MX.sym('Q_d2');
OD_2  = MX.sym('OD_2');
%% Disturbances
T_amb = MX.sym('T_amb');
N_p = MX.sym('N_p');
p1 = MX.sym('p1');
p2 = mean(pr_est); %490;%eur/kWh - district heating price; %MX.sym('p2');
% break
%% Model equations
dT_mt = UA_MT*(T_ind-T_mt)/(M_MT*c_pr_LP) + m_r*(h_gco-h_eo)/(M_MT*c_pr_LP);  
dT_wt = UAc*(T_wt-T_ci)/(M_wt*c_pw)-Q_d/(M_wt*c_pw)-Q_d2/(M_wt*c_pw);
dT_co = m_c*(T_ci-T_co)/(M_rec) + UAc*((v1(1)+v1(2)*Pe_MT)*P_hi - T_co)/(M_rec*c_pw);
dT_ci = m_c*(T_co-T_ci)/(M_rec) - UAc*(T_wt-T_ci)/(M_rec*c_pw);
dh_ho = m_r*OD_1*(P_hi*v2(2)+v2(1) - h_ho)/(M_rec) - UAh*( (v1(1)+v1(2)*Pe_MT)*P_hi - T_co)/(M_rec);

dT_ind = m_sa*(T_sa-T_ind)/M_zone + UA_ind*(T_amb - T_ind)/M_zone + tau*N_p/M_zone + UA_MT*(T_mt-T_ind)/M_zone;
dT_sa  = m_sa*(T_ca-T_sa)/M_air_hc + Q_h/M_air_hc + Q_d/M_air_hc;
dT_ca  = m_sa*(OD_2*T_amb + (1-OD_2)*T_ind - T_ca)/M_air_ca - Q_c/M_air_ca; 

xdot = [dT_mt; dT_wt; dT_co; dT_ci; dh_ho; dT_ind; dT_sa; dT_ca];

x=[T_mt;T_co;T_ci;T_wt;h_ho;T_ind;T_sa;T_ca];
d=[T_amb; N_p; p1];
u=[m_c; m_r; m_gca; P_hi; OD_1; UAc; UAh; Q_d; h_gci; T_gca; Q_h; Q_c; m_sa; Q_d2; OD_2;d];

%% Initial Conditions
m_c_0=0.01305;  m_c_min   =0.0001;     m_c_max=3;
m_r_0=0.04806;   m_r_min   =0.0001;     m_r_max=2;
m_gca_0=0.806; m_gca_min =0.0001;     m_gca_max=10;
P_hi_0=50;      P_hi_min  =45;         P_hi_max=90;
OD_1_0=1.;      OD_1_min  =0.0001;      OD_1_max=1;
UAc_0 =0.3;     UAc_min   =0.0001;     UAc_max=5;
UAh_0=0.3;      UAh_min=0.0001;        UAh_max=5;
Q_d_0= 1;       Q_d_min=0.0001;        Q_d_max=10;
Q_h_0= 5;       Q_h_min=0.0001;        Q_h_max=40;
Q_c_0= 0;       Q_c_min=0.001;        Q_c_max=10;
m_sa_0=1.5;     m_sa_min=0.00001;      m_sa_max=10;
Q_d2_0= 1;      Q_d2_min=0.0001;       Q_d2_max=5;
OD_2_0=1.;      OD_2_min  =0.2;        OD_2_max=1;

T_mt_0=3;       T_mt_min=1;         T_mt_max = 5;
T_co_0=60;      T_co_min=( v1(1)+v1(2)*Pe_MT )*P_hi_min;        T_co_max = ( v1(1)+v1(2)*Pe_MT )*P_hi_max;
T_ci_0=40.3747; T_ci_min=10;        T_ci_max = T_co_max;
T_wt_0=55;      T_wt_min=40;        T_wt_max = T_co_max;
T_gca_0=16;  T_gca_min=max(T_outdoor)+5; T_gca_max= 20;
h_ho_0=420;     h_ho_min = 300;   h_ho_max=P_hi_max*v2(2)+v2(1);
h_gci_0=h_ho_0; h_gci_min= h_gco;   h_gci_max=h_ho_max;
T_ind_0=23;     T_ind_min=21;     T_ind_max = 28;
T_sa_0=25;       T_sa_min=15;        T_sa_max = 35;
T_ca_0=16;       T_ca_min=12.2;        T_ca_max = 28;

%% Objective term
L = (p1*1.25*m_r*(P_hi*v2(2)+v2(1)-h_eo) + p1*0.2*m_gca^2 + p1*10*m_c + p1*0.2*m_sa^2 + p1*1.25*Q_c + p2*1.25*Q_h + p1*0.1*Q_d - p2*Q_d2);
L = L+(OD_2*T_amb + (1-OD_2)*T_ind - T_sa)^2;
% Continuous time dynamics
f = Function('f', {x, u}, {xdot, L});

nx = size(x,1);
nu = size(u,1);
nd = size(d,1);

%% Control discretization
M  = 1; % RK4 steps per interval
DT = T/N/M;
X0 = MX.sym('X0', nx);
U  = MX.sym('U',nu);
X  = X0;
Q  = 0;
%% RK 4 Function
for j=1:M
    [k1, k1_q] = easycall(f, X, U);
    [k2, k2_q] = easycall(f, X + DT/2 * k1, U);
    [k3, k3_q] = easycall(f, X + DT/2 * k2, U);
    [k4, k4_q] = easycall(f, X + DT * k3, U);
    X=X+DT/6*(k1 +2*k2 +2*k3 +k4);
    Q = Q + DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q);
end
F = Function('F', {X0, U}, {X, Q});

m_c_nom=m_c_0; UAc_nom=UAc_0;
m_r_nom=m_r_0; UAh_nom=UAh_0;
x0=[T_mt_0; T_co_0; T_ci_0; T_wt_0; h_ho_0; T_ind_0; T_sa_0; T_ca_0];
u0=[ m_c_0; m_r_0; m_gca_0; P_hi_0; OD_1_0; UAc_0; UAh_0; Q_d_0;  h_gci_0; T_gca_0; Q_h_0; Q_c_0; m_sa_0; Q_d2_0; OD_2_0];


for i=1:Nh
tic
    %% 
% % Start with an empty NLP
w={};
w0 = [];
lbw = [];
ubw = [];
g={};
lbg = [];
ubg = [];

%% "Lift" initial conditions
X0 = MX.sym('X0', nx);
w = {w{:}, X0};


lbw = [lbw; x0]; 
ubw = [ubw; x0]; 
w0 =  [w0;  x0];
% Simple bounds on states
lbx = {};
ubx = {};


% Formulate the NLP
Xk = X0;
J=0;

for k=0:N-1
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)],nu-nd);
    w = {w{:}, Uk};
    lbw = [lbw; m_c_min; m_r_min; m_gca_min; P_hi_min; OD_1_min; UAc_min; UAh_min; Q_d_min;  h_gci_min; T_gca_min;  Q_h_min; Q_c_min; m_sa_min; Q_d2_min; OD_2_min]; %u_min=[m_c; m_r; m_gca; P_hi; OD_1; UAc; UAh; Q_d; d ];
    ubw = [ubw; m_c_max; m_r_max; m_gca_max; P_hi_max; OD_1_max; UAc_max; UAh_max; Q_d_max;  h_gci_max; T_gca_max;  Q_h_max; Q_c_max; m_sa_max; Q_d2_max; OD_2_max];
    w0 = [w0;  u0];

    % Integrate till the end of the interval
    [Xk_end, Jk] = easycall(F, Xk, [Uk;T_outdoor((i+k)*nmin+T_1st); y((i+k)*nmin+np_1st);pr_est((i+k)*nmin+k+p_1st)]);
    J=J+Jk;
    
    % New NLP variable for state at end of interval
    Xk = MX.sym(['X_' num2str(k+1)], nx);
    w = {w{:}, Xk};
    lbw = [lbw; T_mt_min; T_co_min; T_ci_min; T_wt_min; h_ho_min; T_ind_min;T_sa_min;T_ca_min];
    ubw = [ubw; T_mt_max; T_co_max; T_ci_max; T_wt_max; h_ho_max; T_ind_max;T_sa_max;T_ca_max];
    w0 =  [w0;  x0];
        
    % Add equality constraint
    g = {g{:}, Xk_end-Xk};
    lbg = [lbg; 0; 0; 0; 0; 0; 0; 0; 0];
    ubg = [ubg; 0; 0; 0; 0; 0; 0; 0; 0];
    
    g={g{:}, Xk(7)-Xk(8)};
    lbg = [lbg; 0];
    ubg = [ubg; inf];
    
    g={g{:}, Xk(6)-Xk(8)};
    lbg = [lbg; 0];
    ubg = [ubg; inf];
    
    g={g{:}, Uk(9)-(Uk(5)*Xk(5) + (1-Uk(5))*(Uk(4)*v2(2)+v2(1)))};
    lbg = [lbg; 0];
    ubg = [ubg; 0];
%%
    g={g{:}, Uk(6)-UAc_nom*Uk(1)/m_c_nom };
    lbg = [lbg; 0];
    ubg = [ubg; 0];
    g={g{:}, Uk(7)-UAh_nom*Uk(5)*Uk(2)/m_r_nom };
    lbg = [lbg; 0];
    ubg = [ubg; 0];
 %% Energy balance equation across the Heat Exchanger  
    g={g{:},Uk(5)*Uk(2)*(Uk(4)*v2(2)+v2(1)-Xk(5)) - Uk(1)*c_pw*(Xk(2)-Xk(3)) };
    lbg = [lbg; 0];
    ubg = [ubg; 0];
%% Energy balance equation across the Gas Cooler    
    g={g{:},Uk(2)*(Uk(9)-h_gco) - Uk(3)*(Uk(10)-T_outdoor((i+k)*nmin+T_1st)) };
    lbg = [lbg; 0];
    ubg = [ubg; 0];
end


opts=struct;
opts.ipopt.tol=1e-1;
opts.ipopt.fixed_variable_treatment='make_constraint';%'relax_bounds';%
opts.ipopt.bound_push=0.1;
opts.ipopt.dual_inf_tol=0.8;
opts.ipopt.recalc_y='yes';
opts.ipopt.recalc_y_feas_tol=1e-2;
opts.ipopt.compl_inf_tol=0.1;
opts.ipopt.constr_viol_tol=0.01;

prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob,opts);

%% Solve the NLP

arg = struct('x0', w0, 'lbx', lbw, 'ubx', ubw,...
            'lbg', lbg, 'ubg', ubg);

i
sol = solver(arg);

w_opt = full(sol.x);

u_1st=full(sol.x(9:23));




[x_next, J_next] = easycall(F, x0, [u_1st;T_outdoor((i)*nmin+T_1st); y((i)*nmin+np_1st);pr_est(i*nmin+p_1st)]);
x0=x_next;
u0=u_1st;
Xi(:,i)=full(x0);
Ui(:,i)=full(u_1st);
t(i)=toc;
end
T_amb  = T_outdoor(T_1st:nmin:T_1st+Nh*nmin-1);
N_p    = y(T_1st:nmin:T_1st+Nh*nmin-1);
p1 = pr_est(p_1st:nmin:p_1st+Nh*nmin-1);

Dist=[T_amb;N_p;p1];
save('time_2h_warm_whole_sys','t');
save('States_2h_warm_whole_sys','Xi')
save('Ctr_2h_warm_whole_sys','Ui')
save('Solution_2h_warm_whole_sys','sol')
save('Disturbance.mat','Dist')

T_mt =  Xi(1,:);
T_co =  Xi(2,:);
T_ci =  Xi(3,:);
T_wt =  Xi(4,:);
h_ho =  Xi(5,:);
T_ind  =  Xi(6,:);
T_sa  =Xi(7,:);
T_ca  =Xi(8,:);

m_c = Ui(1,:);
m_r = Ui(2,:);
m_gca = Ui(3,:);
P_hi = Ui(4,:);
OD_1 = Ui(5,:);
UA_c = Ui(6,:);
UA_h = Ui(7,:);
Q_d = Ui(8,:);
h_gci= Ui(9,:);
T_gca= Ui(10,:);
Q_h= Ui(11,:);
Q_c= Ui(12,:);
m_sa= Ui(13,:);
Q_d2= Ui(14,:);
OD_2= Ui(15,:);



figure; grid on; hold on;
plot(T_co);plot(T_ci);plot(T_wt);plot(T_gca);
legend('T_{co}','T_{ci}','T_{wt}','T_{gca}');

figure; grid on; hold on;
plot(T_ind);plot(T_sa);plot(T_ca);plot(T_mt);
legend('T_{ind}','T_{sa}','T_{ca}','T_{mt}');
%% Controls 
figure; grid on; hold on;
plot(m_c);plot(m_r);plot(m_gca);plot(m_sa);
legend('m_c','m_r','m_{gca}','m_{sa}');

figure; grid on; hold on;
plot(OD_1);plot(UA_c);plot(UA_h);plot(OD_2);
legend('OD_1','UA_c','UA_h','OD_2');

figure; grid on; hold on;
plot(P_hi); plot(Q_d); plot(Q_c);plot(Q_h);legend('P_{hi}','Q_d','Q_c','Q_h')

figure; grid on; hold on;
plot(h_gci); plot(h_ho);  legend('h_{gci}','h_{ho}');

figure; hold on; 
plot(T_amb); plot(N_p); plot(p1); 
legend('T_{ambient}','Occupancy','Electricity price')
grid on;

figure; plot(t,'or'); legend('Simulation time');grid on;
%% Power consumption
L = p1.*1.25.*m_r.*(P_hi.*v2(2)+v2(1)-h_eo) + p1.*0.2.*m_gca.^2 + p1.*10.*m_c + p1.*0.2.*m_sa.^2 + p1.*1.25.*Q_c + p2.*1.25.*Q_h + p1.*0.1.*Q_d - p2.*Q_d2;
W_comp = 1.25.*m_r.*(P_hi.*v2(2)+v2(1)-h_eo);
W_gcf  = 0.2.*m_gca.^2;
W_saf  = 0.2.*m_sa.^2;
W_cp = 10.*m_c;
W_dp = 0.1.*Q_d;
W_c  = 1.25.*Q_c;
W_h  = 1.25.*Q_h;

figure; grid on; hold on;
plot(L); legend('Objective function')
figure; grid on; hold on;
plot(W_comp); plot(W_gcf); plot(W_saf); plot(W_cp); plot(W_dp); plot(W_c); plot(W_h)
legend('W_comp','W_gcf','W_saf','W_cp','W_dp','W_c','W_h')




