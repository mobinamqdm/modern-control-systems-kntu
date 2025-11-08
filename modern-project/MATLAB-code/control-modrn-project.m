% clc;
% clear all;
% close all;

%% 1) system variables (state and input) and parameters
syms h1 h2 h3 h4 V1 V2
syms k1 k2 a1 a2 a3 a4 g A1 A2 A3 A4 landa1 landa2 kc

%% 2) nonlinear system equations (differential equations)
dh1_dt = (landa1*k1 * V1) / A1 + a3 * sqrt(2 * g * h3) / A1 - a1 * sqrt(2 * g * h1) / A1;
dh2_dt = (landa2*k2 * V2) / A2 + a4 * sqrt(2 * g * h4) / A2 - a2 * sqrt(2 * g * h2) / A2;
dh3_dt = (1 - landa2) * k2 * V2 / A3 - a3 * sqrt(2 * g * h3) / A3;
dh4_dt = (1 - landa1) * k1 * V1 / A4 - a4 * sqrt(2 * g * h4) / A4;

nonlinear_equations = {dh1_dt, dh2_dt, dh3_dt, dh4_dt};
disp('Nonlinear Equations:');
for i = 1:length(nonlinear_equations)
    disp(nonlinear_equations{i});
end

%% 3) output equations (here we consider h1 and h2 only but h3 and h4 can be observed too as the paper has analays them)
y1 = kc * h1;
y2 = kc * h2;

%% 4) state and input vectors 
x = [h1; h2; h3; h4];
u = [V1; V2];        

%% 5) symbolic Jacobian matrices
A_sym = jacobian([dh1_dt; dh2_dt; dh3_dt; dh4_dt], x);
B_sym = jacobian([dh1_dt; dh2_dt; dh3_dt; dh4_dt], u);
C_sym = jacobian([y1; y2], x);
D_sym = jacobian([y1; y2], u);

disp('Symbolic A Matrix:');
disp(A_sym);
disp('Symbolic B Matrix:');
disp(B_sym);
disp('Symbolic C Matrix:');
disp(C_sym);
disp('Symbolic D Matrix:');
disp(D_sym);

%% 6) numerical values for the parameters
k1_val    = 3.33e-5;
k2_val    = 3.35e-5;
a1_val    = 0.071e-4;
a2_val    = 0.057e-4;
a3_val    = 0.071e-4;
a4_val    = 0.057e-4;
A1_val    = 28e-4;
A2_val    = 32e-4;
A3_val    = 28e-4;
A4_val    = 32e-4;
g_val     = 9.81;
landa1_val= 0.7;
landa2_val= 0.6;
kc_val    = 1;
%% 7) Substitute numerical parameter values into the differential equations
dh1_dt = subs(dh1_dt, {k1, k2, a1, a2, a3, a4, A1, A2, A3, A4, g, landa1, landa2, kc}, ...
                     {k1_val, k2_val, a1_val, a2_val, a3_val, a4_val, A1_val, A2_val, A3_val, A4_val, g_val, landa1_val, landa2_val, kc_val});
dh2_dt = subs(dh2_dt, {k1, k2, a1, a2, a3, a4, A1, A2, A3, A4, g, landa1, landa2, kc}, ...
                     {k1_val, k2_val, a1_val, a2_val, a3_val, a4_val, A1_val, A2_val, A3_val, A4_val, g_val, landa1_val, landa2_val, kc_val});
dh3_dt = subs(dh3_dt, {k1, k2, a1, a2, a3, a4, A1, A2, A3, A4, g, landa1, landa2, kc}, ...
                     {k1_val, k2_val, a1_val, a2_val, a3_val, a4_val, A1_val, A2_val, A3_val, A4_val, g_val, landa1_val, landa2_val, kc_val});
dh4_dt = subs(dh4_dt, {k1, k2, a1, a2, a3, a4, A1, A2, A3, A4, g, landa1, landa2, kc}, ...
                     {k1_val, k2_val, a1_val, a2_val, a3_val, a4_val, A1_val, A2_val, A3_val, A4_val, g_val, landa1_val, landa2_val, kc_val});

%% 8) Equilibrium Calculation
% i first set the input v1 and v2 to zero but i faced division by zero
% because all the state and input variables became zero
% according to my research the only vay was to fiund another input other
% than zero 

V1_input = 1;
V2_input = 1;

eqns = [ subs(dh1_dt, {V1,V2}, {V1_input,V2_input}) == 0, ...
         subs(dh2_dt, {V1,V2}, {V1_input,V2_input}) == 0, ...
         subs(dh3_dt, {V1,V2}, {V1_input,V2_input}) == 0, ...
         subs(dh4_dt, {V1,V2}, {V1_input,V2_input}) == 0];

% (For square-root expressions we have these conditions h1,h2,h3,h4 > 0)
assume(h1 > 0);
assume(h2 > 0);
assume(h3 > 0);
assume(h4 > 0);

equilibrium_points = solve(eqns, [h1, h2, h3, h4], 'ReturnConditions', true); 
disp('Validity Conditions for the equilibrium solution:');
disp(equilibrium_points.conditions);

% Convert the computed equilibrium points into numerical values
eq_h1 = double(vpa(equilibrium_points.h1, 10));
eq_h2 = double(vpa(equilibrium_points.h2, 10));
eq_h3 = double(vpa(equilibrium_points.h3, 10));
eq_h4 = double(vpa(equilibrium_points.h4, 10));

disp('Equilibrium Points (numerical):');
disp(['h1 = ' num2str(eq_h1)]);
disp(['h2 = ' num2str(eq_h2)]);
disp(['h3 = ' num2str(eq_h3)]);
disp(['h4 = ' num2str(eq_h4)]);

x_eq = [eq_h1; eq_h2; eq_h3; eq_h4];
u_eq = [V1_input; V2_input];
%% 9) Substitute numeric equilibrium and parameter values into A, B, C, D
subsList = { h1,  h2,  h3,  h4,  k1,    V1,       V2,      k2,   a1,    a2,    a3,    a4,    g,     A1,    A2,    A3,    A4,   landa1,   landa2, kc };
valuesList = { eq_h1, eq_h2, eq_h3, eq_h4, k1_val, V1_input, V2_input, k2_val, a1_val, a2_val, a3_val, a4_val, g_val, A1_val, A2_val, A3_val, A4_val, landa1_val, landa2_val, kc_val };

A_numeric = double(subs(A_sym, subsList, valuesList));
B_numeric = double(subs(B_sym, subsList, valuesList));
C_numeric = double(subs(C_sym, subsList, valuesList));
D_numeric = double(subs(D_sym, subsList, valuesList));

disp('A Matrix at the equilibrium point:');
disp(A_numeric);
disp('B Matrix at the equilibrium point:');
disp(B_numeric);
disp('C Matrix at the equilibrium point:');
disp(C_numeric);
disp('D Matrix at the equilibrium point:');
disp(D_numeric);

%% The stability Analysis
eig_values = eig(A_numeric);
disp('Eigenvalues of A Matrix:');
disp(eig_values);
if all(real(eig_values) < 0)
    disp('The system is asymptotically stable.');
else
    disp('The system is unstable or marginally stable.');
end
%% 10) Controllability, Observability Analysis with controllability matrix
Controllability = ctrb(A_numeric, B_numeric);
Observability   = obsv(A_numeric, C_numeric);
disp('Controllability Matrix:');
disp(Controllability);
disp('Rank of Controllability Matrix:');
disp(rank(Controllability));
rank_C = rank(Controllability);
if rank_C == length(x) 
    disp('The system is controllable with controllability matrix.'); 
else 
    disp('The system is not controllable with controllability matrix.'); 
end 
disp('Observability Matrix:');
disp(Observability);
disp('Rank of Observability Matrix:');
disp(rank(Observability));
rank_O = rank(Observability);
if rank_O == length(x) 
    disp('The system is observable with observability matrix.'); 
else 
    disp('The system is not observable with observability matrix.'); 
end

%%  Controllability, Observability Analysis with PBH Test
fprintf("PBH TEST:\n")
for i=1:4
    phi_o=[eig_values(i,1)*eye(4)-A_numeric;C_numeric];
    r(i)=rank(phi_o);
    fprintf("for the eigen value of the below,the rank of the matrix is %d\n",r(i));
    disp(eig_values(i))
end
if r==[4 ,4 ,4 ,4]
    fprintf("This system is observable with PBH test");
end
for i=1:4
    phi_c=[eig_values(i,1)*eye(4)-A_numeric,B_numeric];
    r(i)=rank(phi_c);
    fprintf("\nfor the eigen value of the below,the rank of the matrix is %d\n",r(i));
    disp(eig_values(i))
end
if r==[4 ,4 ,4 ,4]
    fprintf("This system is controlable with PBH test\n");
end
%% controlability and observability with jordan form Matrix
[v,j]=jordan(A_numeric);
c_new=C_numeric*v;
B_new=inv(v)*B_numeric;
fprintf("Jordan form of the A Matrix:\n");
disp(j);
fprintf("The new matrix of B that is formed by Similarity Transformation:\n");
disp(B_new);
fprintf("The new matrix of C that is formed by Similarity Transformation:\n");
disp(c_new);
%% 11) Linear Simulation using lsim
sys = ss(A_numeric, B_numeric, C_numeric, D_numeric);
t_sim = 0:1:10000;  
%initial_conditions = [1.5 0.8 2 1]; 
%initial_conditions = [0 0 0 0]; 
initial_conditions = [1 2 3 4]; 
u_sim = ones(length(t_sim), 2); 
[y_linear, t_linear] = lsim(sys, u_sim, t_sim,initial_conditions);
figure('Color', 'w','Position', [100, 100, 800, 600]);
subplot(2,1,1);
yyaxis left
plot(t_linear,y_linear(:,1), 'r-', 'LineWidth', 2);
hold on;
plot(t_linear, y_linear(:,2), 'b-', 'LineWidth', 2);
ylabel('Tank Levels (h1 & h2)', 'FontSize', 14);
yyaxis right
pump_input = @(t) 1.0*(t>=0);
plot(t_linear, pump_input(t_linear), 'k--', 'LineWidth', 2);
ylabel('Pump Input', 'FontSize', 14);
xlabel('Time (s)', 'FontSize', 14);
title('Linear System Step Response with Pump Input', 'FontSize', 16);
legend('h1 (Tank 1)', 'h2 (Tank 2)', 'Pump Input', 'Location', 'best');
grid on;
set(gca, 'FontSize', 12);

%% 12) Nonlinear Simulation using ode45
pump_input = @(t) 1.0*(t>=0);
nonlinear_ode = @(t, x) [ (landa1_val*k1_val * pump_input(t)) / A1_val + a3_val * sqrt(2*g_val*x(3)) / A1_val - a1_val * sqrt(2*g_val*x(1)) / A1_val;
                         (landa2_val*k2_val * pump_input(t)) / A2_val + a4_val * sqrt(2*g_val*x(4)) / A2_val - a2_val * sqrt(2*g_val*x(2)) / A2_val;
                         (1 - landa2_val)* k2_val * pump_input(t) / A3_val - a3_val * sqrt(2*g_val*x(3)) / A3_val;
                         (1 - landa1_val)* k1_val * pump_input(t) / A4_val - a4_val * sqrt(2*g_val*x(4)) / A4_val ];
                     

[t_nonlinear, y_nonlinear] = ode45(nonlinear_ode, t_sim, initial_conditions);

subplot(2,1,2);
yyaxis left
plot(t_nonlinear, y_nonlinear(:,1), 'r-', 'LineWidth', 2);
hold on;
plot(t_nonlinear, y_nonlinear(:,2), 'b-', 'LineWidth', 2);
hold on;
ylabel('Tank Levels (h1 & h2)', 'FontSize', 14);
yyaxis right
plot(t_nonlinear, pump_input(t_nonlinear), 'k--', 'LineWidth', 2);
ylabel('Pump Input', 'FontSize', 14);
xlabel('Time (s)', 'FontSize', 14);
title('Nonlinear System Step Response with Pump Input', 'FontSize', 16);
legend('h1 (Tank 1)', 'h2 (Tank 2)','pump input', 'Location', 'best');
grid on;
set(gca, 'FontSize', 12);
%% 13 step response,impulse response and ramp response for nonlinearity system
V1_step  = @(t)  1.0*(t>=0);
V2_step  = @(t) 1.0*(t>=0);
V1_ramp  = @(t) t;  
V2_ramp  = @(t) t;
impulse_duration = 0.1;    
impulse_amplitude = 1/impulse_duration;  
V1_impulse = @(t)  (t < impulse_duration)*impulse_amplitude;
V2_impulse = @(t)  (t < impulse_duration)*impulse_amplitude;
make_ode = @(V1_fun, V2_fun) @(t, x) [
    (landa1_val*k1_val * V1_fun(t)) / A1_val + a3_val * sqrt(2*g_val*x(3)) / A1_val - a1_val * sqrt(2*g_val*x(1)) / A1_val;
    (landa2_val*k2_val * V2_fun(t)) / A2_val + a4_val * sqrt(2*g_val*x(4)) / A2_val - a2_val * sqrt(2*g_val*x(2)) / A2_val;
    (1 - landa2_val)* k2_val * V2_fun(t) / A3_val - a3_val * sqrt(2*g_val*x(3)) / A3_val;
    (1 - landa1_val)* k1_val * V1_fun(t) / A4_val - a4_val * sqrt(2*g_val*x(4)) / A4_val
];
input_types = {'Step','Ramp','Impulse'};
V1_inputs = {V1_step, V1_ramp, V1_impulse};
V2_inputs = {V2_step, V2_ramp, V2_impulse};
figure;
for k = 1:3
    V1_fun = V1_inputs{k};
    V2_fun = V2_inputs{k};
    ode_fun = make_ode(V1_fun, V2_fun);
    [t_nl, y_nl] = ode45(ode_fun, t_sim, [0;0;0;0]);
    V1_values = arrayfun(V1_fun, t_nl);
    V2_values = arrayfun(V2_fun, t_nl);
    subplot(3,2,(k-1)*2+1);
    plot(t_nl, V1_values, 'm', 'LineWidth', 1.5); hold on;
    title([input_types{k}, ' Input Signals']);
    xlabel('Time (s)');
    ylabel('Input (V1, V2)');
    legend('Inputs');
    grid on;
    subplot(3,2,(k-1)*2+2);
    plot(t_nl, y_nl(:,1), 'r', 'LineWidth', 2); hold on;
    plot(t_nl, y_nl(:,2), 'b', 'LineWidth', 2);
    title(['Nonlinear Response to ', input_types{k}]);
    xlabel('Time (s)');
    ylabel('Tank Levels (m)');
    legend('h1 (Tank 1)','h2 (Tank 2)');
    grid on;
end

%% design controller
desired_poles = [-0.6 -0.3 -0.4 -0.2];  
%desired_poles = [-4 -3 -5 -10];
K = place(A_numeric,B_numeric,desired_poles);
disp('state feedback K:');
disp(K);
Acl = A_numeric - B_numeric*K;
disp('closed loop poles:');
disp(eig(Acl));
k1=inv(-C_numeric*inv(Acl)*B_numeric);
 
A_aug = [A_numeric, zeros(4,2); 
         -C_numeric, zeros(2,2)];
B_aug = [B_numeric; zeros(2,2)];
C_aug = [C_numeric, zeros(2,2)];
D_aug = D_numeric;
desired_poles1=[-0.008 -0.0008 -0.0005 -0.0008 -0.0005 -0.007];  

k2 = place(A_aug,B_aug,desired_poles1);
m = A_aug - B_aug*k2;
N=[114.5 -128 139 -209;251.35 -331.226 235.86 -329.4773];
L=[-7.3155 -2.9044 ;-23.16 20.38];
K_i=k2(:,5:6);
k_f=k2(:,1:4);
%% design observar
poles_desire=[-6 -3 -4 -2];
L=place(A_numeric',C_numeric',poles_desire);
l_t=L';