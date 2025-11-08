%% simulink of linear and nonlinear model 
%% nonlinear model
k1 = 3.33e-5;    % Pump constant for V1 (m^3/s per Volt)
k2 = 3.35e-5;    % Pump constant for V2 (m^3/s per Volt)
a1 = 0.071e-4;    % Outflow coefficient for h1 (m^2.5/s)
a2 = 0.057e-4;    % Outflow coefficient for h2 (m^2.5/s)
a3 = 0.071e-4;    % Outflow coefficient for h3 (m^2.5/s)
a4 = 0.057e-4;    % Outflow coefficient for h4 (m^2.5/s)
A1 = 28e-4;     % Cross-sectional area of tank 1 (m^2)
A2 = 32e-4;     % Cross-sectional area of tank 2 (m^2)
A3 = 28e-4;     % Cross-sectional area of tank 3 (m^2)
A4 = 32e-4;     % Cross-sectional area of tank 4 (m^2)
g = 9.81;       % Gravity (m/s^2)
landa1 = 0.7;   % Flow distribution ratio for V1
landa2 = 0.6;   % Flow distribution ratio for V2
V1 = 1;         % Input voltage for pump 1 (Volts)
V2 = 1;    % Input voltage for pump 2 (Volts)
kc = 1;

%% linear model 
A = [ -0.0048         0    0.0132         0
         0   -0.0033         0    0.0100;
         0         0   -0.0132         0;
         0         0         0   -0.0100];
B = [ 0.0083         0;
         0    0.0063;
         0    0.0048;
    0.0031         0];

C = [1    0    0    0;
         0    1    0    0;
         0    0    1    0;
         0    0    0    1];

D = [0 0;
     0 0;
     0 0;
     0 0];
%% simulink of controller
