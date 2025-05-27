clear all
close all
clc

%% Set Python environment
% Change with your path

% Windows
env = pyenv('Version', ...
    'your_path\your_environment\Scripts\python.exe');

% Linux/macOS
% env = pyenv('Version', ...
%    '/Users/your_username/your_environment/bin/python');

%% Instantiate a Patient Object
% Import the simulator module
simulator =...
    py.importlib.import_module('python_anesthesia_simulator.simulator');

% Define patient parameters
age = 18;       % years
height = 170;   % cm
weight = 60;    % kg
gender = 0;     % 0 = female, 1 = male
Ts = 1;         % s

% Instantiate a Patient object
George = simulator.Patient([age, height, weight, gender], ts = Ts);

%% Simulation

T_simu = 3600;             % total simulation time in seconds
N_simu = floor(T_simu/Ts); % total number of simulation steps

% Initialize infusion profiles
propofol_infusion_profile = zeros(N_simu,1);        % mg/s
remifentanil_infusion_profile = zeros(N_simu,1);    % ug/s

% Propofol profile
propofol_infusion_profile(1:floor(50/Ts)) = ...
    2;           % 2 mg/s for 50 seconds
propofol_infusion_profile(floor(121/Ts):end) = ...
    0.2;      % 0.2 mg/s from 121s onward

% Remifentanil profile
remifentanil_infusion_profile(1:floor(30/Ts)) = ...
    1;       % 1 ug/s for 30 seconds
remifentanil_infusion_profile(floor(121/Ts):end) = ...
    0.1;  % 0.1 ug/s from 121s onward


df = George.full_sim(u_propo = propofol_infusion_profile,...
    u_remi = remifentanil_infusion_profile);

% Convert Python DataFrame to MATLAB table
George_tbl = table(df);

%% Plot the results
Time = George_tbl.Time;
% Plot control inputs
figure
subplot(2,1,1)
plot(Time, George_tbl.x_propo_4, 'DisplayName', 'Propofol [ug/ml]') 
hold on
plot(Time, George_tbl.x_remi_4, 'DisplayName', 'Remifentanil [ng/ml]')
xlabel('Time [s]')
ylabel('Effect Site Concentration')
legend('Location','southeast')
grid on

subplot(2,1,2)
plot(Time, George_tbl.u_propo, 'DisplayName', 'Propofol [mg/s]')
hold on
plot(Time, George_tbl.u_remi, 'DisplayName', 'Remifentanil [ug/s]')
ylabel('Infusion Rate')
legend('Location','northeast')
grid on


% Plot physiological responses
figure
subplot(4,1,1)
plot(Time, George_tbl.BIS)
ylabel('BIS')
grid on

subplot(4,1,2)
plot(Time, George_tbl.MAP)
ylabel('MAP')
grid on

subplot(4,1,3)
plot(Time, George_tbl.CO)
ylabel('CO')
grid on

subplot(4,1,4)
plot(Time, George_tbl.TOL)
ylabel('TOL')
grid on
xlabel('Time [s]')