clear all
close all
rng(2025)

%% select data
% data_name = 'artificial';
% data_name = 'cavity';
data_name = 'cylinder';
% data_name = 'dam';    % recommended
% data_name = 'tube';   % recommended
% data_name = 'neuron';   % not recommended
% data_name = 0;
tr = false;

% switch data_name
%     case 'artificial'
%         n = 600;    % dimension
%         steps = 201;    % data size
%         %% generate matrix and data
%         % [A_org, evals, evecs] = rand_mat_real(n);
%         % [A_org, evals, evecs] = rand_mat_sym(n);
%         % [A_org, evals, evecs] = case1(n);
%         [A_org, evals, evecs] = case2(n);
%         distri = ones(n, 1);
%         x0 = evecs * distri;    % control the beginning components
%         data = zeros([n, steps]);
%         data(:, 1) = x0;
%         for i = 2:steps
%             data(:, i) = A_org * data(:, i-1);
%         end
%         % data = data + 1e-8 * eye(size(data));
%         % data = data + 1e-6 * randn(size(data));     % Robust!
% 
%         nx = n; ny = 1;
% 
%     case 'cavity'
%         reyn = 13;     % 13, or 16, 19, 20, 30
%         load_time = 200;
%         file_name = ['Cavity', num2str(reyn), 'k.mat'];
%         load(file_name);
%         nx = VelocityField.N+1; ny = nx;    tr = true;
%         data = VelocityField.Psi(:, 1: load_time);
% 
%     case 'cylinder'
%         file_name = 'CYLINDER_ALL.mat';
%         load(file_name);
%         data = VORTALL;
% 
%     case 'dam'
%         file_name = 'dam_bc_case0000.mat';
%         load(file_name);
%         data = vort;
%         nx = 64; ny = 64;   tr = true;
% 
%     case 'tube'
%         file_name = 'tube_bc_case0000.mat';
%         load(file_name);
%         data = vort;
%         nx = 64; ny = 64;   tr = true;
% 
%     case 'neuron'
%         file_name = 'ecog_window.mat';
%         load(file_name);
%         data = X(:, 1:200);
%         nx = 59; ny = 1;
% 
%     otherwise
%         disp('================================================================')
%         disp('Failed to find data!')
%         disp('================================================================')
%         return
% end 

%% hyper parameters
r = 10;     % lower-rank

% steps = size(data, 2);  % get the time steps
% init = round(steps / 2); w = init;  % recommended setting

% view_data(data, nx, ny, tr)

% %% solve
% [evals_on, vr_on, P_on, B_on] = online_iteration_windowed(data, init, r, w);
% evals = evals_on; vr = vr_on; P = P_on; x0 = data(:, 1); xn = data(:, end);
% save([data_name, '.mat'], 'evals', 'vr', 'P', 'init', 'steps', 'nx', 'ny', 'tr', 'x0', 'xn')
% 
% [evals_dmd, vr_dmd, P_dmd] = baseline_simple(data, steps, r, 'dmd');
% evals = evals_dmd; vr = vr_dmd; P = P_dmd;
% save([data_name, '_dmd.mat'], 'evals', 'vr', 'P', 'init', 'steps', 'nx', 'ny', 'tr', 'x0', 'xn')


view_modes_by_rank(data_name)

%% do mode visualization (by iter)
steps = 100;
view_modes4(data_name, 0.1)
% % view_mode(data_name, 2, 0.1)
% view_mode(data_name, 3, 0.1)

view_pred([data_name, '_dmd.mat'], steps)   % snapshots recovered by dmd
view_pred([data_name, '.mat'], steps)   % snapshots recovered by our alg
% view_data(data, nx, ny, tr)




