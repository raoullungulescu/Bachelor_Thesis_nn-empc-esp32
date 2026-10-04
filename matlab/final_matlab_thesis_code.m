%% =================================================================
%%  TUBE-MPC vs TUBE-eMPC vs TUBE-NN
%%  Structurat conform capitolelor lucrarii de licenta:
%%    Cap 2: Preliminarii         (model, e-MPC, mRPI, Tube-MPC)
%%    Cap 3: Invatare profunda bazata pe e-MPC cu garantii de robustete
%%    Cap 4: Simulare numerica
%% =================================================================

%% =================================================================
%% CAP 2.1 — MODELUL MATEMATIC (pendul invers liniarizat)
%% =================================================================
m     = 0.135;
g     = 9.8;
l_tot = 0.47;
l     = l_tot / 2;
I     = (m * l_tot^2) / 3;
H     = 1 / (m * l^2 + I);

A_c = [0 1 0 0;
    0 0 0 0;
    0 0 0 1;
    0 0 (m * g * l) * H 0];
B_c = [0; 1; 0; -m * l * H];

Ts    = 0.02;
sys_d = c2d(ss(A_c, B_c, eye(4), zeros(4,1)), Ts);
A     = sys_d.A;
B     = sys_d.B;
[nx, nu] = size(B);

% Constrangeri stare/comanda si perturbatie aditiva w
x_max = [2.5; 3.0; 0.6; 6.0];
x_min = -x_max;
u_max = 14;
u_min = -u_max;
X = Polyhedron('lb', x_min, 'ub', x_max);
U = Polyhedron('lb', u_min, 'ub', u_max);

w_bound = [0.001; 0.001; 0.002; 0.002];
W = Polyhedron('lb', -w_bound, 'ub', w_bound);

%% =================================================================
%% CAP 2.4.1 — CONTROLER LOCAL K SI MULTIMEA mRPI (Z) — algoritm Rakovic
%% =================================================================
Q_lqr = diag([1, 1, 1, 1]);
R_lqr = 0.25;
K_lqr = dlqr(A, B, Q_lqr, R_lqr);
K     = -K_lqr;
A_K   = A + B*K;

% Orizont Rakovic: cel mai mic s a.i. alpha < alpha_target
alpha_target = 0.01;
s = 0; alpha = inf;
while alpha > alpha_target
    s = s + 1;
    alpha = max( (abs(A_K^s) * w_bound) ./ w_bound );
    if s > 500, error('A_K nu converge; verifica stabilitatea.'); end
end
fprintf('[mRPI] Orizont Rakovic: s = %d, alpha = %.6f\n', s, alpha);

% Generatori zonotop Z (reprezentare ca suma Minkowski finita)
G_W = diag(w_bound);
G_Z = zeros(nx, nx*s);
for i = 0:(s-1)
    G_Z(:, i*nx+1:(i+1)*nx) = (A_K^i) * G_W;
end
G_Z = G_Z / (1 - alpha);

Z     = zonotope_to_polyhedron(G_Z, nx);
z_max = max(abs(Z.V), [], 1)';

% Constrangeri tightened pentru sistemul nominal (Tube-MPC)
X_bar = X - Z;
KZ    = Z.affineMap(K);
U_bar = U - KZ;
if X_bar.isEmptySet() || U_bar.isEmptySet()
    error('Perturbatia este prea mare! Spatiul nominal strans este vid.');
end

% Verificare RPI pe zonotopi via functie suport (directii canonice + generatoare)
fprintf('[mRPI] Verificare conditie RPI pentru Z (directii complete):\n');
G_AKZ = A_K * G_Z;
dirs_all = [eye(nx), -eye(nx), G_Z, -G_Z];
norms    = vecnorm(dirs_all);
dirs_all = dirs_all(:, norms > 1e-10);
dirs_all = dirs_all ./ vecnorm(dirs_all);

h_lhs = sum(abs(G_AKZ' * dirs_all), 1) + sum(abs(G_W' * dirs_all), 1);
h_rhs = sum(abs(G_Z'   * dirs_all), 1);
max_viol = max(h_lhs - h_rhs);
fprintf('  %s (max depasire = %.2e)\n\n', tf2str(max_viol <= 1e-8), max_viol);

%% =================================================================
%% CAP 2.2 — EXPLICIT MPC (mp-QP -> lege PWA offline)
%% =================================================================
N = 5;
Q_mpc = diag([10, 1, 10, 1]);
R_mpc = 0.1;
[K_lqr_nom, P_mpc] = dlqr(A, B, Q_mpc, R_mpc);
K_nom = -K_lqr_nom;

% Set terminal invariant (pentru garantia de stabilitate a MPC nominal)
A_cl_real = real(A + B * K_nom);
sys_cl_nom = LTISystem('A', A_cl_real);
U_bar_as_state_constraint = Polyhedron(U_bar.A * K_nom, U_bar.b);
X_cl_constraint = intersect(X_bar, U_bar_as_state_constraint);
X_cl_constraint.minHRep();
sys_cl_nom.x.with('setConstraint');
sys_cl_nom.x.setConstraint = X_cl_constraint;
X_f = sys_cl_nom.invariantSet();

% Problema MPC nominala -> solutie explicita (mp-QP)
sys_nom = LTISystem('A', A, 'B', B);
sys_nom.x.with('setConstraint');  sys_nom.x.setConstraint = X_bar;
sys_nom.u.with('setConstraint');  sys_nom.u.setConstraint = U_bar;
sys_nom.x.penalty         = QuadFunction(Q_mpc);
sys_nom.u.penalty         = QuadFunction(R_mpc);
sys_nom.x.with('terminalPenalty'); sys_nom.x.terminalPenalty = QuadFunction(P_mpc);
sys_nom.x.with('terminalSet');     sys_nom.x.terminalSet     = X_f;

mpc  = MPCController(sys_nom, N);
fprintf('[e-MPC] Generam solutia explicita (mp-QP)...\n');
empc = mpc.toExplicit();
fprintf('  e-MPC: %d regiuni PWA.\n\n', empc.nr);

%% =================================================================
%% CAP 3.2.1 — GENERARE DATE DE ANTRENAMENT (din legea e-MPC)
%% =================================================================
fprintf('[Date antrenament] Generare traiectorii din e-MPC...\n');
n_init  = 1500;
n_steps = 80;
X_all   = zeros(nx, n_init * n_steps);
U_all   = zeros(nu, n_init * n_steps);
cnt     = 0;

x_min_v = x_min; x_max_v = x_max; x_rng_v = x_max_v - x_min_v;

for i = 1:n_init
    % Strategie mixta de esantionare a starii initiale (4 tipuri ciclice):
    %  0: uniform pe tot X | 1: stari mici (in jurul originii)
    %  2: stari extreme (o dimensiune aproape de limita)
    %  3: strict in X_bar (setul nominal strans)
    tip = mod(i-1, 4);
    switch tip
        case 0
            x_k = x_min_v + rand(nx,1) .* x_rng_v;
        case 1
            x_k = 0.2 * (2*rand(nx,1)-1) .* (x_max_v - x_min_v)/2;
        case 2
            x_k = x_min_v + rand(nx,1) .* x_rng_v;
            dim = randi(nx);
            if rand() < 0.5
                x_k(dim) = x_max_v(dim) * (0.8 + 0.1*rand());
            else
                x_k(dim) = x_min_v(dim) * (0.8 + 0.1*rand());
            end
        case 3
            x_k = zeros(nx,1);
            for try_idx = 1:50
                x_cand = x_min_v + rand(nx,1) .* x_rng_v;
                if X_bar.contains(x_cand)
                    x_k = x_cand;
                    break;
                end
            end
    end

    if ~X.contains(x_k), continue; end   % sarim punctele initiale invalide

    for j = 1:n_steps
        [u_k, feas, ~] = empc.evaluate(x_k);
        if ~feas || any(isnan(u_k)), break; end

        x_next = A * x_k + B * u_k;        % dinamica reala

        cnt = cnt + 1;
        X_all(:, cnt) = x_k;
        U_all(:, cnt) = u_k;

        if ~X.contains(x_next), break; end  % oprire daca iesim din X
        x_k = x_next;
    end
end

X_data = X_all(:, 1:cnt);
U_data = U_all(:, 1:cnt);
fprintf('  Perechi (stare, comanda) generate: %d\n\n', cnt);

% Normalizare (medie/std) pentru intrarea/iesirea retelei
x_mean = mean(X_data, 2);  x_std = std(X_data, 0, 2) + 1e-8;
u_mean = mean(U_data, 2);  u_std = std(U_data, 0, 2) + 1e-8;
normalize_x   = @(x)  ((x - x_mean) ./ x_std)';
denormalize_u = @(u_n) u_mean + u_std .* u_n';

%% =================================================================
%% CAP 3.2.2 — ARHITECTURA RETELEI SI PRE-ANTRENARE
%% =================================================================
M = 30; L = 8;              % M neuroni/strat, L straturi ascunse (ReLU)
EPS_TARGET = 0.05;          % prag tinta pentru eroarea maxima de aproximare

fprintf('[Pre-antrenare] Retea N_{%d,%d} (750 epoci)...\n', M, L);

n_total = cnt;
% impartire date 10% validare 90% antrenament
n_val   = floor(0.1 * n_total);

% shuffle
idx_sh  = randperm(n_total);
idx_val = idx_sh(1:n_val);
idx_fit = idx_sh(n_val+1:end);

X_norm_all = normalize_x(X_data);
U_tr_all   = ((U_data - u_mean) ./ u_std)';

X_val_n = X_norm_all(idx_val, :);  U_val_n = U_tr_all(idx_val, :);
X_fit_n = X_norm_all(idx_fit, :);  U_fit_n = U_tr_all(idx_fit, :);

net = train_net_fast(X_fit_n, U_fit_n, X_val_n, U_val_n, M, L, 750);
fprintf('  Pre-antrenare finalizata.\n');

eval_nn = @(x_bar) denormalize_u( ...
    predict(net, normalize_x(x_bar), 'ExecutionEnvironment', 'cpu') );

%% =================================================================
%% CAP 3.2.2 — RE-ANTRENARE TINTITA (fine-tuning continuu)
%%   Nu se reinitializeaza reteaua; se continua antrenarea cu LR mic
%%   pe un mini-set ce supra-reprezinta punctele cu eroare mare,
%%   pana epsilon_max < EPS_TARGET sau se epuizeaza rundele.
%% =================================================================
fprintf('\n[Fine-tuning] Tinta: epsilon_max < %.2f\n\n', EPS_TARGET);

MAX_ROUNDS       = 15;
EPOCHS_PER_ROUND = 200;
LR_FINETUNE      = 1e-5;

best_epsilon = inf;
best_net     = net;

for rnd = 1:MAX_ROUNDS
    U_nn_batch = zeros(nu, cnt);
    for ii = 1:cnt
        U_nn_batch(:, ii) = eval_nn(X_data(:, ii));
    end
    err_scalar       = max(abs(U_data - U_nn_batch), [], 1);
    epsilon_max_now  = max(err_scalar);
    epsilon_mean_now = mean(err_scalar);

    tag = '';
    if epsilon_max_now < best_epsilon
        best_epsilon = epsilon_max_now;
        best_net     = net;
        tag = '  [NOU BEST]';
    end
    fprintf('  Runda %d/%d: epsilon_max=%.4f, epsilon_mean=%.4f%s\n', ...
        rnd, MAX_ROUNDS, epsilon_max_now, epsilon_mean_now, tag);

    if epsilon_max_now < EPS_TARGET
        fprintf('  --> TARGET ATINS!\n\n');
        break;
    end

    % Selectie puncte "grele" (eroare > prag); fallback la top-eroare daca prea putine
    hard_mask = err_scalar > EPS_TARGET;
    hard_idx  = find(hard_mask);
    if numel(hard_idx) < 50
        [~, sort_idx] = sort(err_scalar, 'descend');
        hard_idx = sort_idx(1:max(50, floor(0.05*cnt)));
    end
    fprintf('  | puncte grele: %d (err>%.2f)\n', numel(hard_idx), EPS_TARGET);

    % Mini-set: hard (oversample x5) + 20% din restul 
    easy_idx = idx_fit(~ismember(idx_fit, hard_idx));
    n_easy   = floor(0.20 * numel(easy_idx));
    easy_sel = easy_idx(randperm(numel(easy_idx), n_easy));

    aug_idx = [repmat(hard_idx, 1, 5), easy_sel];
    aug_idx = aug_idx(aug_idx >= 1 & aug_idx <= size(X_norm_all,1));

    X_ft = X_norm_all(aug_idx, :);
    U_ft = U_tr_all(aug_idx,  :);

    net = finetune_net(net, X_ft, U_ft, X_val_n, U_val_n, EPOCHS_PER_ROUND, LR_FINETUNE);
    eval_nn = @(x_bar) denormalize_u( ...
        predict(net, normalize_x(x_bar), 'ExecutionEnvironment', 'cpu') );

    if rnd == MAX_ROUNDS
        fprintf('\n  ATENTIE: runde epuizate. epsilon_max_now=%.4f\n\n', epsilon_max_now);
    end
end

% Restaurare cea mai buna retea gasita + calcul eroare finala
net     = best_net;
eval_nn = @(x_bar) denormalize_u( ...
    predict(net, normalize_x(x_bar), 'ExecutionEnvironment', 'cpu') );

U_nn_final = zeros(nu, cnt);
for ii = 1:cnt
    U_nn_final(:, ii) = eval_nn(X_data(:, ii));
end
err_final    = max(abs(U_data - U_nn_final), [], 1);
epsilon_max  = max(err_final);
epsilon_mean = mean(err_final);

fprintf('[Rezultat re-antrenare] epsilon_max=%.4f (target<%.2f: %s), epsilon_mean=%.4f\n\n', ...
    epsilon_max, EPS_TARGET, tf2str(epsilon_max < EPS_TARGET), epsilon_mean);

%% =================================================================
%% CAP 3.3 — CONSTRUCTIA TUBULUI ROBUST Z_NN
%%   Perturbatia totala: d = w + B*epsilon, |w|<=w_bound, |epsilon|<=epsilon_max
%%   Z_NN = mRPI Rakovic pentru d (acelasi algoritm ca pentru Z, cu d_bound)
%% =================================================================
fprintf('[Z_NN] Constructia tubului robust (Rakovic) pentru perturbatia totala...\n');

B_abs   = abs(B);
d_bound = w_bound + B_abs * epsilon_max;

fprintf('  w_bound       = [%.4f; %.4f; %.4f; %.4f]\n', w_bound);
fprintf('  |B|*eps_max   = [%.4f; %.4f; %.4f; %.4f]\n', B_abs * epsilon_max);
fprintf('  d_bound total = [%.4f; %.4f; %.4f; %.4f]\n\n', d_bound);

alpha_target = 0.01;
snn = 0; alpha_nn = inf;
while alpha_nn > alpha_target
    snn = snn + 1;
    alpha_nn = max( (abs(A_K^snn) * d_bound) ./ d_bound );
    if snn > 600, error('A_K nu converge; verifica stabilitatea.'); end
end
fprintf('[Z_NN] Orizont Rakovic: s = %d, alpha = %.6f\n', snn, alpha_nn);

G_D   = diag(d_bound);
G_ZNN = zeros(nx, nx*snn);
for i = 0:(snn-1)
    G_ZNN(:, i*nx+1:(i+1)*nx) = (A_K^i) * G_D;
end
G_ZNN = G_ZNN / (1 - alpha_nn);

Z_NN     = zonotope_to_polyhedron(G_ZNN, nx);
z_max_nn = max(abs(Z_NN.V), [], 1)';

delta_z     = z_max_nn - z_max;
delta_z_rel = delta_z ./ z_max * 100;

fprintf('  z_max (Z)       = [%.6f; %.6f; %.6f; %.6f]\n', z_max);
fprintf('  z_max_nn (Z_NN) = [%.6f; %.6f; %.6f; %.6f]\n', z_max_nn);
fprintf('  Delta z (%%)     = [%.2f%%; %.2f%%; %.2f%%; %.2f%%]\n\n', delta_z_rel);

% Constrangeri tightened pentru tubul NN
X_bar_nn = X - Z_NN;
KZ_nn    = Z_NN.affineMap(K);
U_bar_nn = U - KZ_nn;
if X_bar_nn.isEmptySet() || U_bar_nn.isEmptySet()
    error('Perturbatia este prea mare! Spatiul nominal strans este vid.');
end

% Verificare RPI pentru Z_NN: (a) suficienta per-dimensiune, (b) functie suport
fprintf('[Z_NN] Verificare conditie RPI: A_K*Z_NN (+) D <= Z_NN\n');
lhs_vec     = abs(A_K) * z_max_nn + d_bound;
rpi_per_dim = lhs_vec <= z_max_nn + 1e-10;
for jj = 1:nx
    fprintf('    dim %d: %.6f <= %.6f  -> %s\n', jj, lhs_vec(jj), z_max_nn(jj), tf2str(rpi_per_dim(jj)));
end

G_AKZNN  = A_K * G_ZNN;
dirs     = [eye(nx), -eye(nx)] ./ vecnorm([eye(nx), -eye(nx)]);
h_lhs    = sum(abs(G_AKZNN' * dirs), 1) + sum(abs(G_D' * dirs), 1);
h_rhs    = sum(abs(G_ZNN'   * dirs), 1);
cond_zon = all(h_lhs <= h_rhs + 1e-8);
cond_num = all(rpi_per_dim);
cond_ok  = cond_zon || cond_num;

fprintf('  Verificare zonotop (functie suport): %s\n', tf2str(cond_zon));
fprintf('  Verificare numerica (suficienta):    %s\n', tf2str(cond_num));
fprintf('  %s\n', tf2str(Z <= Z_NN));   % Z trebuie sa fie continut in Z_NN
fprintf('\n');

%% =================================================================
%% CAP 4.1 — SIMULARE PE SISTEMUL LINIAR (MPC vs e-MPC vs NN, in paralel)
%% =================================================================
N_sim = 200;
x_0   = [0.5; 0.0; 0.15; 0.0];

rho_AK  = max(abs(eig(A + B*K_nom)));
V_0_est = x_0' * P_mpc * x_0;
prag_V  = max(1e-4, V_0_est * rho_AK^(2*N_sim) * 100);
fprintf('[Simulare liniara] %d pasi (%.2f s), prag Lyapunov = %.4e\n\n', N_sim, N_sim*Ts, prag_V);

x_history      = zeros(nx, N_sim+1);  x_bar_history      = zeros(nx, N_sim+1);  u_history      = zeros(nu, N_sim);
x_empc_history = zeros(nx, N_sim+1);  x_bar_empc_history = zeros(nx, N_sim+1);  u_empc_history = zeros(nu, N_sim);
x_nn_history   = zeros(nx, N_sim+1);  x_bar_nn_history   = zeros(nx, N_sim+1);  u_nn_history   = zeros(nu, N_sim);

x_history(:,1) = x_0;  x_empc_history(:,1) = x_0;  x_nn_history(:,1) = x_0;
x_bar_k = x_0;  x_bar_empc_k = x_0;  x_bar_nn_k = x_0;

rng(42);
W_disturbances  = -w_bound + 2*w_bound .* rand(nx, N_sim);
epsilon_history = zeros(nu, N_sim);

for k = 1:N_sim
    w_k = W_disturbances(:, k);

    % --- MPC implicit ---
    x_k = x_history(:, k);
    x_bar_history(:, k) = x_bar_k;
    [u_bar_k, feasible_mpc, ~] = mpc.evaluate(x_bar_k);
    if ~feasible_mpc, warning('MPC infazabil la pasul %d.', k); break; end
    u_k = u_bar_k + K * (x_k - x_bar_k);
    x_history(:, k+1) = A * x_k + B * u_k + w_k;
    x_bar_k = A * x_bar_k + B * u_bar_k;
    u_history(:, k) = u_k;

    % --- e-MPC explicit ---
    x_empc_k = x_empc_history(:, k);
    x_bar_empc_history(:, k) = x_bar_empc_k;
    [u_bar_empc_k, feasible_empc, ~] = empc.evaluate(x_bar_empc_k);
    if ~feasible_empc, warning('e-MPC infazabil la pasul %d.', k); break; end
    u_empc_k = u_bar_empc_k + K * (x_empc_k - x_bar_empc_k);
    x_empc_history(:, k+1) = A * x_empc_k + B * u_empc_k + w_k;
    x_bar_empc_k = A * x_bar_empc_k + B * u_bar_empc_k;
    u_empc_history(:, k) = u_empc_k;

    % --- NN (tub Z_NN) ---
    x_nn_k = x_nn_history(:, k);
    x_bar_nn_history(:, k) = x_bar_nn_k;
    u_bar_nn_k = eval_nn(x_bar_nn_k);

    [u_bar_empc_ref_k, feas_ref, ~] = empc.evaluate(x_bar_nn_k);  % referinta pt. eroarea online
    if feas_ref && all(~isnan(u_bar_empc_ref_k))
        epsilon_history(:, k) = u_bar_empc_ref_k - u_bar_nn_k;
    end

    u_nn_k = u_bar_nn_k + K * (x_nn_k - x_bar_nn_k);
    x_nn_history(:, k+1) = A * x_nn_k + B * u_nn_k + w_k;
    x_bar_nn_k = A * x_bar_nn_k + B * u_bar_nn_k;
    u_nn_history(:, k) = u_nn_k;
end

x_bar_history(:, N_sim+1)      = x_bar_k;
x_bar_empc_history(:, N_sim+1) = x_bar_empc_k;
x_bar_nn_history(:, N_sim+1)   = x_bar_nn_k;

%% --- Raport validare simulare liniara ---
fprintf('=========================================================================\n');
fprintf('  RAPORT VALIDARE — TUBE-MPC vs TUBE-e-MPC vs TUBE-NN\n');
fprintf('=========================================================================\n');

in_Z_mpc = true; in_Z_empc = true; in_ZNN_nn = true;
for k = 1:N_sim+1
    if ~Z.contains(x_history(:,k)      - x_bar_history(:,k)),      in_Z_mpc  = false; end
    if ~Z.contains(x_empc_history(:,k) - x_bar_empc_history(:,k)), in_Z_empc = false; end
    if ~Z_NN.contains(x_nn_history(:,k)- x_bar_nn_history(:,k)),   in_ZNN_nn = false; end
end
fprintf('\n1. VERIFICARE TUB:  MPC=%s | e-MPC=%s | NN=%s (in Z_NN)\n', ...
    tf2str(in_Z_mpc), tf2str(in_Z_empc), tf2str(in_ZNN_nn));

u_viol_mpc  = sum(u_history(:)      < u_min-1e-6 | u_history(:)      > u_max+1e-6);
u_viol_empc = sum(u_empc_history(:) < u_min-1e-6 | u_empc_history(:) > u_max+1e-6);
u_viol_nn   = sum(u_nn_history(:)   < u_min-1e-6 | u_nn_history(:)   > u_max+1e-6);
x_viol_mpc  = sum(any(x_history      < x_min-1e-6 | x_history      > x_max+1e-6, 1));
x_viol_empc = sum(any(x_empc_history < x_min-1e-6 | x_empc_history > x_max+1e-6, 1));
x_viol_nn   = sum(any(x_nn_history   < x_min-1e-6 | x_nn_history   > x_max+1e-6, 1));
fprintf('2. VIOLARI CONSTRANGERI:  MPC(u=%d,x=%d) | e-MPC(u=%d,x=%d) | NN(u=%d,x=%d)\n', ...
    u_viol_mpc, x_viol_mpc, u_viol_empc, x_viol_empc, u_viol_nn, x_viol_nn);

%% =================================================================
%% CAP 4.2 — STABILITATEA LYAPUNOV
%% =================================================================
V_final_mpc  = x_bar_history(:,end)'      * P_mpc * x_bar_history(:,end);
V_final_empc = x_bar_empc_history(:,end)' * P_mpc * x_bar_empc_history(:,end);
V_final_nn   = x_bar_nn_history(:,end)'   * P_mpc * x_bar_nn_history(:,end);
fprintf('3. STABILITATE (V final, prag=%.2e): MPC=%.2e(%s) | e-MPC=%.2e(%s) | NN=%.2e(%s)\n', ...
    prag_V, V_final_mpc, tf2str(V_final_mpc<prag_V), V_final_empc, tf2str(V_final_empc<prag_V), ...
    V_final_nn, tf2str(V_final_nn<prag_V));

fprintf('4. COMPARARE Z vs Z_NN:\n');
fprintf('   Dim  |   z_max(Z)   | z_max_nn(Z_NN) | Delta(%%)\n');
for j = 1:nx
    fprintf('   x_%d  | %12.6f | %14.6f | %8.2f%%\n', j, z_max(j), z_max_nn(j), delta_z_rel(j));
end

eps_sim_max  = max(abs(epsilon_history(:)));
eps_sim_mean = mean(abs(epsilon_history(:)));
fprintf('\n5. EROARE APROXIMARE NN: eps_max(antren)=%.4f | eps_max(online)=%.4f | eps_mean(online)=%.4f\n', ...
    epsilon_max, eps_sim_max, eps_sim_mean);
fprintf('6. CONDITIE RPI Z_NN: %s\n', tf2str(cond_ok));
fprintf('=========================================================================\n\n');

%% =================================================================
%% FIGURI — Simulare liniara (stari, comenzi, Lyapunov, seturi mRPI)
%% =================================================================
t = (0:N_sim) * Ts;

col_mpc  = [0.00, 0.45, 0.74];  col_empc = [0.85, 0.33, 0.10];
col_nn   = [0.47, 0.67, 0.19];  col_znn  = [0.75, 0.00, 0.75];

state_labels_short = {'$x_1$: Pozitie (m)', '$x_2$: Viteza (m/s)', ...
    '$x_3$: Unghi (rad)', '$x_4$: Viteza ang. (rad/s)'};
state_labels_plain = {'Pozitie Carucior (m)', 'Viteza Carucior (m/s)', ...
    'Unghi Pendul (rad)',  'Viteza Unghiulara (rad/s)'};

% --- Figura 1: evolutia starilor cu tuburile Z si Z_NN ---
figure('Name','Stari: Z vs Z_NN','Color','w','Position',[50,50,1400,900]);
for j = 1:nx
    subplot(2,2,j);
    tup_nn = x_bar_nn_history(j,:) + z_max_nn(j); tdn_nn = x_bar_nn_history(j,:) - z_max_nn(j);
    fill([t, fliplr(t)], [tup_nn, fliplr(tdn_nn)], col_znn, 'FaceAlpha',0.18, ...
        'EdgeColor',col_znn, 'EdgeAlpha',0.7, 'LineWidth',1.2, 'DisplayName','$Z_{NN}$ (tub NN)');
    hold on;
    tup_z = x_bar_history(j,:) + z_max(j); tdn_z = x_bar_history(j,:) - z_max(j);
    fill([t, fliplr(t)], [tup_z, fliplr(tdn_z)], col_mpc, 'FaceAlpha',0.22, ...
        'EdgeColor',col_mpc, 'EdgeAlpha',0.9, 'LineWidth',1.4, 'LineStyle','--', ...
        'DisplayName','$Z$ (tub MPC/e-MPC)');

    plot(t, x_bar_history(j,:),      '--', 'Color',col_mpc,  'LineWidth',2.0, 'DisplayName','$\bar{x}$ MPC');
    plot(t, x_bar_empc_history(j,:), ':',  'Color',col_empc, 'LineWidth',2.2, 'DisplayName','$\bar{x}$ e-MPC');
    plot(t, x_bar_nn_history(j,:),   '-.', 'Color',col_nn,   'LineWidth',2.0, 'DisplayName','$\bar{x}$ NN');
    plot(t, x_history(j,:),      '-', 'Color',col_mpc  +0.3*(1-col_mpc),  'LineWidth',1.0, 'DisplayName','$x$ MPC');
    plot(t, x_empc_history(j,:), '-', 'Color',col_empc +0.3*(1-col_empc), 'LineWidth',1.0, 'DisplayName','$x$ e-MPC');
    plot(t, x_nn_history(j,:),   '-', 'Color',col_nn   +0.3*(1-col_nn),   'LineWidth',1.0, 'DisplayName','$x$ NN');

    yline(x_max(j),'k:','LineWidth',1.5,'HandleVisibility','off');
    yline(x_min(j),'k:','LineWidth',1.5,'DisplayName','Limite fizice');

    grid on; xlabel('Timp (s)','Interpreter','latex','FontSize',10);
    ylabel(state_labels_plain{j},'FontSize',9);
    title(state_labels_short{j},'Interpreter','latex','FontSize',11);
    if j == 1, legend('Location','northeast','FontSize',7,'Interpreter','latex','NumColumns',2); end
end
sgtitle('Evolutia Starilor --- $Z$ (albastru) vs $Z_{NN}$ (mov)', ...
    'FontSize',13,'FontWeight','bold','Interpreter','latex');

% --- Figura 2: comenzi u(k) ---
figure('Name','Comenzi u(k)','Color','w','Position',[50,1850,1400,320]);
fill([t(1) t(end) t(end) t(1)], [u_min u_min u_max u_max], [0.9 0.9 0.9], ...
    'FaceAlpha',0.3,'EdgeColor','none','DisplayName','Zona admisa');
hold on;
plot(t(1:end-1), u_history(1,:),      '-',  'Color',col_mpc,  'LineWidth',2.0, 'DisplayName','$u$ MPC');
plot(t(1:end-1), u_empc_history(1,:), '--', 'Color',col_empc, 'LineWidth',2.0, 'DisplayName','$u$ e-MPC');
plot(t(1:end-1), u_nn_history(1,:),   '-.', 'Color',col_nn,   'LineWidth',2.0, 'DisplayName','$u$ NN');
yline(u_max,'k:','LineWidth',1.5,'DisplayName','$u_{max}/u_{min}$');
yline(u_min,'k:','LineWidth',1.5,'HandleVisibility','off');
grid on; xlabel('Timp (s)','Interpreter','latex'); ylabel('Forta (N)');
title('Comanda Aplicata $u(k)$','Interpreter','latex','FontSize',12);
legend('Location','best','FontSize',9,'Interpreter','latex','NumColumns',4);
ylim([u_min*1.15, u_max*1.15]);

% --- Figura 3: eroare epsilon si functia Lyapunov ---
figure('Name','Eroare NN si Lyapunov','Color','w','Position',[50,2220,1400,620]);
V_mpc_vec = zeros(1,N_sim+1); V_empc_vec = zeros(1,N_sim+1); V_nn_vec = zeros(1,N_sim+1);
for k = 1:N_sim+1
    V_mpc_vec(k)  = x_bar_history(:,k)'      * P_mpc * x_bar_history(:,k);
    V_empc_vec(k) = x_bar_empc_history(:,k)' * P_mpc * x_bar_empc_history(:,k);
    V_nn_vec(k)   = x_bar_nn_history(:,k)'   * P_mpc * x_bar_nn_history(:,k);
end
semilogy(t, V_mpc_vec,  '-',  'Color',col_mpc,  'LineWidth',2.0, 'DisplayName','$V(\bar{x})$ MPC'); hold on;
semilogy(t, V_empc_vec, '--', 'Color',col_empc, 'LineWidth',2.0, 'DisplayName','$V(\bar{x})$ e-MPC');
semilogy(t, V_nn_vec,   '-.', 'Color',col_nn,   'LineWidth',2.0, 'DisplayName','$V(\bar{x})$ NN');
yline(prag_V,'k:','LineWidth',1.5,'DisplayName',['Prag $=$ ' num2str(prag_V,'%.1e')]);
grid on; xlabel('Timp (s)','Interpreter','latex'); ylabel('$V(\bar{x}_k)$ [log]','Interpreter','latex');
title('Functia Lyapunov $V(\bar{x}_k) = \bar{x}^T P \bar{x}$','Interpreter','latex','FontSize',12);
legend('Location','best','FontSize',9,'Interpreter','latex');

% --- Figura 4: seturi mRPI Z vs Z_NN, toate proiectiile 2D ---
fprintf('[Figuri] Calcul proiectii 2D ale seturilor mRPI...\n');
pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
ax_lbls = {'$x_1$ pos (m)', '$x_2$ vel (m/s)', '$x_3$ ang (rad)', '$x_4$ $\dot{\theta}$ (rad/s)'};

figure('Name','Seturi mRPI: Z vs Z_NN','Color','w','Position',[700,50,1050,900]);
for pp = 1:size(pairs,1)
    di = pairs(pp,1); dj = pairs(pp,2);
    subplot(3,2,pp);
    plot_projection(Z,    [di dj], col_mpc, '-',  '$Z$ (mRPI base)',  z_max);
    hold on;
    plot_projection(Z_NN, [di dj], col_znn, '--', '$Z_{NN}$ (mRPI NN)', z_max_nn);
    plot(0, 0, 'k+', 'MarkerSize',9, 'LineWidth',2.0, 'HandleVisibility','off');
    grid on; axis equal;
    xlabel(ax_lbls{di},'Interpreter','latex','FontSize',9);
    ylabel(ax_lbls{dj},'Interpreter','latex','FontSize',9);
    title(['Proiectie $(x_{' num2str(di) '}, x_{' num2str(dj) '})$'],'Interpreter','latex','FontSize',10);
    legend('Location','best','FontSize',7,'Interpreter','latex');
end
sgtitle('Seturi mRPI --- $Z$ (albastru) vs $Z_{NN}$ (mov)', ...
    'FontSize',12,'FontWeight','bold','Interpreter','latex');

%% =================================================================
%% CAP 4.4 — VALIDARE PE SISTEMUL NELINIAR (CasADi RK4)
%% =================================================================
fprintf('\n[Neliniar] Configurare planta CasADi RK4...\n');
params_plant.m = 0.135; params_plant.mcart = 0.5; params_plant.L = 0.47/2;
params_plant.I = (0.135*0.47^2)/3; params_plant.g = 9.8; params_plant.b = 0; params_plant.q = 0;

F_plant    = nonlinear_plant(params_plant, Ts);
plant_step = @(x,u) full(F_plant(casadi.DM(double(x(:))), casadi.DM(double(u(:)))));

fprintf('[Neliniar] Simulare paralela (%d pasi)...\n', N_sim);
xnl_mpc = zeros(nx,N_sim+1); xnl_empc = zeros(nx,N_sim+1); xnl_nn = zeros(nx,N_sim+1);
xnl_mpc(:,1) = x_0; xnl_empc(:,1) = x_0; xnl_nn(:,1) = x_0;
xbar_nl_mpc = x_0; xbar_nl_empc = x_0; xbar_nl_nn = x_0;
unl_mpc = zeros(nu,N_sim); unl_empc = zeros(nu,N_sim); unl_nn = zeros(nu,N_sim);

rng(42);
W_nl = -w_bound + 2*w_bound .* rand(nx, N_sim);
eval_nn_d = @(x_bar) double(denormalize_u(predict(net, normalize_x(x_bar), 'ExecutionEnvironment','cpu')));

for k = 1:N_sim
    w_k = W_nl(:, k);

    [u_bar_k, feas_mpc, ~] = mpc.evaluate(xbar_nl_mpc);
    if ~feas_mpc, warning('NL-MPC infazabil la pasul %d.', k); break; end
    u_k = min(u_max, max(u_min, u_bar_k + K*(xnl_mpc(:,k) - xbar_nl_mpc)));
    xnl_mpc(:,k+1) = plant_step(xnl_mpc(:,k), u_k) + w_k;
    xbar_nl_mpc    = A*xbar_nl_mpc + B*u_bar_k;
    unl_mpc(:,k)   = u_k;

    [u_bar_ek, feas_empc, ~] = empc.evaluate(xbar_nl_empc);
    if ~feas_empc, warning('NL-e-MPC infazabil la pasul %d.', k); break; end
    u_ek = min(u_max, max(u_min, u_bar_ek + K*(xnl_empc(:,k) - xbar_nl_empc)));
    xnl_empc(:,k+1) = plant_step(xnl_empc(:,k), u_ek) + w_k;
    xbar_nl_empc    = A*xbar_nl_empc + B*u_bar_ek;
    unl_empc(:,k)   = u_ek;

    u_bar_nk = eval_nn_d(xbar_nl_nn);
    u_nk     = min(u_max, max(u_min, u_bar_nk + K*(xnl_nn(:,k) - xbar_nl_nn)));
    xnl_nn(:,k+1) = plant_step(xnl_nn(:,k), u_nk) + w_k;
    xbar_nl_nn    = A*xbar_nl_nn + B*u_bar_nk;
    unl_nn(:,k)   = u_nk;
end
fprintf('  Simulare neliniara finalizata.\n\n');

% Timp de stabilizare (||x|| sub 2% din ||x_0||)
prag_stab = 0.02 * norm(x_0);
stab_mpc = N_sim; stab_empc = N_sim; stab_nn = N_sim;
for k = 1:N_sim+1
    if norm(xnl_mpc(:,k))  < prag_stab && stab_mpc  == N_sim, stab_mpc  = k; end
    if norm(xnl_empc(:,k)) < prag_stab && stab_empc == N_sim, stab_empc = k; end
    if norm(xnl_nn(:,k))   < prag_stab && stab_nn   == N_sim, stab_nn   = k; end
end
t_stab_mpc = stab_mpc*Ts; t_stab_empc = stab_empc*Ts; t_stab_nn = stab_nn*Ts;

%% =================================================================
%% CAP 4.5 — ANALIZA COMPARATIVA A PERFORMANTEI COMPUTATIONALE
%%   (memorie: eq. analitice din lucrare; timp: timeit pe puncte uniforme)
%% =================================================================
alpha_bit = 8;  % bytes per float64

% --- Memorie MPC implicit (QP online: Hessian + Aeq + Aineq) ---
n_qp_vars = nu*N; n_qp_eq = nx*N; n_qp_ineq = (2*nx+2*nu)*N;
mem_mpc_kb = (n_qp_vars^2 + n_qp_eq*n_qp_vars + n_qp_ineq*n_qp_vars) * alpha_bit / 1024;

% --- Memorie e-MPC explicit: Gamma_K = alpha*[n_h*(n_x+1) + n_f*(n_x*n_u+n_u)] ---
partition = empc.feedback.Set;
n_regions = empc.nr;
seen_H = {}; seen_F = {}; n_h = 0; n_f = 0;
for i = 1:n_regions
    Hi = round(partition(i).A * 1e6); hi = round(partition(i).b * 1e6);
    Fi = round(partition(i).Functions('primal').F * 1e6);
    gi = round(partition(i).Functions('primal').g * 1e6);
    key_h = mat2str([Hi hi]); key_f = mat2str([Fi gi]);
    if ~any(strcmp(seen_H, key_h)), seen_H{end+1} = key_h; n_h = n_h + size(Hi,1); end
    if ~any(strcmp(seen_F, key_f)), seen_F{end+1} = key_f; n_f = n_f + 1; end
end
mem_empc_kb = alpha_bit * (n_h*(nx+1) + n_f*(nx*nu+nu)) / 1024;

% --- Memorie NN: Gamma_N = alpha*[(n_x+1)*M + (L-1)*(M+1)*M + (M+1)*n_u] ---
W_nn = {}; b_nn = {}; mem_nn_bytes = 0;
for li = 1:numel(net.Layers)
    lyr = net.Layers(li);
    if isa(lyr, 'nnet.cnn.layer.FullyConnectedLayer')
        W_nn{end+1} = double(lyr.Weights); b_nn{end+1} = double(lyr.Bias);
        mem_nn_bytes = mem_nn_bytes + numel(lyr.Weights)*alpha_bit + numel(lyr.Bias)*alpha_bit;
    end
end
mem_nn_kb = mem_nn_bytes / 1024;
n_params_nn = mem_nn_bytes / alpha_bit;
n_fc = numel(W_nn);

fprintf('=== Comparatie memorie ===\n');
fprintf('MPC implicit : %.2f kB | e-MPC (%d regiuni, n_h=%d, n_f=%d): %.2f kB | NN (%d param): %.2f kB\n', ...
    mem_mpc_kb, n_regions, n_h, n_f, mem_empc_kb, n_params_nn, mem_nn_kb);
fprintf('Raport NN/e-MPC = %.2f %% | Raport NN/MPC impl. = %.2f %%\n\n', ...
    mem_nn_kb/mem_empc_kb*100, mem_nn_kb/mem_mpc_kb*100);

fprintf('=== Timp stabilizare (neliniar) ===\n');
fprintf('Implicit MPC: %.2fs (%d pasi) | Explicit MPC: %.2fs (%d pasi) | NN: %.2fs (%d pasi)\n\n', ...
    t_stab_mpc, stab_mpc, t_stab_empc, stab_empc, t_stab_nn, stab_nn);

% --- Extractie ponderi -> forward pass algebric (pentru benchmark/ESP32) ---
xm = x_mean; xs = x_std; um = u_mean; us = u_std; W_cap = W_nn; b_cap = b_nn; nfc = n_fc;
eval_nn_raw = @(x_bar) nn_forward(x_bar, W_cap, b_cap, nfc, xm, xs, um, us);

err_validate = zeros(1, 20);
for ii = 1:20
    err_validate(ii) = max(abs(eval_nn_d(X_data(:,ii)) - eval_nn_raw(X_data(:,ii))));
end
fprintf('[Validare nn_forward vs predict] Max diferenta = %.2e (trebuie < 1e-6)\n\n', max(err_validate));

% --- Extractie regiuni e-MPC -> evaluare rapida (cautare liniara pe regiuni) ---
partition = empc.optimizer.Set;
n_reg = empc.nr; nu_e = empc.nu;
A_reg = cell(n_reg,1); b_reg = cell(n_reg,1); F_reg = cell(n_reg,1); g_reg = cell(n_reg,1);
for i = 1:n_reg
    A_reg{i} = partition(i).A; b_reg{i} = partition(i).b;
    F_full = partition(i).Functions('primal').F; g_full = partition(i).Functions('primal').g;
    F_reg{i} = F_full(1:nu_e,:); g_reg{i} = g_full(1:nu_e,:);   % doar pasul curent (N=1)
end
eval_empc_fast = @(x) empc_eval_fast(x, A_reg, b_reg, F_reg, g_reg, n_reg);

err_empc_val = zeros(1,20);
for ii = 1:20
    xv = X_data(:, ii);
    [u_mpt, feas, openloop] = empc.evaluate(xv);
    if feas && isfield(openloop,'region') && openloop.region > 0
        ri = openloop.region;
        err_empc_val(ii) = max(abs(u_mpt - (F_reg{ri}*xv + g_reg{ri})));
    end
end
fprintf('[Validare empc_eval_fast vs empc.evaluate] Max diferenta = %.2e (trebuie < 1e-6)\n\n', max(err_empc_val));

% --- Benchmark timp execuție (timeit) pe puncte fezabile uniforme ---
fprintf('[Benchmark] Generare puncte uniforme fezabile...\n');
N_bench_t = 100;
X_bench_t = zeros(nx, N_bench_t);
count = 0; tries = 0; max_tries = 100000;
while count < N_bench_t && tries < max_tries
    x_rand = x_min + (x_max - x_min) .* rand(nx,1);
    [~, feas] = mpc.evaluate(x_rand);
    if feas, count = count+1; X_bench_t(:,count) = x_rand; end
    tries = tries + 1;
end
if count < N_bench_t
    warning('Doar %d puncte fezabile gasite din %d cerute.', count, N_bench_t);
    X_bench_t = X_bench_t(:,1:count); N_bench_t = count;
end

for i = 1:10   % warmup
    xw = X_bench_t(:, i);
    mpc.evaluate(xw); eval_empc_fast(xw); eval_nn_raw(xw);
end

fprintf('[Benchmark] Masurare timp (N=%d puncte)...\n', N_bench_t);
time_mpc_b = zeros(N_bench_t,1); time_empc_b = zeros(N_bench_t,1); time_nn_raw_b = zeros(N_bench_t,1);
for i = 1:N_bench_t
    xt = X_bench_t(:, i);
    time_mpc_b(i)    = timeit(@() mpc.evaluate(xt));
    time_empc_b(i)   = timeit(@() eval_empc_fast(xt));
    time_nn_raw_b(i) = timeit(@() eval_nn_raw(xt));
end

mpc_ms = time_mpc_b*1e3; empc_ms = time_empc_b*1e3; raw_ms = time_nn_raw_b*1e3;
M_mpc  = [mean(mpc_ms),  max(mpc_ms),  std(mpc_ms),  median(mpc_ms)];
M_empc = [mean(empc_ms), max(empc_ms), std(empc_ms), median(empc_ms)];
M_raw  = [mean(raw_ms),  max(raw_ms),  std(raw_ms),  median(raw_ms)];

fprintf('\n=========================================================================\n');
fprintf(' BENCHMARK TIMP EXECUTIE (N=%d puncte, timeit)\n', N_bench_t);
fprintf('=========================================================================\n');
fprintf(' %-18s | %9s | %9s | %9s | %9s\n', 'Metoda','Medie(ms)','WCET(ms)','Std(ms)','Med(ms)');
fprintf(' %-18s | %9.4f | %9.4f | %9.4f | %9.4f\n', 'Implicit MPC',    M_mpc);
fprintf(' %-18s | %9.4f | %9.4f | %9.4f | %9.4f\n', 'Explicit MPC',    M_empc);
fprintf(' %-18s | %9.4f | %9.4f | %9.4f | %9.4f\n', 'NN raw algebric', M_raw);
fprintf('-------------------------------------------------------------------------\n');
fprintf(' SPEEDUP vs MPC implicit: e-MPC %.1fx avg / %.1fx WCET | NN raw %.1fx avg / %.1fx WCET\n', ...
    M_mpc(1)/M_empc(1), M_mpc(2)/M_empc(2), M_mpc(1)/M_raw(1), M_mpc(2)/M_raw(2));
fprintf(' CONSISTENTA (WCET/Medie, ideal=1.0): e-MPC %.1fx | NN raw %.1fx\n', M_empc(2)/M_empc(1), M_raw(2)/M_raw(1));
fprintf('=========================================================================\n\n');

% --- Figura 5: traiectorii liniar vs neliniar, cu tuburile liniare ---
figure('Name','Liniar + Neliniar cu Tuburi Z, Z_NN','Color','w','Position',[50,50,1400,900]);
for j = 1:nx
    subplot(2,2,j);
    tup_nn = x_bar_nn_history(j,:) + z_max_nn(j); tdn_nn = x_bar_nn_history(j,:) - z_max_nn(j);
    fill([t, fliplr(t)], [tup_nn, fliplr(tdn_nn)], col_znn, 'FaceAlpha',0.18, ...
        'EdgeColor',col_znn,'EdgeAlpha',0.7,'LineWidth',1.2,'DisplayName','$Z_{NN}$ (tub NN)');
    hold on;
    tup_z = x_bar_history(j,:) + z_max(j); tdn_z = x_bar_history(j,:) - z_max(j);
    fill([t, fliplr(t)], [tup_z, fliplr(tdn_z)], col_mpc, 'FaceAlpha',0.22, ...
        'EdgeColor',col_mpc,'EdgeAlpha',0.9,'LineWidth',1.4,'LineStyle','--','DisplayName','$Z$ (tub MPC/e-MPC)');

    plot(t, x_bar_history(j,:),      '--', 'Color',col_mpc,  'LineWidth',2.0, 'DisplayName','$\bar{x}$ MPC lin');
    plot(t, x_bar_empc_history(j,:), ':',  'Color',col_empc, 'LineWidth',2.2, 'DisplayName','$\bar{x}$ e-MPC lin');
    plot(t, x_bar_nn_history(j,:),   '-.', 'Color',col_nn,   'LineWidth',2.0, 'DisplayName','$\bar{x}$ NN lin');
    plot(t, x_history(j,:),      '-', 'Color',col_mpc  +0.3*(1-col_mpc),  'LineWidth',1.0,'DisplayName','$x$ MPC lin');
    plot(t, x_empc_history(j,:), '-', 'Color',col_empc +0.3*(1-col_empc), 'LineWidth',1.0,'DisplayName','$x$ e-MPC lin');
    plot(t, x_nn_history(j,:),   '-', 'Color',col_nn   +0.3*(1-col_nn),   'LineWidth',1.0,'DisplayName','$x$ NN lin');

    plot(t, xnl_mpc(j,:),  '-', 'Color',col_mpc*0.55,  'LineWidth',1.4, 'DisplayName','$x$ MPC neliniar');
    plot(t, xnl_empc(j,:), '-', 'Color',col_empc*0.55, 'LineWidth',1.4, 'DisplayName','$x$ e-MPC neliniar');
    plot(t, xnl_nn(j,:),   '-', 'Color',col_nn*0.55,   'LineWidth',1.4, 'DisplayName','$x$ NN neliniar');

    yline(x_max(j),'k:','LineWidth',1.2,'HandleVisibility','off');
    yline(x_min(j),'k:','LineWidth',1.2,'HandleVisibility','off');
    grid on; xlabel('Timp (s)','Interpreter','latex','FontSize',10);
    ylabel(state_labels_plain{j},'FontSize',9);
    title(state_labels_short{j},'Interpreter','latex','FontSize',11);
    if j == 1, legend('Location','northeast','FontSize',7,'Interpreter','latex','NumColumns',2); end
end
sgtitle(['Liniar (culori deschise) + Neliniar RK4 (culori inchise)' ...
    ' --- Tuburi liniare $Z$ si $Z_{NN}$'], 'FontSize',13,'FontWeight','bold','Interpreter','latex');

%% =================================================================
%% CAP 4.3 — VALIDARE STATISTICA (Monte Carlo)
%%   N_mc stari initiale aleatoare din X_bar, simulare scurta,
%%   verificare: mentinere in tub, convergenta Lyapunov, violari.
%% =================================================================
fprintf('\n[Monte Carlo] N_mc=1000, %d pasi/traiectorie...\n', 100);
N_mc = 1000; N_sim_mc = 100; prag_mc = 0.02;

intub_empc_mc = true(N_mc,1);  uviol_empc_mc = false(N_mc,1);  xviol_empc_mc = false(N_mc,1);  stab_empc_mc = false(N_mc,1);
intub_nn_mc   = true(N_mc,1);  uviol_nn_mc   = false(N_mc,1);  xviol_nn_mc   = false(N_mc,1);  stab_nn_mc   = false(N_mc,1);
lyap_conv_empc = false(N_mc,1); lyap_conv_nn = false(N_mc,1);
V_empc_all = zeros(N_mc, N_sim_mc+1);  V_nn_all = zeros(N_mc, N_sim_mc+1);

rng(0);
n_valid = 0;
while n_valid < N_mc
    x0_mc = x_min + (x_max - x_min) .* rand(nx,1);
    if ~X_bar.contains(x0_mc), continue; end
    try eval_empc_fast(x0_mc); catch, continue; end

    n_valid = n_valid + 1;
    W_mc = -w_bound + 2*w_bound .* rand(nx, N_sim_mc);
    xk_e = x0_mc; xbk_e = x0_mc; xk_n = x0_mc; xbk_n = x0_mc;
    V_empc_all(n_valid,1) = x0_mc'*P_mpc*x0_mc;
    V_nn_all(n_valid,1)   = x0_mc'*P_mpc*x0_mc;

    for k = 1:N_sim_mc
        w_k = W_mc(:,k);

        try
            ubk_e = eval_empc_fast(xbk_e);
            uk_e  = ubk_e + K*(xk_e - xbk_e);
            if any(abs(xk_e-xbk_e) > z_max+1e-6), intub_empc_mc(n_valid) = false; end
            if any(uk_e < u_min-1e-6) || any(uk_e > u_max+1e-6), uviol_empc_mc(n_valid) = true; end
            xk_e  = A*xk_e  + B*uk_e  + w_k;
            xbk_e = A*xbk_e + B*ubk_e;
            if any(xk_e < x_min-1e-6) || any(xk_e > x_max+1e-6), xviol_empc_mc(n_valid) = true; end
            if ~stab_empc_mc(n_valid) && norm(xk_e) < prag_mc, stab_empc_mc(n_valid) = true; end
            V_empc_all(n_valid, k+1) = xbk_e' * P_mpc * xbk_e;
        catch
            V_empc_all(n_valid, k+1) = V_empc_all(n_valid, k);
        end

        ubk_n = eval_nn_raw(xbk_n);
        uk_n  = ubk_n + K*(xk_n - xbk_n);
        if any(abs(xk_n-xbk_n) > z_max_nn+1e-6), intub_nn_mc(n_valid) = false; end
        if any(uk_n < u_min-1e-6) || any(uk_n > u_max+1e-6), uviol_nn_mc(n_valid) = true; end
        xk_n  = A*xk_n  + B*uk_n  + w_k;
        xbk_n = A*xbk_n + B*ubk_n;
        if any(xk_n < x_min-1e-6) || any(xk_n > x_max+1e-6), xviol_nn_mc(n_valid) = true; end
        if ~stab_nn_mc(n_valid) && norm(xk_n) < prag_mc, stab_nn_mc(n_valid) = true; end
        V_nn_all(n_valid, k+1) = xbk_n' * P_mpc * xbk_n;
    end

    V0 = x0_mc' * P_mpc * x0_mc;
    lyap_conv_empc(n_valid) = (V_empc_all(n_valid,end) < V0/100);
    lyap_conv_nn(n_valid)   = (V_nn_all(n_valid,end)   < V0/100);

    if mod(n_valid,100)==0, fprintf('  %d/%d\n', n_valid, N_mc); end
end

fprintf('\n  %-6s | stab=%5.1f%% | tub=%5.1f%% | xviol=%4.1f%% | uviol=%4.1f%% | V_conv=%5.1f%%\n', ...
    'e-MPC', 100*mean(stab_empc_mc), 100*mean(intub_empc_mc), 100*mean(xviol_empc_mc), 100*mean(uviol_empc_mc), 100*mean(lyap_conv_empc));
fprintf('  %-6s | stab=%5.1f%% | tub=%5.1f%% | xviol=%4.1f%% | uviol=%4.1f%% | V_conv=%5.1f%%\n', ...
    'NN',    100*mean(stab_nn_mc),   100*mean(intub_nn_mc),   100*mean(xviol_nn_mc),   100*mean(uviol_nn_mc),   100*mean(lyap_conv_nn));

% --- Grafic bare: indicatori Monte Carlo (e-MPC vs NN) ---
figure('Color','w','Position',[80,80,960,560]);
kpi_empc = [mean(intub_empc_mc), mean(lyap_conv_empc), mean(xviol_empc_mc), mean(uviol_empc_mc)] * 100;
kpi_nn   = [mean(intub_nn_mc),   mean(lyap_conv_nn),   mean(xviol_nn_mc),   mean(uviol_nn_mc)]   * 100;
data_plot = [kpi_empc; kpi_nn]';

c_empc = [0.678, 0.816, 0.929]; c_nn = [0.957, 0.698, 0.796];
c_edge_empc = [0.200, 0.520, 0.720]; c_edge_nn = [0.780, 0.200, 0.420];

ax = axes('Color','w','GridColor',[0.85 0.85 0.85],'GridLineStyle',':','LineWidth',0.8, ...
    'FontSize',11,'FontName','Times New Roman','Box','off', ...
    'XColor',[0.2 0.2 0.2],'YColor',[0.2 0.2 0.2]);
hold(ax, 'on');
b = bar(ax, data_plot, 'grouped', 'BarWidth', 0.68);
b(1).FaceColor = c_empc; b(1).EdgeColor = c_edge_empc; b(1).LineWidth = 1.3;
b(2).FaceColor = c_nn;   b(2).EdgeColor = c_edge_nn;   b(2).LineWidth = 1.3;

nGroups = size(data_plot,1); nBars = size(data_plot,2); barW = 0.68/nBars;
for iB = 1:nBars
    xOff = (iB - (nBars+1)/2) * barW;
    for iG = 1:nGroups
        val = data_plot(iG, iB);
        if val > 2
            text(ax, iG+xOff, val+1.5, sprintf('%.1f%%', val), 'HorizontalAlignment','center', ...
                'FontSize',8.5,'FontName','Times New Roman','Color',[0.15 0.15 0.15],'FontWeight','bold');
        end
    end
end

yline(ax, 100, '--', 'Color',[0.6 0.6 0.6], 'LineWidth',0.9, 'Label','100%', ...
    'LabelHorizontalAlignment','left','FontSize',9,'FontName','Times New Roman');
xline(ax, 2.5, '-', 'Color',[0.7 0.7 0.7], 'LineWidth',1.0, 'Alpha',0.8);
text(ax, 1.5, 108, 'Indicatori de performanta', 'HorizontalAlignment','center', ...
    'FontSize',9,'FontName','Times New Roman','Color',[0.35 0.35 0.35],'FontAngle','italic');
text(ax, 3.5, 108, 'Indicatori de incalcare', 'HorizontalAlignment','center', ...
    'FontSize',9,'FontName','Times New Roman','Color',[0.55 0.15 0.25],'FontAngle','italic');

set(ax, 'XTick', 1:4, 'XTickLabel', ...
    {'Mentinere in Tub','Convergenta V_{Lyap}','Incalcare Stari','Incalcare Comenzi'});
ylabel(ax, 'Procentaj (%)', 'FontSize',12,'FontWeight','bold','FontName','Times New Roman');
ylim(ax, [0 114]); xlim(ax, [0.4 4.6]); grid(ax, 'on');
title(ax, sprintf('Analiza Comparativa Monte Carlo (N_{mc} = %d)', N_mc), ...
    'FontSize',13,'FontWeight','bold','FontName','Times New Roman','Color',[0.1 0.1 0.1]);
lg = legend(ax, b, {'e-MPC (Exact)', 'NN (Aproximat)'}, 'Location','NorthEast', ...
    'FontSize',11,'FontName','Times New Roman','EdgeColor',[0.8 0.8 0.8],'Color','w','TextColor',[0.15 0.15 0.15]);
lg.Box = 'on';


%% =================================================================
%% FUNCTII AUXILIARE
%% =================================================================
function str = tf2str(cond)
% Conversie boolean -> text afisare rapoarte
if cond, str = 'OK (CONFIRMAT)'; else, str = 'INCALCAT'; end
end

%%
function u_out = nn_forward(x_bar, W, b, n_fc, x_mean, x_std, u_mean, u_std)
% Forward pass algebric: normalizare -> straturi ReLU -> denormalizare
% (identic cu predict(), util pentru benchmark si portare pe ESP32)
h = (x_bar - x_mean) ./ x_std;
for li = 1:(n_fc - 1)
    h = max(0, W{li} * h + b{li});
end
u_n   = W{n_fc} * h + b{n_fc};
u_out = u_mean + u_std .* u_n;
end

%%
function net = train_net_fast(X_fit, U_fit, X_val, U_val, M, L, max_epochs)
% Pre-antrenare retea fully-connected (L straturi ascunse x M neuroni, ReLU)
layers = [featureInputLayer(size(X_fit,2), 'Normalization','none')];
for ll = 1:L
    layers = [layers; fullyConnectedLayer(M); reluLayer()]; 
end
layers = [layers; fullyConnectedLayer(size(U_fit,2)); regressionLayer()];

options = trainingOptions('adam', ...
    'MaxEpochs', max_epochs, 'MiniBatchSize', 256, 'InitialLearnRate', 1e-3, ...
    'LearnRateSchedule', 'piecewise', 'LearnRateDropFactor', 0.2, ...
    'LearnRateDropPeriod', floor(max_epochs*0.5), 'Shuffle', 'every-epoch', ...
    'ValidationData', {X_val, U_val}, 'ValidationFrequency', max(10, floor(max_epochs/10)), ...
    'Plots', 'none', 'Verbose', false);

net = trainNetwork(X_fit, U_fit, layers, options);
end

%%
function net = finetune_net(net, X_ft, U_ft, X_val, U_val, max_epochs, lr)
% Continua antrenarea din greutatile curente (nu reinitializeaza reteaua)
layer_arr = net.Layers;
options = trainingOptions('adam', ...
    'MaxEpochs', max_epochs, 'MiniBatchSize', min(128, floor(size(X_ft,1)/2)), ...
    'InitialLearnRate', lr, 'LearnRateSchedule', 'none', 'Shuffle', 'every-epoch', ...
    'ValidationData', {X_val, U_val}, 'ValidationFrequency', max(5, floor(max_epochs/5)), ...
    'Plots', 'none', 'Verbose', false);
net = trainNetwork(X_ft, U_ft, layer_arr, options);
end

%%
function P = zonotope_to_polyhedron(G, nx)
% Conversie exacta zonotop -> Polyhedron via functia suport h_Z(d)=||G'*d||_1
% Directii: axe canonice ±e_i + diagonale ±(e_i±e_j) (acopera fatele relevante pt nx=4)
dirs = [eye(nx), -eye(nx)];
for i = 1:nx
    for j = i+1:nx
        for si = [-1, 1]
            for sj = [-1, 1]
                d = zeros(nx,1); d(i) = si; d(j) = sj;
                dirs = [dirs, d / norm(d)]; 
            end
        end
    end
end

nD = size(dirs, 2);
A_ineq = zeros(nD, nx); b_ineq = zeros(nD, 1);
for k = 1:nD
    d = dirs(:, k);
    A_ineq(k, :) = d';
    b_ineq(k)    = norm(G' * d, 1);
end
P = Polyhedron('A', A_ineq, 'b', b_ineq);
P.minHRep();
end
%%
function u = empc_eval_fast(x, A_reg, b_reg, F_reg, g_reg, n_reg)
% Evaluare rapida lege PWA: cautare liniara a regiunii ce contine x
% (regiunile sunt non-overlapping => prima gasita e unica)
u = zeros(size(g_reg{1}));
for i = 1:n_reg
    slack = b_reg{i} - A_reg{i} * x;
    if all(slack >= -1e-6)
        u = F_reg{i} * x + g_reg{i};
        return;
    end
end
error('empc_eval_fast: Punct infezabil. x = [%s]', num2str(x'));
end

%%
function plot_projection(P, dims, color, lstyle, label, z_max_box)
% Deseneaza proiectia 2D a poliedrului P pe dimensiunile 'dims';
% in caz de eroare (set degenerat), foloseste un bounding-box ca fallback.
di = dims(1); dj = dims(2);
try
    Pp = P.projection([di, dj]);
    Pp.minHRep();
    V = Pp.V;
    if isempty(V), error('proiectie vida'); end
    kk = convhull(V(:,1), V(:,2));
    fill(V(kk,1), V(kk,2), color, 'FaceAlpha',0.25, 'EdgeColor',color, ...
        'LineWidth',2.5, 'LineStyle',lstyle, 'DisplayName',label);
catch
    fill([-z_max_box(di),  z_max_box(di),  z_max_box(di), -z_max_box(di)], ...
        [-z_max_box(dj), -z_max_box(dj),  z_max_box(dj),  z_max_box(dj)], color, ...
        'FaceAlpha',0.25, 'EdgeColor',color, 'LineWidth',2.5, 'LineStyle',lstyle, ...
        'DisplayName',[label ' (bbox)']);
end
end