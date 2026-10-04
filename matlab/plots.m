
%%
%%grafice si iar grafice
%% ================================================================
%% HARTA REGIUNILOR eMPC 3D — VIZUALIZARE COMPLEXĂ (CRISTAL)
%% ================================================================
fprintf('\n[GRAFIC] Generare structura 3D pentru partitia eMPC...\n');

figure('Name', 'Partitia e-MPC 3D', 'Color', 'w', 'Position', [100, 100, 950, 700]);

% Taiem DOAR dimensiunea 4 (viteza unghiulara = 0)
% Raman libere: x_pos (Dim 1), v_pos (Dim 2), theta (Dim 3)
partition_3d = empc.partition.slice(4, 0);

% Plotam in 3D cu transparenta (opacity) pentru a vedea poliedrele din interior
% Folosim un contur foarte subtire (linewidth 0.1) ca sa se vada densitatea
partition_3d.plot('colormap', 'hsv', 'opacity', 0.4, 'edgecolor', [0.2 0.2 0.2], 'linewidth', 0.1);

grid on;
ax = gca;
ax.FontSize = 11;
ax.TickLabelInterpreter = 'latex';

% Etichete pentru cele 3 axe ramase
xlabel('Pozitia Caruciorului $x_{pos}$ (m)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Viteza Caruciorului $\dot{x}_{pos}$ (m/s)', 'Interpreter', 'latex', 'FontSize', 12);
zlabel('Unghiul Pendulului $\theta$ (rad)', 'Interpreter', 'latex', 'FontSize', 12);

title(sprintf('Partitia Multidimensionala e-MPC (Total Regiuni: %d) $\\dot{\\theta} = 0$', empc.nr), ...
    'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');

% Rotim automat camera pentru un unghi de impact vizual optim
view(35, 25); 


%%
%% Arhitectura retea
%%

%% ================================================================
%% GENERARE SCHEMĂ INGENEREASCĂ ARHITECTURĂ REȚEA NEURONALĂ DENSĂ
%% Arhitectura: 4 -> [30 x 8] -> 1 (ReLU)
%% ================================================================


% Setări de bază pentru rețea
nr_intrari = 4;
nr_straturi_ascunse = 8;
neuroni_ascunsi = 30;
nr_iesiri = 1;

% Configurație de desenare stilizată (câți neuroni afișăm efectiv pe ecran)
nodes_input = 4;
nodes_hidden = 5; % Se vor afișa 3 sus, puncte, 2 jos pentru a simula 30
nodes_output = 1;

% Coordonatele straturilor pe axa X (Total = 1 + 8 + 1 = 10 straturi)
X_coords = 1 : (1 + nr_straturi_ascunse + 1);
N_layers = length(X_coords);

% Inițializare figură profesională
figure('Name', 'Arhitectura Retelei Neuronale eMPC', 'Color', 'w', 'Position', [100, 100, 1100, 550]);
hold on;

% Culori conform paletei academice (Stil IEEE)
color_input  = [0.00, 0.45, 0.74];  % Albastru
color_hidden = [0.47, 0.67, 0.19];  % Verde (ReLU)
color_output = [0.85, 0.33, 0.10];  % Portocaliu/Roșu
color_line   = [0.82, 0.82, 0.85];  % Gri deschis pentru conexiuni (evită aglomerarea vizuală)

% Celulă în care salvăm coordonatele Y pentru fiecare strat desenat
Y_positions = cell(1, N_layers);

%% 1. Calcularea pozițiilor Y pentru fiecare strat (Centrate în jurul lui 0)
for i = 1:N_layers
    if i == 1
        % Stratul de intrare (4 noduri)
        Y_positions{i} = linspace(1.5, -1.5, nodes_input);
    elseif i == N_layers
        % Stratul de ieșire (1 nod)
        Y_positions{i} = 0;
    else
        % Straturile ascunse (Afișăm indicii 1, 2, 3 și 29, 30 distanțați)
        Y_positions{i} = [1.6, 1.0, 0.4, -0.8, -1.4];
    end
end

%% 2. Desenarea conexiunilor (Synapses) dintre straturi
fprintf('[INFO] Randare conexiuni sinaptice dense...\n');
for i = 1 : (N_layers - 1)
    x1 = X_coords(i);
    x2 = X_coords(i+1);
    y1_set = Y_positions{i};
    y2_set = Y_positions{i+1};
    
    % Conectăm fiecare nod vizibil din stratul curent cu cel din stratul următor
    for j = 1:length(y1_set)
        for k = 1:length(y2_set)
            % Folosim o linie fină, semi-transparentă pentru un efect elegant
            line([x1, x2], [y1_set(j), y2_set(k)], 'Color', color_line, 'LineWidth', 0.6);
        end
    end
end

%% 3. Desenarea nodurilor (Neuroni) și adnotărilor speciale (\vdots)
fprintf('[INFO] Randare noduri neuronale...\n');
for i = 1:N_layers
    x = X_coords(i);
    y_set = Y_positions{i};
    
    if i == 1
        % Plotare noduri intrare
        scatter(ones(size(y_set))*x, y_set, 160, color_input, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
        % Etichete intrări
        labels_in = {'$x_{pos}$', '$\dot{x}_{pos}$', '$\theta$', '$\dot{\theta}$'};
        for j = 1:length(y_set)
            text(x - 0.3, y_set(j), labels_in{j}, 'Interpreter', 'latex', 'FontSize', 11, 'HorizontalAlignment', 'right');
        end
        
    elseif i == N_layers
        % Plotare nod ieșire
        scatter(x, y_set, 180, color_output, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
        text(x + 0.15, y_set, '$u_{NN}$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
        
    else
        % Plotare noduri straturi ascunse
        scatter(ones(size(y_set))*x, y_set, 110, color_hidden, 'filled', 'MarkerEdgeColor', [0.2 0.3 0.1]);
        
        % Introducere puncte de suspensie (\vdots) în mijlocul stratului ascuns
        text(x, -0.1, '$\vdots$', 'Interpreter', 'latex', 'FontSize', 18, ...
            'HorizontalAlignment', 'center', 'Color', [0.3 0.3 0.3]);
        
        % Numerotare simbolică a neuronilor din strat (opțional, pentru rigoare)
        text(x, y_set(1)+0.22, '$n_1$', 'Interpreter', 'latex', 'FontSize', 8, 'HorizontalAlignment', 'center', 'Color', [0.4 0.4 0.4]);
        text(x, y_set(3)+0.22, '$n_3$', 'Interpreter', 'latex', 'FontSize', 8, 'HorizontalAlignment', 'center', 'Color', [0.4 0.4 0.4]);
        text(x, y_set(5)-0.22, '$n_{30}$', 'Interpreter', 'latex', 'FontSize', 8, 'HorizontalAlignment', 'center', 'Color', [0.4 0.4 0.4]);
    end
end

%% 4. Adăugare etichete superioare și inferioare pentru Straturi
% Delimitare zone prin text structural
text(X_coords(1), 2.2, 'Input Layer', 'Interpreter', 'latex', 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(mean(X_coords(2:end-1)), 2.4, 'Hidden Layer (\textit{Hidden Layers}) --- Type: \textit{Fully Connected}', 'Interpreter', 'latex', 'FontSize', 11, 'Color', color_hidden*0.8, 'HorizontalAlignment', 'center');
text(X_coords(end), 2.2, 'Output Layer', 'Interpreter', 'latex', 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Numerotarea fiecărui strat ascuns în parte la baza graficului
for h = 1:nr_straturi_ascunse
    text(X_coords(1+h), 2.1, sprintf('$H_{%d}$', h), 'Interpreter', 'latex', 'FontSize', 10, 'HorizontalAlignment', 'center');
    text(X_coords(1+h), -2.1, '30 $\mathcal{N}$', 'Interpreter', 'latex', 'FontSize', 9, 'HorizontalAlignment', 'center', 'Color', [0.3 0.3 0.3]);
    text(X_coords(1+h), -2.3, 'ReLU', 'Interpreter', 'latex', 'FontSize', 8, 'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5]);
end

text(X_coords(1), -2.1, '4 Inputs', 'Interpreter', 'none', 'FontSize', 9, 'HorizontalAlignment', 'center');
text(X_coords(end), -2.1, '1 Output (liniar)', 'Interpreter', 'latex', 'FontSize', 9, 'HorizontalAlignment', 'center');

%% 5. Finisare aspect estetic (Tuning axe)
title('Neural Network Architecture', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
axis off; % Ascundem axele carteziene standard pentru aspect de diagramă curată
xlim([0.5, N_layers + 0.8]);
ylim([-2.6, 2.6]);

fprintf('[GATA] Diagrama a fost generata cu succes.\n');


%% --- GENERARE GRAFIC LOSS (MSE) REAL ---
%% plot_loss.m - Loss curve, rafinat
% Inlocuieste mse_train, mse_val, iter_train, iter_val cu datele tale


%% Culori din imagine
BLUE   = [0.20 0.45 0.75];   % albastru - Z (mRPI)
PINK   = [0.88 0.35 0.65];   % roz/mov  - Z_NN

%% Date sintetice (inlocuieste cu datele tale reale)
rng(42);
N          = 95000;
iter_train = (1:N)';
base       = 0.5 * exp(-iter_train/3000);
mse_train  = base .* exp(randn(N,1)*1.1) * 0.3;
mse_train(iter_train>50000) = mse_train(iter_train>50000) * 0.25;
mse_train  = max(mse_train, 7e-5 + abs(randn(N,1))*1e-5);

iter_val = (100:100:N)';
mse_val  = movmean(mse_train,500); 
mse_val  = mse_val(iter_val)*2 + abs(randn(length(iter_val),1))*2e-5;
mse_val  = max(mse_val, 9e-5);

%% Medie mobila
W     = 800;
sm    = movmean(mse_train, W);

%% Plot
figure('Color','white','Units','centimeters','Position',[2 2 22 11]);
ax = axes('YScale','log','Box','on','FontName','Times New Roman',...
          'FontSize',11,'XGrid','on','YGrid','on','YMinorGrid','on',...
          'GridColor',[0.85 0.85 0.85],'TickDir','out');
hold on;

% banda incertitudine
win = 1000; nw = floor(N/win);
xm=zeros(nw,1); pl=xm; ph=xm;
for i=1:nw
    idx=(i-1)*win+1:i*win;
    pl(i)=prctile(mse_train(idx),10);
    ph(i)=prctile(mse_train(idx),90);
    xm(i)=mean(iter_train(idx));
end
fill([xm;flipud(xm)],[pl;flipud(ph)],BLUE,...
     'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');

% curba bruta
plot(iter_train, mse_train,'Color',[BLUE 0.15],'LineWidth',0.4,...
     'HandleVisibility','off');

% medie mobila antrenare
h1 = plot(iter_train, sm, 'Color',BLUE,'LineWidth',2,...
          'DisplayName','MSE_{train} (avg)');

% validare
h2 = plot(iter_val, mse_val,'Color',PINK,'LineWidth',1.6,...
          'Marker','o','MarkerSize',3,'MarkerFaceColor',PINK,...
          'MarkerEdgeColor','none','DisplayName','MSE_{val}');

% linie LR drop
xline(50000,'--','Color',[0.4 0.4 0.4],'LineWidth',1.1,...
      'HandleVisibility','off');
text(50500, 5e-1,'LR \downarrow','FontName','Times New Roman',...
     'FontSize',9,'Interpreter','tex','Color',[0.4 0.4 0.4]);

%% Formatare
xlabel('Iteratii','FontName','Times New Roman','FontSize',12);
ylabel('MSE','FontName','Times New Roman','FontSize',12);
title('Evolutia functiei de Loss','FontName','Times New Roman',...
      'FontSize',13,'FontWeight','bold');
legend([h1 h2],'Location','northeast','FontSize',10,...
       'FontName','Times New Roman','Box','on');
xlim([0 N]); ylim([5e-6 2]);

%% Export
exportgraphics(gcf,'loss_curve.png','Resolution',300,'BackgroundColor','white');
exportgraphics(gcf,'loss_curve.pdf','ContentType','vector','BackgroundColor','white');
fprintf('Salvat: loss_curve.png si loss_curve.pdf\n');  

%%
%%

%% ================================================================
%% GRAFIC DETALIU: Eroarea instantanee de aproximare online e_u(k)
%% FĂRĂ clc / clear / close all !
%% ================================================================

e_u = abs(epsilon_history(1, :));
t_sim = (0:N_sim-1) * Ts;

eps_max_retrain = 0.1582;
eps_max_online  = 0.0306;

figure('Name', 'Approximation Error e_u(k)', 'Color', 'w', 'Position', [250, 250, 900, 420]);

fill([0, t_sim(end), t_sim(end), 0], ...
     [0, 0, eps_max_retrain, eps_max_retrain], ...
     [0.47, 0.67, 0.19], 'FaceAlpha', 0.07, 'EdgeColor', 'none', 'HandleVisibility', 'off');
hold on;

plot(t_sim, e_u, '-o', 'Color', [0.47, 0.67, 0.19], 'LineWidth', 1.5, ...
    'MarkerSize', 4, 'MarkerFaceColor', [0.47, 0.67, 0.19], ...
    'DisplayName', '$|e_u(k)|$');

yline(eps_max_retrain, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.8, ...
    'DisplayName', sprintf('$\\epsilon_{\\max}^{\\mathrm{train}} = %.4f$', eps_max_retrain));

yline(eps_max_online, ':', 'Color', [0.0 0.45 0.70], 'LineWidth', 1.6, ...
    'DisplayName', sprintf('$\\epsilon_{\\max}^{\\mathrm{sim}} = %.4f$', eps_max_online));

grid on;
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.7;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 11;

xlabel('Time (s)',         'Interpreter', 'latex', 'FontSize', 12);
ylabel('$|e_u(k)|$ (N)',  'Interpreter', 'latex', 'FontSize', 12);
title(['\textbf{Online approximation error }', ...
       '$e_u(k) = |u_{\mathrm{eMPC}}(\bar{x}_k) - u_{\mathrm{NN}}(\bar{x}_k)|$'], ...
      'Interpreter', 'latex', 'FontSize', 13);

text(t_sim(end)*0.02, eps_max_retrain * 0.78, ...
    'Region guaranteed by robust tube $\mathcal{Z}_{NN}$', ...
    'Interpreter', 'latex', 'FontSize', 10, 'Color', [0.2 0.5 0.1]);

legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 11);
xlim([0, t_sim(end)]);
ylim([0, eps_max_retrain * 1.35]);
%%
%%

nx = 4; nu = 1;
alpha_bytes = 4;

configs = {
    'N(6,12)',   6,  12;
    'N(8,30)',   8,  30;
    'N(10,50)',  10, 50;
    'N(12,64)',  12, 64;
    'N(16,100)', 16, 100;
    'N(20,128)', 20, 128;
};
n = size(configs, 1);
mem_KB  = zeros(n, 1);
n_par   = zeros(n, 1);
labels  = cell(n, 1);

for i = 1:n
    L = configs{i,2}; M = configs{i,3};
    p = (nx+1)*M + (L-1)*(M+1)*M + (M+1)*nu;
    n_par(i)  = p;
    mem_KB(i) = alpha_bytes * p / 1024;
    labels{i} = configs{i,1};
end

ratio_vs_min = mem_KB / mem_KB(1);

% --- Culori: gradient albastru -> rosu in functie de marime ---
cmap = [

    0.00 0.45 0.74;   % N(6,12)

    0.00 0.45 0.74;   % N(8,30)

    0.00 0.45 0.74;   % N(10,50)

    0.00 0.45 0.74;   % N(12,64)

    1.00 0.70 0.70;   % N(16,100) - rosu pal

    1.00 0.60 0.60;    % N(20,128) - rosu aprins

];

figure('Color','w','Position',[80 80 1000 520], ...
    'Name','Amprenta de memorie NN');

ax = axes('Color','w','Box','off','FontSize',10.5, ...
    'GridColor',[0.88 0.88 0.88],'GridLineStyle',':', ...
    'GridAlpha',1,'LineWidth',0.8,'TickDir','out', ...
    'FontName','Times New Roman');
hold(ax,'on');

% --- Bare ---
for i = 1:n
    b = bar(ax, i, mem_KB(i), 'BarWidth',0.62, ...
        'FaceColor',cmap(i,:),'EdgeColor','none');
end

% --- Linie ESP32 SRAM (512 KB) ---
yline(ax, 512, '--', 'Color',[0.55 0.55 0.55], 'LineWidth',1.4, ...
    'Label','  ESP32 SRAM 512 KB', ...
    'LabelHorizontalAlignment','left', ...
    'FontSize',9,'FontName','Times New Roman', ...
    'LabelVerticalAlignment','bottom');

% --- Linie 1 MB ---
yline(ax, 1024, ':', 'Color',[0.80 0.20 0.20], 'LineWidth',1.6, ...
    'Label','  1 MB', ...
    'LabelHorizontalAlignment','left', ...
    'FontSize',9,'FontName','Times New Roman', ...
    'LabelVerticalAlignment','bottom');

% --- Adnotari pe bare ---
offset = max(mem_KB) * 0.025;
for i = 1:n
    % Memorie formatata
    if mem_KB(i) >= 1024
        mem_str = sprintf('%.2f MB', mem_KB(i)/1024);
    else
        mem_str = sprintf('%.1f KB', mem_KB(i));
    end

    % Linie 1: memorie
    text(ax, i, mem_KB(i) + offset, mem_str, ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',9,'FontWeight','bold','Color',[0.10 0.10 0.10], ...
        'FontName','Times New Roman');

    % Linie 2: numar parametri
    text(ax, i, mem_KB(i) + offset*3.8, sprintf('%d par.', n_par(i)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',8,'Color',[0.35 0.35 0.35], ...
        'FontName','Times New Roman');

    % Linie 3: raport fata de cea mai mica arhitectura (doar daca > 1x)
    if ratio_vs_min(i) > 1.01
        text(ax, i, -max(mem_KB)*0.085, sprintf('×%.1f', ratio_vs_min(i)), ...
            'HorizontalAlignment','center','VerticalAlignment','top', ...
            'FontSize',8.5,'Color',cmap(i,:)*0.75, ...
            'FontWeight','bold','FontName','Times New Roman');
    else
        text(ax, i, -max(mem_KB)*0.085, '×1', ...
            'HorizontalAlignment','center','VerticalAlignment','top', ...
            'FontSize',8.5,'Color',[0.45 0.45 0.45], ...
            'FontName','Times New Roman');
    end
end
 v
% Eticheta raport sub axa
text(ax, 0.38, -max(mem_KB)*0.085, 'raport:', ...
    'HorizontalAlignment','right','VerticalAlignment','top', ...
    'FontSize',8,'Color',[0.50 0.50 0.50],'FontAngle','italic', ...
    'FontName','Times New Roman');

% --- Zona "sigura" ESP32 (umbrire) ---
patch(ax, [0.4 n+0.6 n+0.6 0.4], [0 0 512 512], ...
    [0.60 0.85 0.65], 'FaceAlpha',0.07,'EdgeColor','none');
text(ax, n+0.55, 256, 'zona sigura ESP32', ...
    'HorizontalAlignment','right','VerticalAlignment','middle', ...
    'FontSize',8,'Color',[0.25 0.55 0.35],'FontAngle','italic', ...
    'FontName','Times New Roman');

% --- Formatare axe ---
xticks(ax, 1:n);
xticklabels(ax, labels);
set(ax,'XTickLabelRotation',12,'FontName','Times New Roman');
ylabel(ax, 'Memorie \Gamma_N  [KB]', 'FontSize',12, ... §
    'FontName','Times New Roman');
ylim(ax, [-max(mem_KB)*0.13, max(mem_KB)*1.32]);
xlim(ax, [0.35, n+0.65]);
grid(ax, 'on');

title(ax, 'Evolutia memoriei la diverse arhitecturi', ...
    'FontSize',13,'FontWeight','bold','FontName','Times New Roman', ...
    'Color',[0.10 0.10 0.10]);


%% Figura 6 MODIFICATA: Seturi mRPI Z vs Z_NN + verificare invarianta

%% --- Construim zonotopul A_K*Z_NN (+) D ---
G_AKZNN_D = [A_K * G_ZNN, G_D];

dirs_tmp = [eye(nx), -eye(nx)];
for ii_ = 1:nx
    for jj_ = ii_+1:nx
        for si_ = [-1, 1]
            for sj_ = [-1, 1]
                d_ = zeros(nx,1);
                d_(ii_) = si_; d_(jj_) = sj_;
                dirs_tmp = [dirs_tmp, d_ / norm(d_)];
            end
        end
    end
end
nD_tmp     = size(dirs_tmp, 2);
A_ineq_tmp = zeros(nD_tmp, nx);
b_ineq_tmp = zeros(nD_tmp, 1);
for kk_ = 1:nD_tmp
    d_                = dirs_tmp(:, kk_);
    A_ineq_tmp(kk_,:) = d_';
    b_ineq_tmp(kk_)   = norm(G_AKZNN_D' * d_, 1);
end
Z_AKZNN_D = Polyhedron('A', A_ineq_tmp, 'b', b_ineq_tmp);
Z_AKZNN_D.minHRep();

%% --- Culori ---
col_znn = [0.72, 0.18, 0.72];   % mov       — Z_NN (bordura groasa)
col_mpc = [0.18, 0.44, 0.70];   % albastru  — Z
col_inv = [0.80, 0.80, 0.80]; % gri       — A_K*Z_NN (+) D

pairs   = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
n_pairs = size(pairs, 1);
ax_lbls = {'$x_1$ pos (m)', '$x_2$ vel (m/s)', ...
           '$x_3$ ang (rad)', '$x_4$ $\dot{\theta}$ (rad/s)'};

figure('Name','Seturi mRPI: Z vs Z_NN + invarianta', ...
    'Color','w','Position',[100, 50, 1150, 950]);

for pp = 1:n_pairs
    di = pairs(pp,1);
    dj = pairs(pp,2);
    subplot(3,2,pp);
    hold on;

    %% 1) Z_NN — fundal transparent, bordura mov FOARTE GROASA
    try
        Zp_nn = Z_NN.projection([di, dj]); Zp_nn.minHRep();
        V_ZNN = Zp_nn.V;
        if ~isempty(V_ZNN)
            k_nn = convhull(V_ZNN(:,1), V_ZNN(:,2));
            fill(V_ZNN(k_nn,1), V_ZNN(k_nn,2), col_znn, ...
                'FaceAlpha', 0.08, 'EdgeColor', col_znn, ...
                'LineWidth', 5.0,  'LineStyle', '-', ...
                'DisplayName', '$\mathcal{Z}_{NN}$');
        end
    catch
        fill([-z_max_nn(di), z_max_nn(di), z_max_nn(di), -z_max_nn(di)], ...
             [-z_max_nn(dj),-z_max_nn(dj), z_max_nn(dj),  z_max_nn(dj)], col_znn, ...
            'FaceAlpha',0.08,'EdgeColor',col_znn,'LineWidth',5.0,'LineStyle','-', ...
            'DisplayName','$\mathcal{Z}_{NN}$ (bbox)');
    end

    %% 2) A_K*Z_NN (+) D — gri cu bordura inchisa, DEASUPRA fondului mov
    try
        Zp_inv = Z_AKZNN_D.projection([di, dj]); Zp_inv.minHRep();
        V_inv  = Zp_inv.V;
        if ~isempty(V_inv)
            k_inv = convhull(V_inv(:,1), V_inv(:,2));
            fill(V_inv(k_inv,1), V_inv(k_inv,2), col_inv, ...
                'FaceAlpha', 0.55, 'EdgeColor', [0.35, 0.35, 0.35], ...
                'LineWidth', 1.8,  'LineStyle', '--', ...
                'DisplayName', '$A_K\mathcal{Z}_{NN} \oplus \mathcal{D}$');
        end
    catch
        lhs = abs(A_K([di,dj],:)) * z_max_nn + [d_bound(di); d_bound(dj)];
        fill([-lhs(1), lhs(1), lhs(1), -lhs(1)], ...
             [-lhs(2),-lhs(2), lhs(2),  lhs(2)], col_inv, ...
            'FaceAlpha',0.55,'EdgeColor',[0.35,0.35,0.35],'LineWidth',1.8,'LineStyle','--', ...
            'DisplayName','$A_K\mathcal{Z}_{NN} \oplus \mathcal{D}$ (bbox)');
    end

    %% 3) Z — albastru solid, cel mai sus
    try
        Zp = Z.projection([di, dj]); Zp.minHRep();
        V_Z = Zp.V;
        if ~isempty(V_Z)
            k_z = convhull(V_Z(:,1), V_Z(:,2));
            fill(V_Z(k_z,1), V_Z(k_z,2), col_mpc, ...
                'FaceAlpha', 0.50, 'EdgeColor', col_mpc, ...
                'LineWidth', 2.5,  'LineStyle', '-', ...
                'DisplayName', '$\mathcal{Z}$ (mRPI)');
        end
    catch
        fill([-z_max(di), z_max(di), z_max(di), -z_max(di)], ...
             [-z_max(dj),-z_max(dj), z_max(dj),  z_max(dj)], col_mpc, ...
            'FaceAlpha',0.50,'EdgeColor',col_mpc,'LineWidth',2.5,'LineStyle','-', ...
            'DisplayName','$\mathcal{Z}$ (bbox)');
    end

    %% 4) Origine
    plot(0, 0, 'k+', 'MarkerSize', 9, 'LineWidth', 2.0, 'HandleVisibility', 'off');

    grid on; axis equal;
    ax = gca;
    ax.GridLineStyle = ':';
    ax.GridAlpha     = 0.6;
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = 10;
    xlabel(ax_lbls{di}, 'Interpreter','latex','FontSize',10);
    ylabel(ax_lbls{dj}, 'Interpreter','latex','FontSize',10);
    title(['Proiectie $(x_{' num2str(di) '}, x_{' num2str(dj) '})$'], ...
        'Interpreter','latex','FontSize',11);
    legend('Location','best','FontSize',8,'Interpreter','latex');
end
sgtitle(['$\mathcal{Z}$ (albastru) $\subseteq$ $\mathcal{Z}_{NN}$ (mov) ' ...
         '--- $A_K\mathcal{Z}_{NN} \oplus \mathcal{D}$ (gri) ' ...
         '$\subseteq \mathcal{Z}_{NN}$'], ...
    'FontSize',12,'FontWeight','bold','Interpreter','latex');



%%
%%

%% ================================================================
%% PASUL 14: Grafice — tubes vizibile + set mRPI Z vs Z_NN
%% ================================================================
t = (0:N_sim) * Ts;
col_mpc  = [0.00, 0.45, 0.74];
col_empc = [0.85, 0.33, 0.10];
col_nn   = [0.47, 0.67, 0.19];
col_znn  = [0.75, 0.00, 0.75];

% Culori vii pentru traiectorii reale — vizibile pe fundal albastru/mov
col_mpc_real  = [1.00, 1.00, 0.00];   % galben
col_empc_real = [1.00, 0.40, 0.00];   % portocaliu aprins
col_nn_real   = [0.00, 1.00, 0.80];   % cyan/turcoaz

state_labels_y = {'Pozitie Carucior (m)', 'Viteza Carucior (m/s)', ...
                  'Unghi Pendul (rad)',    'Viteza Unghiulara (rad/s)'};
state_titles   = {'$x_1$ -- Pozitie Carucior', '$x_2$ -- Viteza Carucior', ...
                  '$x_3$ -- Unghi Pendul',      '$x_4$ -- Viteza Unghiulara'};

figure('Name','Stari: Z vs Z_NN', 'Color','w', 'Position',[50,30,1100,1100]);

% for j = 1:nx
%     ax = subplot(4,1,j);
% end
%% ================================================================
%% PASUL 14: Grafice — tubes vizibile + set mRPI Z vs Z_NN
%% ================================================================
t = (0:N_sim) * Ts;
col_mpc  = [0.00, 0.45, 0.74];
col_empc = [0.85, 0.33, 0.10];
col_nn   = [0.47, 0.67, 0.19];
col_znn  = [0.75, 0.00, 0.75];

% Culori vii pentru traiectorii reale — vizibile pe fundal albastru/mov
col_mpc_real  = [1.00, 1.00, 0.00];   % galben
col_empc_real = [1.00, 0.40, 0.00];   % portocaliu aprins
col_nn_real   = [0.00, 1.00, 0.80];   % cyan/turcoaz

state_labels_y = {'Pozitie Carucior (m)', 'Viteza Carucior (m/s)', ...
                  'Unghi Pendul (rad)',    'Viteza Unghiulara (rad/s)'};
state_titles   = {'$x_1$ -- Pozitie Carucior', '$x_2$ -- Viteza Carucior', ...
                  '$x_3$ -- Unghi Pendul',      '$x_4$ -- Viteza Unghiulara'};

figure('Name','Stari: Z vs Z_NN', 'Color','w', 'Position',[50,30,1100,1100]);

for j = 1:nx
    ax = subplot(4,1,j);

    %% Tuburi
    % Z_NN — fundal mov
    tup_nn = x_bar_nn_history(j,:) + z_max_nn(j);
    tdn_nn = x_bar_nn_history(j,:) - z_max_nn(j);
     fill([t, fliplr(t)], [tup_nn, fliplr(tdn_nn)], col_znn, ...
        'FaceAlpha',0.15, 'EdgeColor',col_znn, 'EdgeAlpha',0.6, ...
        'LineWidth',1.0,  'LineStyle','-', ...
        'DisplayName','$Z_{NN}$ (tub NN)');

    hold on;

    % Z — albastru
    tup_z = x_bar_history(j,:) + z_max(j);
    tdn_z = x_bar_history(j,:) - z_max(j);
    fill([t, fliplr(t)], [tup_z, fliplr(tdn_z)], col_mpc, ...
        'FaceAlpha',0.18, 'EdgeColor',col_mpc, 'EdgeAlpha',0.7, ...
        'LineWidth',1.0,  'LineStyle','--', ...
        'DisplayName','$Z$ (tub MPC)');
    %% Traiectorii reale — culori vii, stiluri diferite, groase
    plot(t, x_history(j,:),      '--', 'Color',col_mpc_real,  'LineWidth',2.0, ...
        'DisplayName','$x$ MPC');
    plot(t, x_empc_history(j,:), ':',  'Color',col_empc_real, 'LineWidth',2.5, ...
        'DisplayName','$x$ eMPC');
    plot(t, x_nn_history(j,:),   '-.', 'Color',col_nn_real,   'LineWidth',2.0, ...
        'DisplayName','$x$ NN');

    %% Limite fizice
    yline(x_max(j), 'k:', 'LineWidth',1.8, 'HandleVisibility','off');
    yline(x_min(j), 'k:', 'LineWidth',1.8, 'DisplayName','Limite fizice');

    %% Aspect
    grid on; box on;
    set(ax, 'FontSize',11, 'TickLabelInterpreter','latex');
    ylabel(state_labels_y{j}, 'Interpreter','none', 'FontSize',11);
    title(state_titles{j}, 'Interpreter','latex', 'FontSize',13);

    if j < nx
        set(ax, 'XTickLabel', {});
    else
        xlabel('Timp (s)', 'Interpreter','latex', 'FontSize',12);
    end

    legend('Location','northeast', 'FontSize',9, 'Interpreter','latex', ...
           'NumColumns',2, 'Box','on');

    y_range = max(x_max(j) - x_min(j), 1e-3);
    ylim([x_min(j) - 0.08*y_range, x_max(j) + 0.08*y_range]);
end

sgtitle('Evolutia Starilor -- $Z$ (albastru) vs $Z_{NN}$ (mov)', ...
    'FontSize',14, 'FontWeight','bold', 'Interpreter','latex');