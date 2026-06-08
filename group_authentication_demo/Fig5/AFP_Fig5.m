clc;
clear;
close all;

%% ============================================================
% Load data
%% ============================================================

file_list = {
    'CHT-PLA/CHT-PLA.mat'
    'MCT-PLA/MCT-PLA.mat'
    'MTE-PLA/MTE-PLA.mat'
    'SIT-PLA/SIT-PLA.mat'
    'SLO-PLA/SLO-PLA.mat'
    'TDT-PLA/TDT-PLA.mat'
    };

legend_text = {
    'CHT-PLA'
    'MCT-PLA'
    'MTE-PLA'
    'SIT-PLA'
    'SLO-PLA'
    'TDT-PLA'
    };

num_file = length(file_list);

AFP = cell(num_file,4);

for file_index = 1:num_file

    data = load(file_list{file_index});

    AFP(file_index,:) = data.AFP;

    if file_index == 1

        RF_chain_cell = data.RF_chain_cell;
        K_0 = data.K_0;

    end

end

%% ============================================================
% Plot
%% ============================================================

figure;

set(gcf,...
    'unit','normalized',...
    'position',[0.1,0.1,0.2*1.2*2,0.3*1.2*2]);

panel_label = {'(a)','(b)','(c)','(d)'};

color_style = {...
    '-or', ...    % CHT
    '-sg', ...    % MCT
    '-dc', ...    % MTE
    '-k*', ...    % SIT
    '-sb', ...    % SLO
    '-m^'  ...    % TDT
    };

for index_K = 1:length(K_0)

    subplot(2,2,index_K);

    hold on;
    box on;
    grid on;

    x_data = RF_chain_cell{index_K};

    for scheme_index = 1:num_file

        y_data = AFP{scheme_index,index_K};

        loglog(x_data,...
               y_data,...
               color_style{scheme_index},...
               'LineWidth',1,...
               'MarkerSize',6);

    end

    %% Force log-log axes

    set(gca,...
        'XScale','log',...
        'YScale','log',...
        'FontName','Times New Roman',...
        'FontSize',14,...
        'LineWidth',1);

    %% Axis range

    xlim([min(x_data),max(x_data)]);

    y_all = [];

    for scheme_index = 1:num_file

        y_all = [y_all; AFP{scheme_index,index_K}(:)];

    end

    y_min = min(y_all(y_all>0));

    ylim([0.9*y_min,1]);

    if isempty(y_min)

        y_min = 1e-5;

    end

    %ylim([10^floor(log10(y_min)),1]);

    %% Labels

    xlabel('$M$',...
        'FontSize',14,...
        'FontName','Times New Roman',...
        'Interpreter','latex');

    ylabel('AFP',...
        'FontSize',14,...
        'FontName','Times New Roman');

    % title(['$K=',num2str(K_0(index_K)),'$'],...
    %     'FontSize',14,...
    %     'FontName','Times New Roman',...
    %     'Interpreter','latex');

    %% Legend

    legend(legend_text,...
        'Location','southwest',...
        'FontSize',12,...
        'FontName','Times New Roman');

    %% Panel label

    text(0.5,...
         -0.25,...
         panel_label{index_K},...
         'Units','normalized',...
         'HorizontalAlignment','center',...
         'FontSize',14,...
         'FontName','Times New Roman');

end
%exportgraphics(gcf,'AFP_Fig5.pdf','Resolution',600);