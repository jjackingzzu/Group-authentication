clc; clear all; tic;

% System parameters
K_0 = [2 3 4 5];

RF_chain_cell = {
    [4 8 12 16 24 32], ...
    [64 100 128 128+64 256], ...
    [128 128+64 256 256+128 512], ...
    [512:(512)/4:1024] ...
};

theta_opt_deg = [45 90 45 76];   % Replace with the optimized angles in the paper
theta_opt = theta_opt_deg/180*pi;


frequency_central = 3.5;                   % GHz
noise_figure = 9;                          % dB
Band_width = 1e6;                         % Hz

Noise_thermal = 10^(-174/10)/1000;         % W/Hz
noise_sigma_square = Noise_thermal * Band_width * 10^(noise_figure/10);

distance_range = [400 600];
time_slots = 5;
message_length = 160;
block = message_length/(time_slots-1);
run_time = [1e5*4 1e4 1e4 1e4];
%run_time = [1e4 1e3 1e3 1e3];

% Authentication threshold from false alarm probability constraint
Pf_target = 1e-3;
L = message_length;

pf_cdf = zeros(1,L+1);
pmf = 2^(-L);
pf_cdf(1) = pmf;

for e = 1:L
    pmf = pmf * (L-e+1) / e;
    pf_cdf(e+1) = pf_cdf(e) + pmf;
end

eta = find(pf_cdf <= Pf_target,1,'last') - 1;
Pf_actual = pf_cdf(eta+1);

fprintf('Pf target = %.3e, eta = %d bits, Pf actual = %.3e\n', ...
        Pf_target, eta, Pf_actual);

% Storage
AFP = cell(1,length(K_0));
AFP_auth_only = cell(1,length(K_0));
SER_message = cell(1,length(K_0));
SER_tag = cell(1,length(K_0));

total_cases = 0;
for index_K = 1:length(K_0)
    total_cases = total_cases + length(RF_chain_cell{index_K});
end
case_count = 0;

for index_K = 1:length(K_0)

    K = K_0(index_K);
    theta = theta_opt(index_K);
    run_times=run_time(K-1);

    Constellation_group_selection;

    Constellation_sum = Constellation_table(:,K+1);
    Constellation_table_rotated = Constellation_table * exp(-1i*theta);
    Constellation_sum_rotated = Constellation_table_rotated(:,K+1);

    distance_link = linspace(distance_range(1), distance_range(2), K);

    P_L = 35.3*log10(distance_link) + 22.4 + 21.3*log10(frequency_central);
    P_L = 10.^(P_L/10);
    D = diag(P_L);
    power_transmit = 10^(30/10)/1000 / mean(P_L); % actual P_T=30dBm
    % The path loss P_L and the sqrt(D) pre-compensation in the encoding process cancel each other in the received signal.

    M_list = RF_chain_cell{index_K};

    AFP{index_K} = zeros(1,length(M_list));
    AFP_auth_only{index_K} = zeros(1,length(M_list));
    SER_message{index_K} = zeros(1,length(M_list));
    SER_tag{index_K} = zeros(1,length(M_list));

    for index_RF = 1:length(M_list)

        RF_chain = M_list(index_RF);

        error_message_temp = zeros(1,run_times);
        error_tag_temp = zeros(1,run_times);

        parfor times = 1:run_times

            data_bit = round(rand(K,time_slots-1,block));

            data_1 = data_bit .* Constellation_table(end,1:K).' ...
                   + (1-data_bit) .* Constellation_table(1,1:K).';

            data_2 = data_bit .* Constellation_table_rotated(end,1:K).' ...
                   + (1-data_bit) .* Constellation_table_rotated(1,1:K).';

            tag_bit = round(rand(1,time_slots-1,block));

            data_rotated = data_1 .* tag_bit + data_2 .* (1-tag_bit);

            data_sent = reshape(sum(data_rotated,1),1,message_length);

            data = [ones(K,1,block), data_rotated];

            H = randn(RF_chain,K,block)/sqrt(2) ...
              + 1i*randn(RF_chain,K,block)/sqrt(2);

            Y = sqrt(power_transmit) * pagemtimes(H,data) ...
              + (randn(RF_chain,time_slots,block) ...
              + 1i*randn(RF_chain,time_slots,block))/sqrt(2) ...
              * sqrt(noise_sigma_square);

            Y_tilde = Y / sqrt(power_transmit);

            Y_tilde1 = zeros(RF_chain,2*(time_slots-1),block);
            Y_tilde1(:,1:2:2*time_slots-3,:) = repmat(Y_tilde(:,1,:),1,time_slots-1,1);
            Y_tilde1(:,2:2:2*time_slots-2,:) = Y_tilde(:,2:time_slots,:);
            Y_tilde1 = reshape(Y_tilde1,RF_chain,2,message_length);

            R_12 = reshape(sum(conj(Y_tilde1(:,1,:)) .* Y_tilde1(:,2,:),1) ...
                   / RF_chain, 1, message_length);

            d_1 = abs(R_12 - Constellation_sum);
            [d_min_1,index_1] = min(d_1,[],1);
            data_estimated_1 = Constellation_sum(index_1).';

            d_2 = abs(R_12 - Constellation_sum_rotated);
            [d_min_2,index_2] = min(d_2,[],1);
            data_estimated_2 = Constellation_sum_rotated(index_2).';

            data_estimated = reshape(data_estimated_1,1,time_slots-1,block) .* tag_bit ...
                           + reshape(data_estimated_2,1,time_slots-1,block) .* (1-tag_bit);

            error_message_temp(times) = nnz(round(reshape(data_estimated,size(data_sent)) - data_sent,12));

            tag_bit_vec = reshape(tag_bit,1,message_length);
            tag_bit_hat_vec = d_min_1 < d_min_2;

            error_tag_temp(times) = nnz(tag_bit_hat_vec ~= tag_bit_vec);

        end

        AFP{index_K}(index_RF) = mean((error_message_temp > 0) | (error_tag_temp > eta));


        SER_message{index_K}(index_RF) = sum(error_message_temp) / (run_times * message_length);
        SER_tag{index_K}(index_RF) = sum(error_tag_temp) / (run_times * message_length);

        case_count = case_count + 1;
        progress_percent = case_count / total_cases * 100;

        fprintf('Progress: %.2f%% completed, K = %d, M = %d, TFP = %.3e\n', ...
                progress_percent, K, RF_chain, AFP{index_K}(index_RF));

    end
end

% Plot TFP versus M
figure;
set(gcf,'unit','normalized','position',[0.1,0.1,0.2*1.2*2,0.3*1.2*2]);

panel_label = {'(a)','(b)','(c)','(d)'};

for index_K = 1:length(K_0)

    subplot(2,2,index_K);

    x_data = RF_chain_cell{index_K};
    y_data = AFP{index_K};

    semilogy(x_data, y_data, '-or', ...
             'LineWidth',1,'MarkerSize',5);

    grid on; box on;
    set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);

    xlim([min(x_data), max(x_data)]);

    y_min = min(y_data(y_data > 0));
    if isempty(y_min)
        y_min = 1e-5;
    end
    ylim([y_min, 1]);

    xlabel('$M$','FontSize',14,'FontName','Times New Roman','Interpreter','latex');
    ylabel('AFP','FontSize',14,'FontName','Times New Roman','Interpreter','tex');

    title(['$K=',num2str(K_0(index_K)),'$'], ...
          'FontSize',14,'FontName','Times New Roman','Interpreter','latex');

    text(0.5, -0.25, panel_label{index_K}, ...
         'Units','normalized', ...
         'HorizontalAlignment','center', ...
         'FontSize',14, ...
         'FontName','Times New Roman');

end

fprintf('Simulation time = %.6f s\n', toc);