clc; clear all; tic;

% CoTag-PLA: multiple users cooperatively transmit one group tag
% and perform unified group authentication.

K_0 = [2 3 4 5];

RF_chain_cell = {
    [4 8 12 16 24 32], ...
    [64 100 128 128+64 256], ...
    [128 128+64 256 256+128 512], ...
    [512:(512)/4:1024] ...
};

rho_cell = {
    0.001:0.001:0.03, ...
    0.0002:0.0002:0.01, ...
    0.0002:0.0002:0.01, ...
    0.0002:0.0002:0.005 ...
};

theta = 45/180*pi;                         % Tag direction


frequency_central = 3.5;                   % GHz
noise_figure = 9;                          % dB
Band_width = 1e6;                         % Hz

Noise_thermal = 10^(-174/10)/1000;         % W/Hz
noise_sigma_square = Noise_thermal * Band_width * 10^(noise_figure/10);

distance_range = [400 600];
time_slots = 5;
message_length = 160;
block = message_length/(time_slots-1);
run_time = [1e5 1e4 1e4 1e4];

% Unified authentication threshold
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

AFP_rho = cell(1,length(K_0));
AFP = cell(1,length(K_0));
rho_opt = cell(1,length(K_0));

total_cases = 0;
for index_K = 1:length(K_0)
    total_cases = total_cases + length(RF_chain_cell{index_K}) * length(rho_cell{index_K});
end

case_count = 0;

for index_K = 1:length(K_0)

    K = K_0(index_K);
    run_times = run_time(index_K);
    rho_0 = rho_cell{index_K};

    Constellation_group_selection;

    Constellation_sum = Constellation_table(:,K+1);

    distance_link = linspace(distance_range(1), distance_range(2), K);

    P_L = 35.3*log10(distance_link) + 22.4 + 21.3*log10(frequency_central);
    P_L = 10.^(P_L/10);
    D = diag(P_L);
    power_transmit = 10^(30/10)/1000 / mean(P_L); % actual P_T=30dBm
    % The path loss P_L and the sqrt(D) pre-compensation cancel each other.

    M_list = RF_chain_cell{index_K};

    AFP{index_K} = zeros(1,length(M_list));
    AFP_rho{index_K} = zeros(length(rho_0),length(M_list));
    rho_opt{index_K} = zeros(1,length(M_list));

    for index_RF = 1:length(M_list)

        RF_chain = M_list(index_RF);

        for index_rho = 1:length(rho_0)

            rho = rho_0(index_rho);

            error_message_temp = zeros(1,run_times);
            error_tag_temp = zeros(1,run_times);

            parfor times = 1:run_times

                data_bit = round(rand(K,time_slots-1,block));

                data_symbol = data_bit .* Constellation_table(end,1:K).' ...
                            + (1-data_bit) .* Constellation_table(1,1:K).';

                data_sent = reshape(sum(data_symbol,1),1,message_length);

                tag_bit = round(rand(1,time_slots-1,block));
                tag_symbol = 2*tag_bit - 1;

                tx_symbol = sqrt(1-rho) * data_symbol ...
                          + sqrt(rho) * repmat(tag_symbol,K,1,1) * exp(1i*theta);

                data = [ones(K,1,block), tx_symbol];

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

                % Message detection under the message-power scaling
                d_message = abs(R_12 - sqrt(1-rho) * Constellation_sum);
                [~,index_message] = min(d_message,[],1);

                data_estimated_sum = Constellation_sum(index_message).';
                error_message_temp(times) = nnz(round(data_estimated_sum - data_sent,12));

                % Group tag detection from the residual signal
                residual_signal = R_12 - sqrt(1-rho) * data_estimated_sum;

                tag_stat = real(residual_signal * exp(-1i*theta));
                tag_hat = tag_stat > 0;

                tag_bit_vec = reshape(tag_bit,1,message_length);
                error_tag_temp(times) = nnz(tag_hat ~= tag_bit_vec);

            end

            message_fail = error_message_temp > 0;
            authentication_fail = error_tag_temp > eta;

            AFP_rho{index_K}(index_rho,index_RF) = ...
                mean(message_fail | authentication_fail);

            case_count = case_count + 1;
            progress_percent = case_count / total_cases * 100;

            fprintf('Progress: %.2f%%, K = %d, M = %d, rho = %.4f, AFP = %.3e\n', ...
                    progress_percent, K, RF_chain, rho, ...
                    AFP_rho{index_K}(index_rho,index_RF));

        end

        [AFP{index_K}(index_RF),rho_index] = min(AFP_rho{index_K}(:,index_RF));
        rho_opt{index_K}(index_RF) = rho_0(rho_index);

    end
end

% Plot AFP versus M with optimized rho
figure;
set(gcf,'unit','normalized','position',[0.1,0.1,0.2*1.2*2,0.3*1.2*2]);

panel_label = {'(a)','(b)','(c)','(d)'};

for index_K = 1:length(K_0)

    subplot(2,2,index_K);

    x_data = RF_chain_cell{index_K};
    y_data = AFP{index_K};

    semilogy(x_data,y_data,'-or','LineWidth',1,'MarkerSize',5);
    set(gca,'YScale','log');

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

% Plot AFP versus rho for each K and M
figure;
set(gcf,'unit','normalized','position',[0.1,0.1,0.2*1.2*2,0.3*1.2*2]);

for index_K = 1:length(K_0)

    subplot(2,2,index_K);

    rho_0 = rho_cell{index_K};
    M_list = RF_chain_cell{index_K};

    hold on; grid on; box on;

    for index_RF = 1:length(M_list)
        semilogy(rho_0, AFP_rho{index_K}(:,index_RF), ...
                 'LineWidth',1);
    end

    set(gca,'YScale','log');
    set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);

    xlim([min(rho_0), max(rho_0)]);

    y_data = AFP_rho{index_K};
    y_min = min(y_data(y_data > 0));
    if isempty(y_min)
        y_min = 1e-5;
    end
    ylim([y_min, 1]);

    xlabel('$\rho$','FontSize',14,'FontName','Times New Roman','Interpreter','latex');
    ylabel('AFP','FontSize',14,'FontName','Times New Roman','Interpreter','tex');

    title(['$K=',num2str(K_0(index_K)),'$'], ...
          'FontSize',14,'FontName','Times New Roman','Interpreter','latex');

    legend("M = " + string(M_list), ...
           'FontSize',9, ...
           'FontName','Times New Roman', ...
           'Location','best');

    text(0.5, -0.25, panel_label{index_K}, ...
         'Units','normalized', ...
         'HorizontalAlignment','center', ...
         'FontSize',14, ...
         'FontName','Times New Roman');

end

fprintf('Simulation time = %.6f s\n', toc);