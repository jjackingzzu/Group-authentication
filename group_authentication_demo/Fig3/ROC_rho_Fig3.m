clc; clear all; tic;

% System parameters
theta_0 = (1:2:7)/180*pi;                         % Rotation angle

K_0 = [2 3 4 5];                           % Number of users
RF_chain_0 = [32 256 256 256];                 % Number of antennas/RF chains


% Channel and noise parameters
frequency_central = 3.5;                   % Carrier frequency in GHz
noise_figure = 9;                          % Noise figure in dB
Band_width = 10^6;                      % Bandwidth in Hz
Noise_thermal = 10^(-174/10)/1000;         % Thermal noise PSD in W/Hz
noise_sigma_square = Noise_thermal * Band_width * 10^(noise_figure/10);

% Distance and frame parameters
distance_range = [400 600];               % User distance range in meters
time_slots = 5;                            % One pilot slot plus four data slots
message_length = 160;                      % Number of message symbols
block = message_length/(time_slots-1);     % Number of blocks
run_times=1e5;

% error storage
error_message = zeros(length(K_0), length(theta_0), run_times);
error_tag = zeros(length(K_0), length(theta_0), run_times);

total_cases = length(K_0) * length(theta_0);
case_count = 0;

for index_K = 1:length(K_0)

    K = K_0(index_K);
    RF_chain=RF_chain_0(index_K);

    % Generate constellation table for the current number of users
    Constellation_group_selection;

    % Sum constellation points for non-rotated and rotated cases
    Constellation_sum = Constellation_table(:,K+1);
    

    % Uniformly distribute users within the distance range
    distance_link = linspace(distance_range(1), distance_range(2), K);

    % Path loss calculation
    P_L = 35.3*log10(distance_link) + 22.4 + 21.3*log10(frequency_central);
    P_L = 10.^(P_L/10);
    D = diag(P_L);% The path loss P_L and the sqrt(D) pre-compensation in the encoding process cancel each other in the received signal.
    power_transmit = 10^(30/10)/1000 / mean(P_L); % actual P_T=30dBm
    for index_theta = 1:length(theta_0)

        theta = theta_0(index_theta);

        Constellation_table_rotated = Constellation_table * exp(-1i*theta);
        Constellation_sum_rotated = Constellation_table_rotated(:,K+1);

        % Temporary SER arrays for parfor slicing
        error_message_temp = zeros(1,run_times);
        error_tag_temp = zeros(1,run_times);

        parfor times = 1:run_times

            % Generate random user data bits
            data_bit = round(rand(K,time_slots-1,block));

            % Map bits to non-rotated constellation symbols
            data_1 = data_bit .* Constellation_table(end,1:K).' ...
                   + (1-data_bit) .* Constellation_table(1,1:K).';

            % Map bits to rotated constellation symbols
            data_2 = data_bit .* Constellation_table_rotated(end,1:K).' ...
                   + (1-data_bit) .* Constellation_table_rotated(1,1:K).';

            % Generate one tag bit for each message symbol
            tag_bit = round(rand(1,time_slots-1,block));

            % Select rotated or non-rotated symbols according to tag bits
            data_rotated = data_1 .* tag_bit + data_2 .* (1-tag_bit);

            % Superimposed transmitted message symbols
            data_sent = reshape(sum(data_rotated,1),1,message_length);

            % Add the pilot slot before data slots
            data = [ones(K,1,block), data_rotated];

            % Generate Rayleigh fading channel
            H = randn(RF_chain,K,block)/sqrt(2) ...
              + 1i*randn(RF_chain,K,block)/sqrt(2);

            % Generate received signal with AWGN
            Y = sqrt(power_transmit) * pagemtimes(H,data) ...
              + (randn(RF_chain,time_slots,block) ...
              + 1i*randn(RF_chain,time_slots,block))/sqrt(2) ...
              * sqrt(noise_sigma_square);

            % Normalize the received signal by transmit power
            Y_tilde = Y / sqrt(power_transmit);

            % Pair the pilot slot with each data slot
            Y_tilde1 = zeros(RF_chain,2*(time_slots-1),block);
            Y_tilde1(:,1:2:2*time_slots-3,:) = repmat(Y_tilde(:,1,:),1,time_slots-1,1);
            Y_tilde1(:,2:2:2*time_slots-2,:) = Y_tilde(:,2:time_slots,:);
            Y_tilde1 = reshape(Y_tilde1,RF_chain,2,message_length);

            % Compute only the R(1,2) correlation term
            R_12 = reshape(sum(conj(Y_tilde1(:,1,:)) .* Y_tilde1(:,2,:),1) ...
                   / RF_chain, 1, message_length);

            % Nearest-neighbor detection under non-rotated hypothesis
            d_1 = abs(R_12 - Constellation_sum);
            [d_min_1,index_1] = min(d_1,[],1);
            data_estimated_1 = Constellation_sum(index_1).';

            % Nearest-neighbor detection under rotated hypothesis
            d_2 = abs(R_12 - Constellation_sum_rotated);
            [d_min_2,index_2] = min(d_2,[],1);
            data_estimated_2 = Constellation_sum_rotated(index_2).';

            % Select the corresponding estimate using the known tag bits
            data_estimated = reshape(data_estimated_1,1,time_slots-1,block) .* tag_bit ...
                           + reshape(data_estimated_2,1,time_slots-1,block) .* (1-tag_bit);

            % Message symbol error rate
            error_message_temp(times) = nnz(round(reshape(data_estimated,size(data_sent)) - data_sent,12));

            % Tag detection by comparing two minimum distances
            tag_bit_vec = reshape(tag_bit,1,message_length);
            tag_bit_hat_vec = d_min_1 < d_min_2;

            % Tag symbol error rate
            error_tag_temp(times) = nnz(tag_bit_hat_vec ~= tag_bit_vec);

        end

        % error over Monte Carlo trials
        error_message(index_K,index_theta,:) = error_message_temp;
        error_tag(index_K,index_theta,:) = error_tag_temp;

        case_count = case_count + 1;
        progress_percent = case_count / total_cases * 100;

        fprintf('Progress: %.2f%% completed, K = %d, theta = %d\n', ...
                progress_percent, K, theta/pi*180);

    end
end

% ROC plotting based on the already computed error_tag array

threshold_0 = 0:message_length;

P_d = zeros(length(K_0), length(theta_0), length(threshold_0));
P_f = zeros(1, length(threshold_0));

% False authentication probability for random tag guessing
% attack_error ~ Binomial(message_length, 0.5)
P_f_density = zeros(1, length(threshold_0));

for i = 1:length(threshold_0)

    eta_threshold = threshold_0(i);

    % Legitimate user authentication success probability:
    % success if tag error bits <= threshold
    P_d(:,:,i) = mean(error_tag <= eta_threshold, 3);

    % Attacker authentication success probability:
    % random guessing has exactly i-1 errors with probability C(N,i-1)/2^N
    P_f_density(i) = nchoosek(message_length, i-1);
    P_f(i) = sum(P_f_density) / 2^message_length;

end

figure;
set(gcf,'unit','normalized','position',[0.1,0.1,0.2*1.2*2,0.3*1.2*2]);

panel_label = {'(a)','(b)','(c)','(d)'};

for index_K = 1:length(K_0)

    subplot(2,2,index_K);

    plot(P_f, squeeze(P_d(index_K,1,:)),'-sb','LineWidth',1, ...
         'MarkerIndices',65:5:90); hold on;
    plot(P_f, squeeze(P_d(index_K,2,:)),'-^','LineWidth',1, ...
         'MarkerIndices',65:5:90);
    plot(P_f, squeeze(P_d(index_K,3,:)),'-or','LineWidth',1, ...
         'MarkerIndices',65:5:90);
    plot(P_f, squeeze(P_d(index_K,4,:)),'-k*','LineWidth',1, ...
         'MarkerIndices',65:5:90);

    grid on; box on;
    set(gca,'FontName','Times New Roman','FontSize',14,'LineWidth',0.5);

    legend('$\theta=1^{\circ}$','$\theta=3^{\circ}$', ...
           '$\theta=5^{\circ}$','$\theta=7^{\circ}$', ...
           'FontSize',14,'FontName','Times New Roman', ...
           'Location','southeast','Interpreter','latex');

    xlabel('$P_{\rm{f}}$','FontSize',14, ...
           'FontName','Times New Roman','Interpreter','latex');

    ylabel('$P_{\rm{d}}$','FontSize',14, ...
           'FontName','Times New Roman','Interpreter','latex');


    text(0.5, -0.25, panel_label{index_K}, ...
         'Units','normalized', ...
         'HorizontalAlignment','center', ...
         'FontSize',14, ...
         'FontName','Times New Roman');

end
%exportgraphics(gcf,'ROC_theta_Fig3.pdf','Resolution',600);