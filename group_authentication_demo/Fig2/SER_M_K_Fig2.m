clc; clear all; tic;

% System parameters
theta = 45/180*pi;                         % Rotation angle
K_0 = [2 3 4 5];                           % Number of users
RF_chain_0 = 2*2.^(2:1:7);                 % Number of antennas/RF chains


% Channel and noise parameters
frequency_central = 3.5;                   % Carrier frequency in GHz
noise_figure = 9;                          % Noise figure in dB
Band_width = 1*10^6;                      % Bandwidth in Hz
Noise_thermal = 10^(-174/10)/1000;         % Thermal noise PSD in W/Hz
noise_sigma_square = Noise_thermal * Band_width * 10^(noise_figure/10);

% Distance and frame parameters
distance_range = [400 600];               % User distance range in meters
time_slots = 5;                            % One pilot slot plus four data slots
message_length = 160;                      % Number of message symbols
block = message_length/(time_slots-1);     % Number of blocks
run_times=1e5;

% SER storage
error_message_rate = zeros(length(K_0), length(RF_chain_0), run_times);
error_tag_rate = zeros(length(K_0), length(RF_chain_0), run_times);
SER_message = zeros(length(K_0), length(RF_chain_0));
SER_tag = zeros(length(K_0), length(RF_chain_0));

total_cases = length(K_0) * length(RF_chain_0);
case_count = 0;

for index_K = 1:length(K_0)

    K = K_0(index_K);

    % Generate constellation table for the current number of users
    Constellation_group_selection;

    % Sum constellation points for non-rotated and rotated cases
    Constellation_sum = Constellation_table(:,K+1);
    Constellation_table_rotated = Constellation_table * exp(-1i*theta);
    Constellation_sum_rotated = Constellation_table_rotated(:,K+1);

    % Uniformly distribute users within the distance range
    distance_link = linspace(distance_range(1), distance_range(2), K);

    % Path loss calculation
    P_L = 35.3*log10(distance_link) + 22.4 + 21.3*log10(frequency_central);
    P_L = 10.^(P_L/10);
    D = diag(P_L);   % The path loss P_L and the sqrt(D) pre-compensation in the encoding process cancel each other in the received signal.
    power_transmit = 10^(30/10)/1000 / mean(P_L); % actual P_T=30dBm
    for index_RF = 1:length(RF_chain_0)

        RF_chain = RF_chain_0(index_RF);

        % Temporary SER arrays for parfor slicing
        error_message_rate_temp = zeros(1,run_times);
        error_tag_rate_temp = zeros(1,run_times);

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
            error_message = round(reshape(data_estimated,size(data_sent)) - data_sent,12);
            error_message_rate_temp(times) = sum(abs(sign(error_message)),2) / message_length;

            % Tag detection by comparing two minimum distances
            tag_bit_vec = reshape(tag_bit,1,message_length);
            tag_bit_hat_vec = d_min_1 < d_min_2;

            % Tag symbol error rate
            error_tag = double(tag_bit_hat_vec) - double(tag_bit_vec);
            error_tag_rate_temp(times) = sum(abs(error_tag),2) / message_length;

        end

        % Average SER over Monte Carlo trials
        SER_message(index_K,index_RF) = mean(error_message_rate_temp);
        SER_tag(index_K,index_RF) = mean(error_tag_rate_temp);

        case_count = case_count + 1;
        progress_percent = case_count / total_cases * 100;

        fprintf('Progress: %.2f%% completed, K = %d, M = %d\n', ...
                progress_percent, K, RF_chain);

    end
end

fprintf('Plotting\n');

% Remove zero values from the log-scale plot
RF_chain_plot = RF_chain_0;
SER_message(SER_message == 0) = NaN;
SER_tag(SER_tag == 0) = NaN;

% Plot SER curves
figure;
set(gcf,'unit','normalized','position',[0.1,0.1,0.3*1.2,0.3*1.2]);
set(gca,'FontName','Times New Roman','FontSize',10,'LineWidth',1);
set(gca,'XScale','log','YScale','log');
grid on; hold on; box on;

loglog(RF_chain_plot,SER_message(1,:),'-or','LineWidth',1);
loglog(RF_chain_plot,SER_message(2,:),'-m^','LineWidth',1);
loglog(RF_chain_plot,SER_message(3,:),'-sb','LineWidth',1);
loglog(RF_chain_plot,SER_message(4,:),'-k*','LineWidth',1);

loglog(RF_chain_plot,SER_tag(1,:),'--or','LineWidth',1);
loglog(RF_chain_plot,SER_tag(2,:),'--m^','LineWidth',1);
loglog(RF_chain_plot,SER_tag(3,:),'--sb','LineWidth',1);
loglog(RF_chain_plot,SER_tag(4,:),'--k*','LineWidth',1);

ylim([1e-4,1]);

xlabel('$M$','FontSize',14,'FontName','Times New Roman','Interpreter','latex');
ylabel('SER','FontSize',14,'FontName','Times New Roman','Interpreter','tex');

% Add and manually shift the legend
lgd = legend('$K=2$,~Mes','$K=3$,~Mes','$K=4$,~Mes','$K=5$,~Mes', ...
       '$K=2$,~Tag','$K=3$,~Tag','$K=4$,~Tag','$K=5$,~Tag', ...
       'FontSize',14,'FontName','Times New Roman', ...
       'Location','southeast','Interpreter','latex');

set(lgd,'Box','on');
pos = get(lgd,'Position');
pos(1) = pos(1) - 0.17;
set(lgd,'Position',pos);

fprintf('Simulation time = %.6f s\n', toc);
%exportgraphics(gcf,'SER_M_K_Fig2.pdf','Resolution',600);