
clc; clear; tic;

%% ============================================================
%  CHT-PLA SNR Verification
%
%  Numerical SNR:
%       E{|s|^2}/E{|z-s|^2}
%
%  Theoretical SNR:
%       E{|s|^2}/Var(z)
%
%  Detection statistic:
%       z = [R]_{1,2}
%       R = Y^H Y / M
%% ============================================================

%% System parameters

K_set = 2:5;
M_set = 2*2.^(1:7);

      % W
frequency_central = 3.5;               % GHz
noise_figure = 9;                      % dB
bandwidth = 1e6;                      % Hz

message_length = 160;
time_slots = 5;
block = message_length/(time_slots-1);

run_times = 1e4;

%% Noise power

Noise_thermal = 10^((-174-30)/10);     % W/Hz

noise_sigma_square = ...
    Noise_thermal * bandwidth * ...
    10^(noise_figure/10);

%% Pre-allocation

SNR_theory = zeros(length(K_set),length(M_set));
SNR_numerical = zeros(length(K_set),length(M_set));
SER_Eu_d = zeros(length(K_set),length(M_set));

total_cases = length(K_set)*length(M_set);
case_count = 0;

%% Main loop

for index_K = 1:length(K_set)

    K = K_set(index_K);

    Constellation_group_selection;

    Constellation_sum = Constellation_table(:,K+1);

    %% Distance setting

    distance_link = linspace(400,600,K);

    P_L = 35.3*log10(distance_link) ...
        + 22.4 ...
        + 21.3*log10(frequency_central);

    P_L = 10.^(P_L/10);
    P_t=30; % dBm actual P_T
    power_transmit = 10^(P_t/10)/1000 / mean(P_L); 
    D = diag(1./P_L);

    % Large-scale fading is normalized and therefore
    % does not affect the resulting detection SNR.

    for index_M = 1:length(M_set)

        M = M_set(index_M);

        %% =====================================================
        % Theoretical SNR
        %% =====================================================

        signal_power_theory = ...
            mean(abs(Constellation_sum).^2);

        Var_theory = ...
            1/M * ...
            ( ...
            K * sum(abs(Constellation_table(:,1:K)).^2,2) ...
            + (noise_sigma_square/power_transmit)^2 ...
            + K*noise_sigma_square/power_transmit ...
            + noise_sigma_square/power_transmit ...
            .*sum(abs(Constellation_table(:,1:K)).^2,2) ...
            );

        SNR_theory(index_K,index_M) = ...
            10*log10( ...
            signal_power_theory/mean(Var_theory) ...
            );

        %% =====================================================
        % Numerical SNR
        %% =====================================================

        SNR_temp = zeros(1,run_times);
        SER_temp = zeros(1,run_times);

        parfor trial = 1:run_times

            %% Generate message

            data_bit = round(rand(K,time_slots-1,block));

            data_symbol = ...
                data_bit ...
                .* Constellation_table(end,1:K).' ...
                + (1-data_bit) ...
                .* Constellation_table(1,1:K).';

            data_sent = ...
                reshape(sum(data_symbol,1),1,message_length);

            pilot_symbol = ones(K,1,block);

            tx_data = [pilot_symbol,data_symbol];

            %% Channel

            H = ...
                randn(M,K,block)/sqrt(2) ...
                + 1i*randn(M,K,block)/sqrt(2);

            %% Received signal

            noise = ...
                (randn(M,time_slots,block) ...
                + 1i*randn(M,time_slots,block))/sqrt(2) ...
                * sqrt(noise_sigma_square);

            Y = ...
                sqrt(power_transmit) ...
                * pagemtimes(H,tx_data) ...
                + noise;

            %% Normalization

            Y_tilde = Y/sqrt(power_transmit);

            %% Construct matrix

            Y_pair = zeros(M,2*(time_slots-1),block);

            Y_pair(:,1:2:end,:) = ...
                repmat(Y_tilde(:,1,:),1,time_slots-1,1);

            Y_pair(:,2:2:end,:) = ...
                Y_tilde(:,2:time_slots,:);

            Y_pair = ...
                reshape(Y_pair,M,2,message_length);

            %% Detection statistic

            R12 = reshape( ...
                sum( ...
                conj(Y_pair(:,1,:)) ...
                .* Y_pair(:,2,:), ...
                1)/M,...
                1,...
                message_length);

            %% ML detection

            distance_metric = ...
                abs( ...
                repmat(R12,2^K,1) ...
                - repmat(Constellation_sum,1,message_length));

            [~,index_min] = min(distance_metric,[],1);

            data_estimated = ...
                Constellation_sum(index_min).';

            %% SNR

            SNR_temp(trial) = ...
                mean(abs(data_sent).^2) ...
                / mean(abs(R12-data_sent).^2);

            %% SER

            error_symbol = ...
                data_estimated - data_sent;

            SER_temp(trial) = ...
                sum(abs(sign(error_symbol))) ...
                / message_length;

        end

        SNR_numerical(index_K,index_M) = ...
            10*log10(mean(SNR_temp));

        SER_Eu_d(index_K,index_M) = ...
            mean(SER_temp);

        case_count = case_count + 1;

        fprintf( ...
            'Progress %.2f%% | K=%d | M=%d\n',...
            100*case_count/total_cases,...
            K,...
            M);

    end
end

%% Plot SER

figure;

semilogy(M_set,SER_Eu_d.','LineWidth',1.2);

grid on;
box on;

xlabel('$M$','Interpreter','latex');
ylabel('SER','Interpreter','latex');

legend('$K=2$','$K=3$','$K=4$','$K=5$', ...
       'Interpreter','latex');

fprintf('Simulation time = %.2f s\n',toc);

filename = 'SNR.xlsx';

% first page：Theoretical SNR
writematrix(SNR_theory,...
            filename,...
            'Sheet','Theoretical SNR');

% second：Numerical SNR
writematrix(SNR_numerical,...
            filename,...
            'Sheet','Numerical SNR');