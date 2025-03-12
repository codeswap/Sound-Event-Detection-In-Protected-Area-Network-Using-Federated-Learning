% solve the error of the mismatching of phi and data matrix [DONE]
% make static sink node [DONE]
% plotting [there are some more plotting remaining in the code] [DONE] 
% algorithm for the sound file compression [DONE]
% read a file and convert it to a vector [DONE]
% choosing a block size [DONE]
% changing compression percentages [DONE]
% initializing compressed matrices [DONE]
% actual compression [DONE]
% plotting audio signals [DONE]
% expanded view of audio signals [DONE]
% spectrogram of audio signals [DONE]

%Parameters

% Parameters
xm = 100;
ym = 100;
n = 100;
p = 0.05;
Eo = 0.050;
ETX = 50 * 0.000000001;
ERX = 50 * 0.000000001;
Efs = 10 * 0.000000000001;
Emp = 0.0013 * 0.000000000001;
EDA = 5 * 0.000000001;
rmax = 30;
do = sqrt(Efs / Emp);
m = 10;
a = 0.1;

% Read audio file and convert it to a vector
sig = audioread('music.mp3');
sig = sig(5500:100000);

% Set sample rate for audio
f = 44100; % Set your sample rate accordingly

% Calculate data size for each sensor node
data_size = floor(length(sig) / n);

% Compression function
compressionFactor = 4; % You can change this value as needed

windowSize = 8192;

switch compressionFactor
    case 2
        samplesToKeep = windowSize / 2;
    case 4
        samplesToKeep = windowSize / 4;
    case 8
        samplesToKeep = windowSize / 8;
    otherwise
        error('Invalid compression factor. Please choose 2, 4, or 8.');
end

% Initialize figure for plotting
figure(1);
hold off;

% Plotting the Region of Interest (ROI) with red color and static sink nodes at the vertices
vertices = [0 0; 0 ym; xm 0; xm ym];
plot(vertices(:, 1), vertices(:, 2), 'r-', 'LineWidth', 2);
hold on;

% Define the coordinates of the static sink nodes at the vertices
sink_coords = [0, 0; 0, ym; xm, 0; xm, ym];

% Plot static sink nodes at the vertices of the plot
for i = 1:size(sink_coords, 1)
    plot(sink_coords(i, 1), sink_coords(i, 2), 's', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
end

% Place sensor nodes randomly within the region of interest
XR = zeros(1, n);
YR = zeros(1, n);
for i = 1:n
    % Randomly place sensor nodes
    S(i).xd = rand(1, 1) * xm;
    XR(i) = S(i).xd;
    S(i).yd = rand(1, 1) * ym;
    YR(i) = S(i).yd;
    S(i).G = 0;
    S(i).type = 'N';

    % Assign audio data to each sensor node
    start_idx = (i-1) * data_size + 1;
    end_idx = i * data_size;
    S(i).data = sig(start_idx:end_idx);

    % Initialize neighbor data and compression data
    S(i).neigh_data = zeros(20 * data_size, 1);
    S(i).compress_data = zeros(data_size, 1);
    S(i).count = 1;

    % Energy initialization based on your conditions
    temp_rnd0 = rand;
    if (temp_rnd0 >= m * n + 1)
        S(i).E = Eo;
        S(i).ENERGY = 0;
        plot(S(i).xd, S(i).yd, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');
        hold on;
    else
        S(i).E = Eo * (1 + a);
        S(i).ENERGY = 1;
        plot(S(i).xd, S(i).yd, '+', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');
        hold on;
    end
end

% Place sink nodes at the vertices
for i = 1:size(sink_coords, 1)
    S(n+i).xd = sink_coords(i, 1);
    S(n+i).yd = sink_coords(i, 2);
end

countCHs = 0;
rcountCHs = 0;
cluster = 1;

countCHs;
rcountCHs = rcountCHs + countCHs;
flag_first_dead = 0;
first_dead_time = 0;

start_time = tic; % Start the timer

% Initialize arrays to store statistics
DEAD = zeros(1, rmax+1);
DEAD_N = zeros(1, rmax+1);
DEAD_A = zeros(1, rmax+1);
packets_to_CH = zeros(1, rmax+1);
packets_to_BS = zeros(1, rmax+1);
CLUSTERHS = zeros(1, rmax+1);
Energy_disp = zeros(1, rmax+1);
network_lifetime = zeros(1, rmax+1);
delay_with_compression = zeros(1, rmax+1);
delay_without_compression = zeros(1, rmax+1);
first_dead_round = -1;
last_dead_round = -1;

THROUGHPUT = zeros(1, rmax+1); % Initialize throughput array
RESIDUAL_ENERGY = zeros(1, rmax+1); % Initialize residual energy array

STATISTICS = struct();

% Initialize an array to store energy consumption of each node over rounds
energy_per_node = zeros(n, rmax+1);

for r = 0:rmax
    r;
    if (mod(r, round(1/p)) == 0)
        for i = 1:n
            S(i).G = 0;
            S(i).cl = 0;
        end
    end

    hold off;
    dead = 0;
    dead_a = 0;
    dead_n = 0;
    packets_to_BS(r+1) = 0;
    packets_to_CH(r+1) = 0;
    figure(1);

    for i = 1:n
        if (S(i).E <= 0)
            plot(S(i).xd, S(i).yd, '^', 'LineWidth', 1, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'MarkerSize', 8);
            dead = dead + 1;
            if (S(i).ENERGY == 1)
                dead_a = dead_a + 1;
            end
            if (S(i).ENERGY == 0)
                dead_n = dead_n + 1;
            end
            hold on;
        end

        if S(i).E > 0
            S(i).type = 'N';
            if (S(i).ENERGY == 0)
                plot(S(i).xd, S(i).yd, 'o', 'LineWidth', 1, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
            end
            if (S(i).ENERGY == 1)
                plot(S(i).xd, S(i).yd, '+', 'LineWidth', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
            end
            % Counting on end nodes and cluster heads
            text(S(i).xd, S(i).yd, num2str(S(i).count), 'Color', 'k', 'FontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            hold on;
        end

        % Store energy consumption for each node
        energy_per_node(i, r+1) = S(i).E;
    end

    STATISTICS(r+1).DEAD = dead;
    DEAD(r+1) = dead;
    DEAD_N(r+1) = dead_n;
    DEAD_A(r+1) = dead_a;

    if (dead == 1 && first_dead_round == -1)
        first_dead_round = r;
        first_dead_time = toc(start_time); % Record time when the first node dies
    end

    countCHs = 0;
    cluster = 1;
    for i = 1:n
        if (S(i).E > 0)
            temp_rand = rand;
            if ((S(i).G) <= 0)
                if (temp_rand <= (p / (1 - p * mod(r, round(1/p)))))
                    countCHs = countCHs + 1;
                    packets_to_BS(r+1) = packets_to_BS(r+1) + 1; % Increment packets sent to BS
                    S(i).type = 'C';
                    S(i).G = round(1/p) - 1;
                    C(cluster).xd = S(i).xd;
                    C(cluster).yd = S(i).yd;
                    plot(S(i).xd, S(i).yd, 'k*');

                    distance = sqrt((S(i).xd - S(n+1).xd)^2 + (S(i).yd - S(n+1).yd)^2);
                    C(cluster).distance = distance;
                    C(cluster).id = i;
                    X(cluster) = S(i).xd;
                    Y(cluster) = S(i).yd;
                    cluster = cluster + 1;

                    distance;
                    if (distance > do)
                        S(i).E = S(i).E - ((ETX + EDA) * (4000) + Emp * 4000 * (distance^4));
                    end
                    if (distance <= do)
                        S(i).E = S(i).E - ((ETX + EDA) * (4000) + Efs * 4000 * (distance^2));
                    end
                    Energy_disp(r+1) = S(i).E;
                end
            end
        end
    end

    STATISTICS(r+1).CLUSTERHEADS = cluster - 1;
    CLUSTERHS(r+1) = cluster - 1;

    % Plot lines from end nodes to cluster heads in yellow-green color
    for i = 1:n
        if (S(i).type == 'N' && S(i).E > 0)
            % Find the nearest cluster head
            min_distance = inf;
            nearest_CH = [];
            for j = 1:cluster-1
                distance_to_CH = sqrt((S(i).xd - C(j).xd)^2 + (S(i).yd - C(j).yd)^2);
                if distance_to_CH < min_distance
                    min_distance = distance_to_CH;
                    nearest_CH = j;
                end
            end
            % Plot dotted line from end node to the nearest cluster head
            line([S(i).xd, C(nearest_CH).xd], [S(i).yd, C(nearest_CH).yd], 'Color', 'g', 'LineStyle', '-.');
        end
    end

    % Plot lines from cluster heads to sink nodes in red color
    for j = 1:cluster-1
        % Find the nearest sink node
        min_distance_sink = inf;
        nearest_sink = [];
        for i = 1:size(sink_coords, 1)
            distance_to_sink = sqrt((C(j).xd - sink_coords(i, 1))^2 + (C(j).yd - sink_coords(i, 2))^2);
            if distance_to_sink < min_distance_sink
                min_distance_sink = distance_to_sink;
                nearest_sink = i;
            end
        end
        % Plot line from cluster head to the nearest sink node
        line([C(j).xd, sink_coords(nearest_sink, 1)], [C(j).yd, sink_coords(nearest_sink, 2)], 'Color', 'r');
    end

    axis([0 xm 0 ym]);

    % Calculate network lifetime
    network_lifetime(r+1) = toc(start_time);

    % Calculate throughput for the current round
    THROUGHPUT(r+1) = packets_to_BS(r+1) + packets_to_CH(r+1);

    % Calculate residual energy for the current round
    RESIDUAL_ENERGY(r+1) = sum([S(1:n).E]);

    % Check if all nodes are dead
    if dead == n && last_dead_round == -1
        last_dead_round = r; % Record round when the last node dies
        last_dead_time = toc(start_time);
    end

    % Calculate delay with and without compression
    if r > 0
        if first_dead_round ~= -1
            delay_with_compression(r+1) = network_lifetime(r+1) - network_lifetime(first_dead_round);
            delay_without_compression(r+1) = network_lifetime(r+1) - network_lifetime(first_dead_round);
        else
            % If first_dead_round is not set, delay_with_compression and delay_without_compression will remain zero
            delay_with_compression(r+1) = 0;
            delay_without_compression(r+1) = 0;
        end
    end

    % Print statistics for the current round
    disp(['Round ' num2str(r)]);
    disp(['  Alive Nodes: ' num2str(n - dead)]);
    disp(['  Dead Nodes: ' num2str(dead)]);
    disp(['  Packets to CH: ' num2str(packets_to_CH(r+1))]);
    disp(['  Packets to BS: ' num2str(packets_to_BS(r+1))]);
    disp(['  Total Energy Consumption: ' num2str(sum(energy_per_node(:, r+1)))]);
    disp('  Cluster Heads:');
    for k = 1:cluster-1
        disp(['    CH ' num2str(k) ': (' num2str(C(k).xd) ', ' num2str(C(k).yd) ')']);
    end
end

% Plotting graphs

% Plotting the number of dead nodes over rounds
figure;
plot(DEAD);
title('Number of Dead Nodes Over Rounds');
xlabel('Rounds');
ylabel('Number of Dead Nodes');

% Plotting the number of cluster heads over rounds
figure;
plot(CLUSTERHS);
title('Number of Cluster Heads Over Rounds');
xlabel('Rounds');
ylabel('Number of Cluster Heads');

% Plotting energy consumption
figure;
plot(Energy_disp);
title('Energy Consumption Over Rounds');
xlabel('Rounds');
ylabel('Energy');

% Plotting network lifetime
figure;
plot(network_lifetime);
title('Network Lifetime');
xlabel('Rounds');
ylabel('Time (s)');

% Plotting throughput over rounds
figure;
plot(THROUGHPUT);
title('Throughput Over Rounds');
xlabel('Rounds');
ylabel('Throughput');

% Plotting residual energy over rounds
figure;
plot(RESIDUAL_ENERGY);
title('Residual Energy Over Rounds');
xlabel('Rounds');
ylabel('Residual Energy (J)');

% Plotting the round at which the last node dies
if last_dead_round ~= -1
    figure;
    plot(last_dead_round, n, 'ro', 'MarkerSize', 10);
    title('Round of Last Node Death');
    xlabel('Round');
    ylabel('Number of Nodes');
else
    disp('No round found where all nodes are dead.');
end

% Plotting the round at which the first node dies
if first_dead_round ~= -1
    figure;
    plot(first_dead_round, 1, 'ro', 'MarkerSize', 10);
    title('Round of First Node Death');
    xlabel('Round');
    ylabel('Number of Nodes');
else
    disp('No round found where the first node dies.');
end

% Display time when the first node dies and when the last node dies
disp(['Round when first node dies: ', num2str(first_dead_round)]);
disp(['Round when last node dies: ', num2str(last_dead_round)]);
disp(['Time when first node dies: ', num2str(first_dead_time), ' seconds']);
if exist('last_dead_time', 'var')
    disp(['Time when last node dies: ', num2str(last_dead_time), ' seconds']);
else
    disp('No last_dead_time recorded.');
end

% Implement sound compression using DCT
for i = 1:n
    % Perform DCT compression for each sensor node
    windowDCT = dct(S(i).data);
    if samplesToKeep <= length(windowDCT)
        S(i).compress_data = idct(windowDCT(1:samplesToKeep), length(S(i).data));
    else
        % Handle the case where samplesToKeep exceeds the length of windowDCT
        S(i).compress_data = S(i).data;
    end
end

% Plot original and compressed audio signals
figure(2);
subplot(2,1,1);
plot(sig), title('Original Waveform');
subplot(2,1,2);
plot(S(1).compress_data), title(['Compressed Waveform - Compression Factor ', num2str(compressionFactor)]);

% Plot frequency spectrum of original waveform
figure(3);
plot_frequency_spectrum(sig, f, 'Original Waveform Frequency Spectrum');

% Plot frequency spectrum of compressed waveform
figure(4);
plot_frequency_spectrum(S(1).compress_data, f, ['Compressed Waveform Frequency Spectrum - Compression Factor ', num2str(compressionFactor)]);

% Plotting delay with and without compression
figure;
plot(delay_with_compression, 'b', 'LineWidth', 2);
hold on;
plot(delay_without_compression, 'r', 'LineWidth', 2);
title('Delay with versus without Compression');
xlabel('Rounds');
ylabel('Time (s)');
legend('With Compression', 'Without Compression');

% Plotting data size factor with 2, 4, and 8 with bar graph
data_size_factors = [2, 4, 8];
compressed_data_sizes = zeros(1, length(data_size_factors));
for i = 1:length(data_size_factors)
    compressed_data_sizes(i) = data_size / data_size_factors(i);
end
figure;
bar(data_size_factors, compressed_data_sizes);
title('Data Size Factor with Different Compression Factors');
xlabel('Compression Factor');
ylabel('Compressed Data Size');

% Plotting energy consumption per node over rounds
figure;
for i = 1:n
    if mod(i, 4) == 0 % Plot every 4th node to reduce clutter
        plot(0:rmax, energy_per_node(i, :), 'DisplayName', ['Node ' num2str(i)]);
        hold on;
    end
end
title('Energy Consumption per Node Over Rounds');
xlabel('Rounds');
ylabel('Energy');
legend('show', 'Location', 'northeastoutside');

% Calculate and plot the quality of the compressed signal using MSE and PSNR
compression_factors = [2, 4, 8];
mse_values = zeros(1, length(compression_factors));
psnr_values = zeros(1, length(compression_factors));

for i = 1:length(compression_factors)
    compressionFactor = compression_factors(i);
    samplesToKeep = windowSize / compressionFactor;
    compressed_signal = zeros(size(sig));
    for j = 1:n
        windowDCT = dct(S(j).data);
        if samplesToKeep <= length(windowDCT)
            compressed_signal((j-1)*data_size + 1:j*data_size) = idct(windowDCT(1:samplesToKeep), length(S(j).data));
        else
            compressed_signal((j-1)*data_size + 1:j*data_size) = S(j).data;
        end
    end
    mse_values(i) = immse(sig, compressed_signal);
    psnr_values(i) = psnr(sig, compressed_signal);
end

figure;
subplot(2, 1, 1);
bar(compression_factors, mse_values);
title('MSE for Different Compression Factors');
xlabel('Compression Factor');
ylabel('MSE');

subplot(2, 1, 2);
bar(compression_factors, psnr_values);
title('PSNR for Different Compression Factors');
xlabel('Compression Factor');
ylabel('PSNR (dB)');

% Define a function to plot frequency spectrum
function plot_frequency_spectrum(data, Fs, titleText)
    N = length(data);
    Y = fft(data);
    f = (0:N-1)*(Fs/N);
    power = abs(Y).^2/N;
    plot(f(1:floor(N/2)), power(1:floor(N/2)));
    xlabel('Frequency (Hz)');
    ylabel('Power');
    title(titleText);
end
