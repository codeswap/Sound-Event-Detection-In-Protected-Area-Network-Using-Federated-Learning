% Initialize common parameters
xm = 100;
ym = 100;
n = 150;
p = 0.05;
Eo = 0.050;
ETX = 50 * 0.000000001;
ERX = 50 * 0.000000001;
Efs = 10 * 0.000000000001;
Emp = 0.0013 * 0.000000000001;
EDA = 5 * 0.000000001;
rmax = 60;
do = sqrt(Efs / Emp);
m = 10;
a = 0.1;
data_size = 400;
sig = audioread('music.mp3');
sig = sig(5500:100000);
Q_STATISTICS = struct();
LEACH_STATISTICS = struct();
GAME_STATISTICS = struct();
PSO_STATISTICS = struct();
HEED_STATISTICS = struct();
ACO_STATISTICS = struct();

% Initialize nodes for all algorithms
S_leach = initialize_nodes(n, xm, ym, Eo, m, a, data_size, sig);
S_q = initialize_nodes(n, xm, ym, Eo, m, a, data_size, sig);
S_game = initialize_nodes(n, xm, ym, Eo, m, a, data_size, sig);
S_pso = initialize_nodes(n, xm, ym, Eo, m, a, data_size, sig);
S_heed = initialize_nodes(n, xm, ym, Eo, m, a, data_size, sig);
S_aco = initialize_nodes(n, xm, ym, Eo, m, a, data_size, sig);

% Run LEACH clustering
[LEACH_STATISTICS, sink_leach] = run_leach_clustering(S_leach, n, rmax, p, xm, ym, ETX, ERX, Efs, Emp, EDA, do, data_size);

% Run Q-learning based clustering
[Q_STATISTICS, sink_q] = run_q_learning_clustering(S_q, n, rmax, p, xm, ym, ETX, ERX, Efs, Emp, EDA, do, data_size);

% Run Game Theory based clustering
[GAME_STATISTICS, sink_game] = run_game_theory_clustering(S_game, n, rmax, p, xm, ym, ETX, ERX, Efs, Emp, EDA, do, data_size, Eo);

% Run PSO-based clustering
[PSO_STATISTICS, sink_pso] = run_pso_clustering(S_pso, n, rmax, p, xm, ym, ETX, ERX, Efs, Emp, EDA, do, data_size, Eo);

% Run HEED clustering
[HEED_STATISTICS, sink_heed] = run_heed_clustering(S_heed, n, rmax, p, xm, ym, ETX, ERX, Efs, Emp, EDA, do, data_size);

% Run ACO-based clustering
[ACO_STATISTICS, sink_aco] = run_aco_clustering(S_aco, n, rmax, p, xm, ym, ETX, ERX, Efs, Emp, EDA, do, data_size);

% Plot comparisons
plot_comparisons(LEACH_STATISTICS, Q_STATISTICS, GAME_STATISTICS, PSO_STATISTICS, HEED_STATISTICS, ACO_STATISTICS, rmax);

function S = initialize_nodes(n, xm, ym, Eo, m, a, data_size, sig)
    for i = 1:n
        S(i).xd = rand(1, 1) * xm;
        S(i).yd = rand(1, 1) * ym;
        S(i).G = 0;
        S(i).type = 'N';
        S(i).data = sig(((i - 1) * data_size) + 1:(i * data_size));
        S(i).neigh_data = zeros(20 * data_size, 1);
        S(i).compress_data = zeros(data_size, 1);
        S(i).count = 1;
        S(i).min_dis_cluster = 0;
        S(i).min_dis = 0;
        temp_rnd0 = i;
        if (temp_rnd0 >= m * n + 1)
            S(i).E = Eo;
            S(i).ENERGY = 0;
        else
            S(i).E = Eo * (1 + a);
            S(i).ENERGY = 1;
        end
    end
    S(n + 1).xd = xm / 2;
    S(n + 1).yd = ym / 2;
end

function particles = initialize_particles(n, xm, ym, num_particles)
    particles(num_particles).position = [];
    particles(num_particles).velocity = [];
    particles(num_particles).pbest = [];
    particles(num_particles).pbest_fitness = [];
    particles(num_particles).gbest = [];
    particles(num_particles).gbest_fitness = [];

    for i = 1:num_particles
        particles(i).position = [rand(1) * xm, rand(1) * ym];
        particles(i).velocity = [rand(1) - 0.5, rand(1) - 0.5];
        particles(i).pbest = particles(i).position;
        particles(i).pbest_fitness = -inf;
        particles(i).gbest = particles(i).position;
        particles(i).gbest_fitness = -inf;
    end
end

function fitness = evaluate_particle(position, energy, Eo)
    fitness = energy / Eo; % Fitness based on remaining energy
end

function [STATISTICS, sink] = run_leach_clustering(S, n, rmax, p, xm, ym, ETX, ERX, Efs, Emp, EDA, do, data_size)
    countCHs = 0;
    rcountCHs = 0;
    cluster = 1;
    flag_first_dead = 0;
    STATISTICS = struct();
    sink.x = xm / 2;
    sink.y = ym / 2;
    sink.data = zeros(data_size * n, 3);

    for r = 0:rmax
        if (mod(r, round(1 / p)) == 0)
            for i = 1:n
                S(i).G = 0;
                S(i).cl = 0;
            end
        end
        dead = 0;
        dead_a = 0;
        dead_n = 0;
        packets_TO_BS = 0;
        packets_TO_CH = 0;

        for i = 1:n
            if (S(i).E <= 0)
                dead = dead + 1;
                if (S(i).ENERGY == 1)
                    dead_a = dead_a + 1;
                else
                    dead_n = dead_n + 1;
                end
            end
        end
        STATISTICS(r + 1).DEAD = dead;
        STATISTICS(r + 1).DEAD_N = dead_n;
        STATISTICS(r + 1).DEAD_A = dead_a;

        if (dead == 1 && flag_first_dead == 0)
            first_dead = r;
            flag_first_dead = 1;
        end

        countCHs = 0;
        cluster = 1;
        for i = 1:n
            if (S(i).E > 0)
                temp_rand = rand;
                if ((S(i).G) <= 0)
                    if (temp_rand <= (p / (1 - p * mod(r, round(1 / p)))))
                        countCHs = countCHs + 1;
                        packets_TO_BS = packets_TO_BS + 1;
                        STATISTICS(r + 1).PACKETS_TO_BS = packets_TO_BS;
                        S(i).type = 'C';
                        S(i).G = round(1 / p) - 1;
                        C(cluster).xd = S(i).xd;
                        C(cluster).yd = S(i).yd;
                        distance = sqrt((S(i).xd - sink.x)^2 + (S(i).yd - sink.y)^2);
                        C(cluster).distance = distance;
                        C(cluster).id = i;
                        cluster = cluster + 1;
                        if (distance > do)
                            S(i).E = S(i).E - ((ETX + EDA) * (4000) + Emp * 4000 * (distance^4));
                        else
                            S(i).E = S(i).E - ((ETX + EDA) * (4000) + Efs * 4000 * (distance^2));
                        end
                    end
                end
            end
        end
        STATISTICS(r + 1).CLUSTERHEADS = cluster - 1;

        for i = 1:n
            if (S(i).type == 'N' && S(i).E > 0)
                if (cluster - 1 >= 1)
                    min_dis = sqrt((S(i).xd - sink.x)^2 + (S(i).yd - sink.y)^2);
                    min_dis_cluster = 1;
                    for c = 1:cluster - 1
                        temp = min(min_dis, sqrt((S(i).xd - C(c).xd)^2 + (S(i).yd - C(c).yd)^2));
                        if (temp < min_dis)
                            min_dis = temp;
                            min_dis_cluster = c;
                        end
                    end
                    if (min_dis > do)
                        S(i).E = S(i).E - (ETX * (4000) + Emp * 4000 * (min_dis^4));
                    else
                        S(i).E = S(i).E - (ETX * (4000) + Efs * 4000 * (min_dis^2));
                    end
                    if (min_dis > 0)
                        distance = sqrt((S(C(min_dis_cluster).id).xd - sink.x)^2 + (S(C(min_dis_cluster).id).yd - sink.y)^2);
                        S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ERX + EDA) * 4000);
                        if (distance > do)
                            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ETX + EDA) * (4000) + Emp * 4000 * (distance^4));
                        else
                            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ETX + EDA) * (4000) + Efs * 4000 * (distance^2));
                        end
                        packets_TO_CH = n - dead - cluster + 1;
                    end
                    S(i).min_dis = min_dis;
                    S(i).min_dis_cluster = min_dis_cluster;
                end
            end
        end
        STATISTICS(r + 1).PACKETS_TO_CH = packets_TO_CH;

        su = 0;
        for i = 1:n
            if (S(i).E > 0)
                su = su + S(i).E;
            end
        end
        avg = su / n;
        STATISTICS(r + 1).AVG = avg;

        for i = 1:n
            if S(i).min_dis_cluster
                pt = C(S(i).min_dis_cluster).id;
            else
                pt = i;
            end
            if S(pt).count
                N = 400;
                M = 50;
                M = ceil(M);
                if N > M
                    rng(i);
                    Phi = randn(M, N);
                    S(i).compress_data = Phi * S(i).data';
                else
                    S(i).compress_data = S(i).data;
                end
                len = size(S(i).compress_data);
                S(pt).neigh_data(S(pt).count:S(pt).count + len - 1, 1) = S(i).compress_data;
                S(pt).neigh_data(S(pt).count + data_size - 1, 1) = i;
                S(pt).count = S(pt).count + data_size;
            end
        end

        count = 1;
        for i = 0:length(C) - 1
            clust_point = C(i + 1).id;
            data = S(clust_point).neigh_data;
            for i = 1:data_size:S(clust_point).count - 400
                node_id = data(399 + i, 1);
                sink.data(((node_id - 1) * data_size) + 1:(node_id * data_size), 2) = data(i:i + 399, 1);
            end
        end

        stoping_point = find(sink.data(:, 1) == -1);
        start = 1;
        for j = 1:n
            data = sink.data(start:start + M - 1, 2);
            node_id = sink.data(start + data_size - 1, 2);
            re_x = reconstruct(data, N, node_id);
            sink.data(start:start + data_size - 1, 3) = re_x;
            start = start + data_size;
        end
    end
end

function [STATISTICS, sink] = run_q_learning_clustering(S, n, rmax, p, xm, ym, ETX, ERX, Efs, Emp, EDA, do, data_size)
    % Q-learning parameters
    alpha = 0.5;
    beta = 0.5;
    STATISTICS = struct();
    sink.x = xm / 2;
    sink.y = ym / 2;
    sink.data = zeros(data_size * n, 3);

    for i = 1:n
        Q(i).qi = rand();
        Q(i).last_transmission_round = 0;
    end

    for r = 0:rmax
        if (mod(r, round(1 / p)) == 0)
            for i = 1:n
                S(i).G = 0;
                S(i).cl = 0;
            end
        end
        dead = 0;
        dead_a = 0;
        dead_n = 0;
        packets_TO_BS = 0;
        packets_TO_CH = 0;

        for i = 1:n
            if (S(i).E <= 0)
                dead = dead + 1;
                if (S(i).ENERGY == 1)
                    dead_a = dead_a + 1;
                else
                    dead_n = dead_n + 1;
                end
            end
        end
        STATISTICS(r + 1).DEAD = dead;
        STATISTICS(r + 1).DEAD_N = dead_n;
        STATISTICS(r + 1).DEAD_A = dead_a;

        if (dead == 1 && flag_first_dead == 0)
            first_dead = r;
            flag_first_dead = 1;
        end

        countCHs = 0;
        cluster = 1;
        for i = 1:n
            if (S(i).E > 0)
                temp_rand = rand;
                if ((S(i).G) <= 0)
                    if (temp_rand <= Q(i).qi)
                        countCHs = countCHs + 1;
                        packets_TO_BS = packets_TO_BS + 1;
                        STATISTICS(r + 1).PACKETS_TO_BS = packets_TO_BS;
                        S(i).type = 'C';
                        S(i).G = round(1 / p) - 1;
                        C(cluster).xd = S(i).xd;
                        C(cluster).yd = S(i).yd;
                        distance = sqrt((S(i).xd - sink.x)^2 + (S(i).yd - sink.y)^2);
                        C(cluster).distance = distance;
                        C(cluster).id = i;
                        cluster = cluster + 1;
                        if (distance > do)
                            S(i).E = S(i).E - ((ETX + EDA) * (4000) + Emp * 4000 * (distance^4));
                        else
                            S(i).E = S(i).E - ((ETX + EDA) * (4000) + Efs * 4000 * (distance^2));
                        end
                        Q(i).qi = alpha * Q(i).qi + beta * (1 - Q(i).qi);
                    end
                end
            end
        end
        STATISTICS(r + 1).CLUSTERHEADS = cluster - 1;

        for i = 1:n
            if (S(i).type == 'N' && S(i).E > 0)
                if (cluster - 1 >= 1)
                    min_dis = sqrt((S(i).xd - sink.x)^2 + (S(i).yd - sink.y)^2);
                    min_dis_cluster = 1;
                    for c = 1:cluster - 1
                        temp = min(min_dis, sqrt((S(i).xd - C(c).xd)^2 + (S(i).yd - C(c).yd)^2));
                        if (temp < min_dis)
                            min_dis = temp;
                            min_dis_cluster = c;
                        end
                    end
                    if (min_dis > do)
                        S(i).E = S(i).E - (ETX * (4000) + Emp * 4000 * (min_dis^4));
                    else
                        S(i).E = S(i).E - (ETX * (4000) + Efs * 4000 * (min_dis^2));
                    end
                    if (min_dis > 0)
                        distance = sqrt((S(C(min_dis_cluster).id).xd - sink.x)^2 + (S(C(min_dis_cluster).id).yd - sink.y)^2);
                        S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ERX + EDA) * 4000);
                        if (distance > do)
                            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ETX + EDA) * (4000) + Emp * 4000 * (distance^4));
                        else
                            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ETX + EDA) * (4000) + Efs * 4000 * (distance^2));
                        end
                        packets_TO_CH = n - dead - cluster + 1;
                    end
                    S(i).min_dis = min_dis;
                    S(i).min_dis_cluster = min_dis_cluster;
                end
            end
        end
        STATISTICS(r + 1).PACKETS_TO_CH = packets_TO_CH;

        su = 0;
        for i = 1:n
            if (S(i).E > 0)
                su = su + S(i).E;
            end
        end
        avg = su / n;
        STATISTICS(r + 1).AVG = avg;

        for i = 1:n
            if S(i).min_dis_cluster
                pt = C(S(i).min_dis_cluster).id;
            else
                pt = i;
            end
            if S(pt).count
                N = 400;
                M = 50;
                M = ceil(M);
                if N > M
                    rng(i);
                    Phi = randn(M, N);
                    S(i).compress_data = Phi * S(i).data';
                else
                    S(i).compress_data = S(i).data;
                end
                len = size(S(i).compress_data);
                S(pt).neigh_data(S(pt).count:S(pt).count + len - 1, 1) = S(i).compress_data;
                S(pt).neigh_data(S(pt).count + data_size - 1, 1) = i;
                S(pt).count = S(pt).count + data_size;
            end
        end

        count = 1;
        for i = 0:length(C) - 1
            clust_point = C(i + 1).id;
            data = S(clust_point).neigh_data;
            for i = 1:data_size:S(clust_point).count - 400
                node_id = data(399 + i, 1);
                sink.data(((node_id - 1) * data_size) + 1:(node_id * data_size), 2) = data(i:i + 399, 1);
            end
        end

        stoping_point = find(sink.data(:, 1) == -1);
        start = 1;
        for j = 1:n
            data = sink.data(start:start + M - 1, 2);
            node_id = sink.data(start + data_size - 1, 2);
            re_x = reconstruct(data, N, node_id);
            sink.data(start:start + data_size - 1, 3) = re_x;
            start = start + data_size;
        end
    end
end

function [STATISTICS, sink] = run_game_theory_clustering(S, n, rmax, p, xm, ym, ETX, ERX, Efs, Emp, EDA, do, data_size, Eo)
    STATISTICS = struct();
    sink.x = xm / 2;
    sink.y = ym / 2;
    sink.data = zeros(data_size * n, 3);

    for r = 0:rmax
        if (mod(r, round(1 / p)) == 0)
            for i = 1:n
                S(i).G = 0;
                S(i).cl = 0;
            end
        end
        dead = 0;
        dead_a = 0;
        dead_n = 0;
        packets_TO_BS = 0;
        packets_TO_CH = 0;

        for i = 1:n
            if (S(i).E <= 0)
                dead = dead + 1;
                if (S(i).ENERGY == 1)
                    dead_a = dead_a + 1;
                else
                    dead_n = dead_n + 1;
                end
            end
        end
        STATISTICS(r + 1).DEAD = dead;
        STATISTICS(r + 1).DEAD_N = dead_n;
        STATISTICS(r + 1).DEAD_A = dead_a;

        countCHs = 0;
        cluster = 1;
        for i = 1:n
            if (S(i).E > 0)
                % Game theory decision-making for cluster head selection
                temp_rand = rand;
                utility = S(i).E / Eo; % Simple utility function based on energy level
                if ((S(i).G) <= 0)
                    if (temp_rand <= (p / (1 - p * mod(r, round(1/p)))) * utility)
                        countCHs = countCHs + 1;
                        packets_TO_BS = packets_TO_BS + 1;
                        STATISTICS(r + 1).PACKETS_TO_BS = packets_TO_BS;
                        S(i).type = 'C';
                        S(i).G = round(1 / p) - 1;
                        C(cluster).xd = S(i).xd;
                        C(cluster).yd = S(i).yd;
                        distance = sqrt((S(i).xd - sink.x)^2 + (S(i).yd - sink.y)^2);
                        C(cluster).distance = distance;
                        C(cluster).id = i;
                        cluster = cluster + 1;
                        if (distance > do)
                            S(i).E = S(i).E - ((ETX + EDA) * (4000) + Emp * 4000 * (distance^4));
                        else
                            S(i).E = S(i).E - ((ETX + EDA) * (4000) + Efs * 4000 * (distance^2));
                        end
                    end
                end
            end
        end
        STATISTICS(r + 1).CLUSTERHEADS = cluster - 1;

        for i = 1:n
            if (S(i).type == 'N' && S(i).E > 0)
                if (cluster - 1 >= 1)
                    min_dis = sqrt((S(i).xd - sink.x)^2 + (S(i).yd - sink.y)^2);
                    min_dis_cluster = 1;
                    for c = 1:cluster - 1
                        temp = min(min_dis, sqrt((S(i).xd - C(c).xd)^2 + (S(i).yd - C(c).yd)^2));
                        if (temp < min_dis)
                            min_dis = temp;
                            min_dis_cluster = c;
                        end
                    end
                    if (min_dis > do)
                        S(i).E = S(i).E - (ETX * (4000) + Emp * 4000 * (min_dis^4));
                    else
                        S(i).E = S(i).E - (ETX * (4000) + Efs * 4000 * (min_dis^2));
                    end
                    if (min_dis > 0)
                        distance = sqrt((S(C(min_dis_cluster).id).xd - sink.x)^2 + (S(C(min_dis_cluster).id).yd - sink.y)^2);
                        S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ERX + EDA) * 4000);
                        if (distance > do)
                            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ETX + EDA) * (4000) + Emp * 4000 * (distance^4));
                        else
                            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ETX + EDA) * (4000) + Efs * 4000 * (distance^2));
                        end
                        packets_TO_CH = n - dead - cluster + 1;
                    end
                    S(i).min_dis = min_dis;
                    S(i).min_dis_cluster = min_dis_cluster;
                end
            end
        end
        STATISTICS(r + 1).PACKETS_TO_CH = packets_TO_CH;

        su = 0;
        for i = 1:n
            if (S(i).E > 0)
                su = su + S(i).E;
            end
        end
        avg = su / n;
        STATISTICS(r + 1).AVG = avg;

        for i = 1:n
            if S(i).min_dis_cluster
                pt = C(S(i).min_dis_cluster).id;
            else
                pt = i;
            end
            if S(pt).count
                N = 400;
                M = 50;
                M = ceil(M);
                if N > M
                    rng(i);
                    Phi = randn(M, N);
                    S(i).compress_data = Phi * S(i).data';
                else
                    S(i).compress_data = S(i).data;
                end
                len = size(S(i).compress_data);
                S(pt).neigh_data(S(pt).count:S(pt).count + len - 1, 1) = S(i).compress_data;
                S(pt).neigh_data(S(pt).count + data_size - 1, 1) = i;
                S(pt).count = S(pt).count + data_size;
            end
        end

        count = 1;
        for i = 0:length(C) - 1
            clust_point = C(i + 1).id;
            data = S(clust_point).neigh_data;
            for i = 1:data_size:S(clust_point).count - 400
                node_id = data(399 + i, 1);
                sink.data(((node_id - 1) * data_size) + 1:(node_id * data_size), 2) = data(i:i + 399, 1);
            end
        end

        stoping_point = find(sink.data(:, 1) == -1);
        start = 1;
        for j = 1:n
            data = sink.data(start:start + M - 1, 2);
            node_id = sink.data(start + data_size - 1, 2);
            re_x = reconstruct(data, N, node_id);
            sink.data(start:start + data_size - 1, 3) = re_x;
            start = start + data_size;
        end
    end
end

function [STATISTICS, sink] = run_pso_clustering(S, n, rmax, p, xm, ym, ETX, ERX, Efs, Emp, EDA, do, data_size, Eo)
    % PSO parameters
    num_particles = n;
    w = 0.5; % Inertia weight
    c1 = 1.5; % Cognitive (particle) weight
    c2 = 1.5; % Social (swarm) weight

    STATISTICS = struct();
    sink.x = xm / 2;
    sink.y = ym / 2;
    sink.data = zeros(data_size * n, 3);

    % Initialize particles
    particles = initialize_particles(n, xm, ym, num_particles);

    for r = 0:rmax
        if (mod(r, round(1 / p)) == 0)
            for i = 1:n
                S(i).G = 0;
                S(i).cl = 0;
            end
        end
        dead = 0;
        dead_a = 0;
        dead_n = 0;
        packets_TO_BS = 0;
        packets_TO_CH = 0;

        for i = 1:n
            if (S(i).E <= 0)
                dead = dead + 1;
                if (S(i).ENERGY == 1)
                    dead_a = dead_a + 1;
                else
                    dead_n = dead_n + 1;
                end
            end
        end
        STATISTICS(r + 1).DEAD = dead;
        STATISTICS(r + 1).DEAD_N = dead_n;
        STATISTICS(r + 1).DEAD_A = dead_a;

        countCHs = 0;
        cluster = 1;
        for i = 1:n
            if (S(i).E > 0)
                % PSO-based cluster head selection
                for p = 1:num_particles
                    % Update velocity
                    particles(p).velocity = w * particles(p).velocity + c1 * rand * (particles(p).pbest - particles(p).position) + c2 * rand * (particles(p).gbest - particles(p).position);
                    % Update position
                    particles(p).position = particles(p).position + particles(p).velocity;
                end

                % Evaluate particles
                for p = 1:num_particles
                    particles(p).fitness = evaluate_particle(particles(p).position, S(i).E, Eo);
                    if particles(p).fitness > particles(p).pbest_fitness
                        particles(p).pbest = particles(p).position;
                        particles(p).pbest_fitness = particles(p).fitness;
                    end
                    if particles(p).fitness > particles(p).gbest_fitness
                        particles(p).gbest = particles(p).position;
                        particles(p).gbest_fitness = particles(p).fitness;
                    end
                end

                % Select cluster head based on gbest
                if particles(1).gbest_fitness > rand
                    countCHs = countCHs + 1;
                    packets_TO_BS = packets_TO_BS + 1;
                    STATISTICS(r + 1).PACKETS_TO_BS = packets_TO_BS;
                    S(i).type = 'C';
                    S(i).G = round(1 / p) - 1;
                    C(cluster).xd = S(i).xd;
                    C(cluster).yd = S(i).yd;
                    distance = sqrt((S(i).xd - sink.x)^2 + (S(i).yd - sink.y)^2);
                    C(cluster).distance = distance;
                    C(cluster).id = i;
                    cluster = cluster + 1;
                    if (distance > do)
                        S(i).E = S(i).E - ((ETX + EDA) * (4000) + Emp * 4000 * (distance^4));
                    else
                        S(i).E = S(i).E - ((ETX + EDA) * (4000) + Efs * 4000 * (distance^2));
                    end
                end
            end
        end
        STATISTICS(r + 1).CLUSTERHEADS = cluster - 1;

        for i = 1:n
            if (S(i).type == 'N' && S(i).E > 0)
                if (cluster - 1 >= 1)
                    min_dis = sqrt((S(i).xd - sink.x)^2 + (S(i).yd - sink.y)^2);
                    min_dis_cluster = 1;
                    for c = 1:cluster - 1
                        temp = min(min_dis, sqrt((S(i).xd - C(c).xd)^2 + (S(i).yd - C(c).yd)^2));
                        if (temp < min_dis)
                            min_dis = temp;
                            min_dis_cluster = c;
                        end
                    end
                    if (min_dis > do)
                        S(i).E = S(i).E - (ETX * (4000) + Emp * 4000 * (min_dis^4));
                    else
                        S(i).E = S(i).E - (ETX * (4000) + Efs * 4000 * (min_dis^2));
                    end
                    if (min_dis > 0)
                        distance = sqrt((S(C(min_dis_cluster).id).xd - sink.x)^2 + (S(C(min_dis_cluster).id).yd - sink.y)^2);
                        S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ERX + EDA) * 4000);
                        if (distance > do)
                            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ETX + EDA) * (4000) + Emp * 4000 * (distance^4));
                        else
                            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ETX + EDA) * (4000) + Efs * 4000 * (distance^2));
                        end
                        packets_TO_CH = n - dead - cluster + 1;
                    end
                    S(i).min_dis = min_dis;
                    S(i).min_dis_cluster = min_dis_cluster;
                end
            end
        end
        STATISTICS(r + 1).PACKETS_TO_CH = packets_TO_CH;

        su = 0;
        for i = 1:n
            if (S(i).E > 0)
                su = su + S(i).E;
            end
        end
        avg = su / n;
        STATISTICS(r + 1).AVG = avg;

        for i = 1:n
            if S(i).min_dis_cluster
                pt = C(S(i).min_dis_cluster).id;
            else
                pt = i;
            end
            if S(pt).count
                N = 400;
                M = 50;
                M = ceil(M);
                if N > M
                    rng(i);
                    Phi = randn(M, N);
                    S(i).compress_data = Phi * S(i).data';
                else
                    S(i).compress_data = S(i).data;
                end
                len = size(S(i).compress_data);
                S(pt).neigh_data(S(pt).count:S(pt).count + len - 1, 1) = S(i).compress_data;
                S(pt).neigh_data(S(pt).count + data_size - 1, 1) = i;
                S(pt).count = S(pt).count + data_size;
            end
        end

        count = 1;
        for i = 0:length(C) - 1
            clust_point = C(i + 1).id;
            data = S(clust_point).neigh_data;
            for i = 1:data_size:S(clust_point).count - 400
                node_id = data(399 + i, 1);
                sink.data(((node_id - 1) * data_size) + 1:(node_id * data_size), 2) = data(i:i + 399, 1);
            end
        end

        stoping_point = find(sink.data(:, 1) == -1);
        start = 1;
        for j = 1:n
            data = sink.data(start:start + M - 1, 2);
            node_id = sink.data(start + data_size - 1, 2);
            re_x = reconstruct(data, N, node_id);
            sink.data(start:start + data_size - 1, 3) = re_x;
            start = start + data_size;
        end
    end
end

function [STATISTICS, sink] = run_heed_clustering(S, n, rmax, p, xm, ym, ETX, ERX, Efs, Emp, EDA, do, data_size)
    % HEED parameters
    cprob = 0.05; % Initial cluster head probability
    max_iter = 10; % Maximum number of iterations

    STATISTICS = struct();
    sink.x = xm / 2;
    sink.y = ym / 2;
    sink.data = zeros(data_size * n, 3);

    for r = 0:rmax
        if (mod(r, round(1 / p)) == 0)
            for i = 1:n
                S(i).G = 0;
                S(i).cl = 0;
            end
        end
        dead = 0;
        dead_a = 0;
        dead_n = 0;
        packets_TO_BS = 0;
        packets_TO_CH = 0;

        for i = 1:n
            if (S(i).E <= 0)
                dead = dead + 1;
                if (S(i).ENERGY == 1)
                    dead_a = dead_a + 1;
                else
                    dead_n = dead_n + 1;
                end
            end
        end
        STATISTICS(r + 1).DEAD = dead;
        STATISTICS(r + 1).DEAD_N = dead_n;
        STATISTICS(r + 1).DEAD_A = dead_a;

        countCHs = 0;
        cluster = 1;
        iter = 0;
        prob = cprob;

        while countCHs < 1 && iter < max_iter
            countCHs = 0;
            for i = 1:n
                if (S(i).E > 0)
                    temp_rand = rand;
                    if (S(i).G <= 0 && temp_rand <= prob)
                        countCHs = countCHs + 1;
                        packets_TO_BS = packets_TO_BS + 1;
                        STATISTICS(r + 1).PACKETS_TO_BS = packets_TO_BS;
                        S(i).type = 'C';
                        S(i).G = round(1 / p) - 1;
                        C(cluster).xd = S(i).xd;
                        C(cluster).yd = S(i).yd;
                        distance = sqrt((S(i).xd - sink.x)^2 + (S(i).yd - sink.y)^2);
                        C(cluster).distance = distance;
                        C(cluster).id = i;
                        cluster = cluster + 1;
                        if (distance > do)
                            S(i).E = S(i).E - ((ETX + EDA) * (4000) + Emp * 4000 * (distance^4));
                        else
                            S(i).E = S(i).E - ((ETX + EDA) * (4000) + Efs * 4000 * (distance^2));
                        end
                    end
                end
            end
            prob = prob * 2;
            iter = iter + 1;
        end
        STATISTICS(r + 1).CLUSTERHEADS = cluster - 1;

        for i = 1:n
            if (S(i).type == 'N' && S(i).E > 0)
                if (cluster - 1 >= 1)
                    min_dis = sqrt((S(i).xd - sink.x)^2 + (S(i).yd - sink.y)^2);
                    min_dis_cluster = 1;
                    for c = 1:cluster - 1
                        temp = min(min_dis, sqrt((S(i).xd - C(c).xd)^2 + (S(i).yd - C(c).yd)^2));
                        if (temp < min_dis)
                            min_dis = temp;
                            min_dis_cluster = c;
                        end
                    end
                    if (min_dis > do)
                        S(i).E = S(i).E - (ETX * (4000) + Emp * 4000 * (min_dis^4));
                    else
                        S(i).E = S(i).E - (ETX * (4000) + Efs * 4000 * (min_dis^2));
                    end
                    if (min_dis > 0)
                        distance = sqrt((S(C(min_dis_cluster).id).xd - sink.x)^2 + (S(C(min_dis_cluster).id).yd - sink.y)^2);
                        S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ERX + EDA) * 4000);
                        if (distance > do)
                            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ETX + EDA) * (4000) + Emp * 4000 * (distance^4));
                        else
                            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ETX + EDA) * (4000) + Efs * 4000 * (distance^2));
                        end
                        packets_TO_CH = n - dead - cluster + 1;
                    end
                    S(i).min_dis = min_dis;
                    S(i).min_dis_cluster = min_dis_cluster;
                end
            end
        end
        STATISTICS(r + 1).PACKETS_TO_CH = packets_TO_CH;

        su = 0;
        for i = 1:n
            if (S(i).E > 0)
                su = su + S(i).E;
            end
        end
        avg = su / n;
        STATISTICS(r + 1).AVG = avg;

        for i = 1:n
            if S(i).min_dis_cluster
                pt = C(S(i).min_dis_cluster).id;
            else
                pt = i;
            end
            if S(pt).count
                N = 400;
                M = 50;
                M = ceil(M);
                if N > M
                    rng(i);
                    Phi = randn(M, N);
                    S(i).compress_data = Phi * S(i).data';
                else
                    S(i).compress_data = S(i).data;
                end
                len = size(S(i).compress_data);
                S(pt).neigh_data(S(pt).count:S(pt).count + len - 1, 1) = S(i).compress_data;
                S(pt).neigh_data(S(pt).count + data_size - 1, 1) = i;
                S(pt).count = S(pt).count + data_size;
            end
        end

        count = 1;
        for i = 0:length(C) - 1
            clust_point = C(i + 1).id;
            data = S(clust_point).neigh_data;
            for i = 1:data_size:S(clust_point).count - 400
                node_id = data(399 + i, 1);
                sink.data(((node_id - 1) * data_size) + 1:(node_id * data_size), 2) = data(i:i + 399, 1);
            end
        end

        stoping_point = find(sink.data(:, 1) == -1);
        start = 1;
        for j = 1:n
            data = sink.data(start:start + M - 1, 2);
            node_id = sink.data(start + data_size - 1, 2);
            re_x = reconstruct(data, N, node_id);
            sink.data(start:start + data_size - 1, 3) = re_x;
            start = start + data_size;
        end
    end
end

function [STATISTICS, sink] = run_aco_clustering(S, n, rmax, p, xm, ym, ETX, ERX, Efs, Emp, EDA, do, data_size)
    % ACO parameters
    num_ants = n;
    alpha = 1; % Pheromone importance
    beta = 2; % Heuristic importance
    rho = 0.5; % Pheromone evaporation rate

    STATISTICS = struct();
    sink.x = xm / 2;
    sink.y = ym / 2;
    sink.data = zeros(data_size * n, 3);

    pheromone = ones(n, 1); % Initialize pheromone levels

    for r = 0:rmax
        if (mod(r, round(1 / p)) == 0)
            for i = 1:n
                S(i).G = 0;
                S(i).cl = 0;
            end
        end
        dead = 0;
        dead_a = 0;
        dead_n = 0;
        packets_TO_BS = 0;
        packets_TO_CH = 0;

        for i = 1:n
            if (S(i).E <= 0)
                dead = dead + 1;
                if (S(i).ENERGY == 1)
                    dead_a = dead_a + 1;
                else
                    dead_n = dead_n + 1;
                end
            end
        end
        STATISTICS(r + 1).DEAD = dead;
        STATISTICS(r + 1).DEAD_N = dead_n;
        STATISTICS(r + 1).DEAD_A = dead_a;

        countCHs = 0;
        cluster = 1;

        for i = 1:n
            if (S(i).E > 0)
                % ACO-based cluster head selection
                prob = (pheromone(i)^alpha) * ((1 / S(i).E)^beta);
                if rand < prob
                    countCHs = countCHs + 1;
                    packets_TO_BS = packets_TO_BS + 1;
                    STATISTICS(r + 1).PACKETS_TO_BS = packets_TO_BS;
                    S(i).type = 'C';
                    S(i).G = round(1 / p) - 1;
                    C(cluster).xd = S(i).xd;
                    C(cluster).yd = S(i).yd;
                    distance = sqrt((S(i).xd - sink.x)^2 + (S(i).yd - sink.y)^2);
                    C(cluster).distance = distance;
                    C(cluster).id = i;
                    cluster = cluster + 1;
                    if (distance > do)
                        S(i).E = S(i).E - ((ETX + EDA) * (4000) + Emp * 4000 * (distance^4));
                    else
                        S(i).E = S(i).E - ((ETX + EDA) * (4000) + Efs * 4000 * (distance^2));
                    end
                    pheromone(i) = (1 - rho) * pheromone(i) + rho * (1 / S(i).E);
                end
            end
        end
        STATISTICS(r + 1).CLUSTERHEADS = cluster - 1;

        for i = 1:n
            if (S(i).type == 'N' && S(i).E > 0)
                if (cluster - 1 >= 1)
                    min_dis = sqrt((S(i).xd - sink.x)^2 + (S(i).yd - sink.y)^2);
                    min_dis_cluster = 1;
                    for c = 1:cluster - 1
                        temp = min(min_dis, sqrt((S(i).xd - C(c).xd)^2 + (S(i).yd - C(c).yd)^2));
                        if (temp < min_dis)
                            min_dis = temp;
                            min_dis_cluster = c;
                        end
                    end
                    if (min_dis > do)
                        S(i).E = S(i).E - (ETX * (4000) + Emp * 4000 * (min_dis^4));
                    else
                        S(i).E = S(i).E - (ETX * (4000) + Efs * 4000 * (min_dis^2));
                    end
                    if (min_dis > 0)
                        distance = sqrt((S(C(min_dis_cluster).id).xd - sink.x)^2 + (S(C(min_dis_cluster).id).yd - sink.y)^2);
                        S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ERX + EDA) * 4000);
                        if (distance > do)
                            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ETX + EDA) * (4000) + Emp * 4000 * (distance^4));
                        else
                            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ETX + EDA) * (4000) + Efs * 4000 * (distance^2));
                        end
                        packets_TO_CH = n - dead - cluster + 1;
                    end
                    S(i).min_dis = min_dis;
                    S(i).min_dis_cluster = min_dis_cluster;
                end
            end
        end
        STATISTICS(r + 1).PACKETS_TO_CH = packets_TO_CH;

        su = 0;
        for i = 1:n
            if (S(i).E > 0)
                su = su + S(i).E;
            end
        end
        avg = su / n;
        STATISTICS(r + 1).AVG = avg;

        for i = 1:n
            if S(i).min_dis_cluster
                pt = C(S(i).min_dis_cluster).id;
            else
                pt = i;
            end
            if S(pt).count
                N = 400;
                M = 50;
                M = ceil(M);
                if N > M
                    rng(i);
                    Phi = randn(M, N);
                    S(i).compress_data = Phi * S(i).data';
                else
                    S(i).compress_data = S(i).data;
                end
                len = size(S(i).compress_data);
                S(pt).neigh_data(S(pt).count:S(pt).count + len - 1, 1) = S(i).compress_data;
                S(pt).neigh_data(S(pt).count + data_size - 1, 1) = i;
                S(pt).count = S(pt).count + data_size;
            end
        end

        count = 1;
        for i = 0:length(C) - 1
            clust_point = C(i + 1).id;
            data = S(clust_point).neigh_data;
            for i = 1:data_size:S(clust_point).count - 400
                node_id = data(399 + i, 1);
                sink.data(((node_id - 1) * data_size) + 1:(node_id * data_size), 2) = data(i:i + 399, 1);
            end
        end

        stoping_point = find(sink.data(:, 1) == -1);
        start = 1;
        for j = 1:n
            data = sink.data(start:start + M - 1, 2);
            node_id = sink.data(start + data_size - 1, 2);
            re_x = reconstruct(data, N, node_id);
            sink.data(start:start + data_size - 1, 3) = re_x;
            start = start + data_size;
        end
    end
end

function re_x = reconstruct(data, N, node_id)
    M = length(data);
    Phi = randn(M, N);
    re_x = pinv(Phi) * data;
end

function plot_comparisons(LEACH_STATISTICS, Q_STATISTICS, GAME_STATISTICS, PSO_STATISTICS, HEED_STATISTICS, ACO_STATISTICS, rmax)
    % Plot the number of dead nodes over rounds for all methods
    figure;
    plot([LEACH_STATISTICS.DEAD], 'r');
    hold on;
    plot([Q_STATISTICS.DEAD], 'b');
    plot([GAME_STATISTICS.DEAD], 'g');
    plot([PSO_STATISTICS.DEAD], 'm');
    plot([HEED_STATISTICS.DEAD], 'c');
    plot([ACO_STATISTICS.DEAD], 'k');
    title('Number of Dead Nodes Over Rounds');
    xlabel('Rounds');
    ylabel('Number of Dead Nodes');
    legend('LEACH', 'Q-learning', 'Game Theory', 'PSO', 'HEED', 'ACO');

    % Plot the number of cluster heads over rounds for all methods
    figure;
    plot([LEACH_STATISTICS.CLUSTERHEADS], 'r');
    hold on;
    plot([Q_STATISTICS.CLUSTERHEADS], 'b');
    plot([GAME_STATISTICS.CLUSTERHEADS], 'g');
    plot([PSO_STATISTICS.CLUSTERHEADS], 'm');
    plot([HEED_STATISTICS.CLUSTERHEADS], 'c');
    plot([ACO_STATISTICS.CLUSTERHEADS], 'k');
    title('Number of Cluster Heads Over Rounds');
    xlabel('Rounds');
    ylabel('Number of Cluster Heads');
    legend('LEACH', 'Q-learning', 'Game Theory', 'PSO', 'HEED', 'ACO');

    % Plot the energy consumption over rounds for all methods
    figure;
    plot([LEACH_STATISTICS.AVG], 'r');
    hold on;
    plot([Q_STATISTICS.AVG], 'b');
    plot([GAME_STATISTICS.AVG], 'g');
    plot([PSO_STATISTICS.AVG], 'm');
    plot([HEED_STATISTICS.AVG], 'c');
    plot([ACO_STATISTICS.AVG], 'k');
    title('Average Energy Consumption Over Rounds');
    xlabel('Rounds');
    ylabel('Average Energy');
    legend('LEACH', 'Q-learning', 'Game Theory', 'PSO', 'HEED', 'ACO');
end
