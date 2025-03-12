clc
clear all;
close all;
xm = 300;
ym = 300;
sink.x = 100; % Adjust sink coordinates
sink.y = 75;

n = 10;
p = 0.1;
Eo = 0.5;
ETX = 50 * 0.000000001;
ERX = 50 * 0.000000001;
Efs = 10e-12;
Emp = 0.0013e-12;
EDA = 5 * 0.000000001;
rmax = 1000;
do = sqrt(Efs / Emp);
Et = 0;

% Q-value parameters
alpha = 0.5;
beta = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Q-value parameters
for i = 1:n
    Q(i).qi = rand(); % Initial Q-value for each node
    Q(i).last_transmission_round = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Actual Q-value based clustering method
for r = 0:rmax
    % Initialize network and nodes
    S(n+1).xd = sink.x;
    S(n+1).yd = sink.y;
    Et = 0;

    for i = 1:n
        S(i).xd = rand(1, 1) * xm;
        XR(i) = S(i).xd;
        S(i).yd = rand(1, 1) * ym;
        YR(i) = S(i).yd;
        distance = sqrt((S(i).xd - S(n+1).xd)^2 + (S(i).yd - S(n+1).yd)^2);
        S(i).distance = distance;
        S(i).G = 0;
        S(i).type = 'N';
        S(i).E = Eo;
        Et = Et + S(i).E;
        figure(1)
        plot(S(i).xd, S(i).yd, 'bo');
        text(S(i).xd + 1, S(i).yd - 0.5, num2str(i));
        hold on;
    end

    plot(S(n+1).xd, S(n+1).yd, 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    text(S(n+1).xd + 1, S(n+1).yd - 0.5, num2str(n+1));
    hold off;

    % Q-value based clustering
    countCHs = 0;
    cluster = 1;

    for i = 1:n
        if (S(i).E > 0)
            temp_rand = rand;
            if ((S(i).G) <= 0)
                if (temp_rand <= Q(i).qi)
                    countCHs = countCHs + 1;
                    S(i).type = 'C';
                    S(i).G = round(1 / p) - 1;

                    % Q-value update
                    Q(i).qi = alpha * Q(i).qi + beta * (1 - Q(i).qi);
                    
                    % Energy consumption (similar to LEACH)
                    distance = sqrt((S(i).xd - S(n+1).xd)^2 + (S(i).yd - S(n+1).yd)^2);
                    if (distance > do)
                        S(i).E = S(i).E - ((ETX + EDA) * (4000) + Emp * 4000 * (distance^4));
                    end
                    if (distance <= do)
                        S(i).E = S(i).E - ((ETX + EDA) * (4000) + Efs * 4000 * (distance^2));
                    end

                    % Cluster information
                    C(cluster).xd = S(i).xd;
                    C(cluster).yd = S(i).yd;
                    C(cluster).distance = distance;
                    C(cluster).id = i;
                    X(cluster) = S(i).xd;
                    Y(cluster) = S(i).yd;
                    cluster = cluster + 1;
                end
            end
        end
    end

    % Q-value update for normal nodes
    for i = 1:n
        if (S(i).type == 'N' && S(i).E > 0)
            Q(i).qi = alpha * Q(i).qi + beta * (1 - Q(i).qi);
        end
    end
    
    % Number of dead nodes
    dead = sum([S.E] <= 0);

    % Statistics collection (similar to LEACH)
    STATISTICS.DEAD(r+1) = dead;
    STATISTICS.ALLIVE(r+1) = n - dead;
    STATISTICS.PACKETS_TO_BS(r+1) = countCHs;
    STATISTICS.PACKETS_TO_BS_PER_ROUND(r+1) = countCHs;
    STATISTICS.PACKETS_TO_CH(r+1) = countCHs;
    
    % Throughput calculation
    STATISTICS.THROUGHPUT(r+1) = STATISTICS.PACKETS_TO_BS(r+1) + STATISTICS.PACKETS_TO_CH(r+1);
    
    % Average Residual Energy
    En = sum([S.E]);
    STATISTICS.ENERGY(r+1) = En;

    % Plotting
    figure(2)
    plot(S(n+1).xd, S(n+1).yd, 'x');
    hold on;
    plot(X, Y, 'k*');
    axis([0 xm 0 ym]);
    title(['Q-value Clustering (Round ', num2str(r), ')']);
    hold off;
end

% Plotting statistics
r = 0:rmax;
figure(3)
plot(r, STATISTICS.DEAD);
title('Dead Nodes')

figure(4)
plot(r, STATISTICS.ALLIVE);
title('Live Nodes')

figure(5)
plot(r, STATISTICS.PACKETS_TO_BS);
title('Packets to BS')

figure(6)
plot(r, STATISTICS.PACKETS_TO_BS_PER_ROUND);
title('Packets to BS per Round')

figure(7)
plot(r, STATISTICS.PACKETS_TO_CH);
title('Packets to CH')

figure(8)
plot(r, STATISTICS.THROUGHPUT);
title('Throughput')

figure(9)
plot(r, STATISTICS.ENERGY);
title('Average Residual Energy')
