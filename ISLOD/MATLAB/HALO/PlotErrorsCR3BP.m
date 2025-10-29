nRelay = 2;

hrs = 300;                             % simulation span in hrs
dt_p = 0.000157 / 5;                   % propagation step (0.2 min)

if nRelay == 1
    load('output/dro1(1)');
    load('output/elfo(1)');
elseif nRelay == 2
    load('output/dro1(2)');
    load('output/dro2(2)');
    load('output/elfo(2)');
elseif nRelay == 3
    load('output/dro1(3)');
    load('output/dro2(3)');
    load('output/dro3(3)');
    load('output/elfo(3)');
end

% % % DRO
% figure;
% subplot(2, 2, 1);
% plot_position_errors(dro2.errors, dt_p, hrs, 'Ground Meausred DRO 2', 51);
% subplot(2, 2, 2);
% plot_velocity_errors(dro2.errors, dt_p, hrs, 'Ground Measured DRO 2', 43);
% subplot(2, 2, 3);
% plot_rErrorsNorm(dro1.rErrorsNorm, dro2.rErrorsNorm, [], dt_p, hrs, 49, 51, []);
% subplot(2, 2, 4);
% plot_vErrorsNorm(dro1.vErrorsNorm, dro2.vErrorsNorm, [], dt_p, hrs, 40, 43, [])
% sgtitle('Ground-Based OD Errors of DRO Satellite');
% 
% % % % ELFO
% % figure;
% % subplot(2, 2, 1);
% % plot_position_errors(elfo.errors, dt_p, hrs, 'Ground Meausred ELFO', 59);
% % subplot(2, 2, 2);
% % plot_velocity_errors(elfo.errors, dt_p, hrs, 'Ground Measured ELFO', 35);
% % subplot(2, 2, 3);
% % plot_rErrorsNorm(elfo.rErrorsNorm, [], [], dt_p, hrs, 59);
% % subplot(2, 2, 4);
% % plot_vErrorsNorm(elfo.vErrorsNorm, [], [], dt_p, hrs, 35);
% % sgtitle('Ground-Based OD Errors of ELFO Satellite');
% 
% % ELFO of SST
% figure;
% subplot(2, 2, 1);
% plot_position_errors(elfo.errors_ISL, dt_p, hrs, 'SST Measured ELFO', 162);
% subplot(2, 2, 2);
% plot_velocity_errors(elfo.errors_ISL, dt_p, hrs, 'SST Measured ELFO', 171);
% subplot(2, 2, 3);
% plot_rErrorsNorm(elfo.rErrorsNorm_ISL, [], [], dt_p, hrs, 162);
% subplot(2, 2, 4);
% plot_vErrorsNorm(elfo.vErrorsNorm_ISL, [], [], dt_p, hrs, 171);
% sgtitle('One-Way SST-Based OD Errors of ELFO Satellite');
% 
% % Visibility
% % figure;
% % subplot(3,1,1);
% % plot_SiteVisib(dro1.visArrays(:,:,1), dt_p, hrs, 'DRO 1');
% % subplot(3,1,2);
% % plot_SiteVisib(elfo.visArrays(:,:,1), dt_p, hrs, 'ELFO');
% % subplot(3,1,3);
% % plot_ISLVisib(elfo.ISLvisArrays(:,:,1), dt_p, hrs);
% 
% % Ellipsoid
% figure;
% subplot(3,2,1);
% plot_PositionEllipsoid(dro2.Ellipsoid_r, [], [], dt_p, hrs, 'Ground Measured DRO 2', 51);
% subplot(3,2,2);
% plot_VelocityEllipsoid(dro2.Ellipsoid_v, [], [], dt_p, hrs, 'Ground Measured DRO 2', 43);
% subplot(3,2,3);
% plot_PositionEllipsoid(elfo.Ellipsoid_r, [], [], dt_p, hrs, 'Ground Measured ELFO', 59);
% subplot(3,2,4);
% plot_VelocityEllipsoid(elfo.Ellipsoid_v, [], [], dt_p, hrs, 'Ground Measured ELFO', 35);
% subplot(3,2,5);
% plot_PositionEllipsoid(elfo.ISL_Ellipsoid_r, [], [], dt_p, hrs, 'SST Measured ELFO', 162);
% subplot(3,2,6);
% plot_VelocityEllipsoid(elfo.ISL_Ellipsoid_v, [], [], dt_p, hrs, 'SST Measured ELFO', 171);

% NEES and NIS
% figure;
% subplot(2,3,1);
% plot_NEES(dro3.NEES, dt_p, kRuns, 6, hrs, 'Ground-Measured DRO 3',100);
% subplot(2,3,2);
% plot_NEES(elfo.NEES, dt_p, kRuns, 6, hrs, 'Ground-Measured ELFO',100);
% subplot(2,3,3);
% plot_NEES(elfo.NEES_sst, dt_p, kRuns, 6, hrs, 'SST-Measured ELFO',100);
% subplot(2,3,4);
% plot_NIS(dro3.NIS, dt_p, kRuns, 4, hrs, 'Ground-Measured DRO 3',100);
% subplot(2,3,5);
% plot_NIS(elfo.NIS, dt_p, kRuns, 4, hrs, 'Ground-Measured ELFO',100);
% subplot(2,3,6);
% plot_NIS(elfo.NIS_sst, dt_p, kRuns, 6, hrs, 'SST-Measured ELFO',100);


elfo1 = load('output/elfo(1)');
elfo2 = load('output/elfo(2)');
elfo3 = load('output/elfo(3)');

figure;
subplot(2,2,1)
compare3casesPos(elfo1.elfo.rErrorsNorm_ISL, elfo2.elfo.rErrorsNorm_ISL, elfo3.elfo.rErrorsNorm_ISL, dt_p);
subplot(2,2,2)
compare3casesVel(elfo1.elfo.vErrorsNorm_ISL, elfo2.elfo.vErrorsNorm_ISL, elfo3.elfo.vErrorsNorm_ISL, dt_p);

function plot_position_errors(array, dt_p, xlimit, string, startHr)

    % Time vector calculation (scaled by 4.43)
    time_vector = (0:dt_p:(length(array) - 1) * dt_p) * 4.43*24;
    hold on; grid on;

    n = size(array, 1);
    start_index = round(startHr/300*n) + 1;

    last_elements = array(start_index:n, 1);
    average = mean(last_elements) * 389703 * 1000;
    disp(['Average x position = ', num2str(average), ' ', string])

    last_elements = array(start_index:n, 2);
    average = mean(last_elements) * 389703 * 1000;
    disp(['Average y position = ', num2str(average), ' ', string])

    last_elements = array(start_index:n, 3);
    average = mean(last_elements) * 389703 * 1000;
    disp(['Average z position = ', num2str(average), ' ', string])

    % Plot x, y, z position errors
    semilogy(time_vector, abs(array(:, 1)) * 389703, 'b', 'LineWidth', 0.5, 'DisplayName', 'x', 'LineWidth', 1.5);
    semilogy(time_vector, abs(array(:, 2)) * 389703, 'r', 'LineWidth', 0.5, 'DisplayName', 'y', 'LineWidth', 1.5);
    semilogy(time_vector, abs(array(:, 3)) * 389703, 'g', 'LineWidth', 0.5, 'DisplayName', 'z', 'LineWidth', 1.5);

    set(gca, 'MinorGridLineStyle', 'none');
    xlim([0 xlimit]);
    xlabel('Time (hr)');
    ylabel('Position Error (km)');
    title('Position OD Errors in x, y and z Directions');
    legend('Location', 'best');
    hold off;
    ax = gca; 
    set(ax,'YScale','log'); 
    ytickformat(ax,'%.0e'); 
end

function plot_velocity_errors(array, dt_p, xlimit, string, startHr)

    % Compute the time vector (scaled by 4.43 to convert to days)
    time_vector = (0:dt_p:(length(array) - 1) * dt_p) * 4.43*24;
    hold on; grid on;

    n = size(array, 1);
    start_index = round(startHr/300*n) + 1;

    last_elements = array(start_index:n, 4);
    average = mean(last_elements) * 1.01755178455328 * 1e6;
    disp(['Average x velocity = ', num2str(average), ' ', string])

    last_elements = array(start_index:n, 5);
    average = mean(last_elements) * 1.01755178455328 * 1e6;
    disp(['Average y velocity = ', num2str(average), ' ', string])

    last_elements = array(start_index:n, 6);
    average = mean(last_elements) * 1.01755178455328 * 1e6;
    disp(['Average z velocity = ', num2str(average), ' ', string])

    % Plot velocity errors in x, y, z directions
    semilogy(time_vector, abs(array(:, 4)) * 1.01755178455328 * 1e3, 'b', 'DisplayName', 'x', 'LineWidth', 1.5);
    semilogy(time_vector, abs(array(:, 5)) * 1.01755178455328 * 1e3, 'r', 'DisplayName', 'y', 'LineWidth', 1.5);
    semilogy(time_vector, abs(array(:, 6)) * 1.01755178455328 * 1e3, 'g', 'DisplayName', 'z', 'LineWidth', 1.5);

    set(gca, 'MinorGridLineStyle', 'none');
    xlim([0 xlimit]);
    xlabel('Time (hr)');
    ylabel('Velocity Error (m/s)');
    title('Velocity OD Errors in x, y and z Directions');
    legend('Location', 'best');
    hold off;
    ax = gca;
    set(ax,'YScale','log');
    ytickformat(ax,'%.0e');
end

function plot_rErrorsNorm(array1, array2, array3, dt_p, xlimit, startHr1, startHr2, startHr3)

    % Compute the time vector (scaled by 4.43 to represent time in days)
    time_vector = (0:dt_p:(length(array1) - 1) * dt_p) * 4.43*24;
    hold on; grid on;

    n = size(array1);

    if isempty(array2) && isempty(array3)

        start_index = round(startHr1/300*n) + 1;
        last_elements = array1(start_index:n);
        average1 = mean(last_elements)*1000;

        semilogy(time_vector, array1, 'b', 'DisplayName', ['C.Avg.=', num2str(average1, '%.1f'), 'm'], 'LineWidth', 1.5);

    elseif isempty(array3)

        start_index = round(startHr1/300*n) + 1;
        last_elements = array1(start_index:n);
        average1 = mean(last_elements)*1000;

        start_index = round(startHr2/300*n) + 1;
        last_elements = array2(start_index:n);
        average2 = mean(last_elements)*1000;

        semilogy(time_vector, array1, 'b', 'DisplayName', ['DRO1, C.Avg.=', num2str(average1, '%.1f'), 'm'], 'LineWidth', 1.5);
        semilogy(time_vector, array2, 'r', 'DisplayName', ['DRO2, C.Avg.=', num2str(average2, '%.1f'), 'm'], 'LineWidth', 1.5);

    else

        start_index = round(startHr1/300*n) + 1;
        last_elements = array1(start_index:n);
        average1 = mean(last_elements)*1000;

        start_index = round(startHr2/300*n) + 1;
        last_elements = array2(start_index:n);
        average2 = mean(last_elements)*1000;

        start_index = round(startHr3/300*n) + 1;
        last_elements = array3(start_index:n);
        average3 = mean(last_elements)*1000;

        semilogy(time_vector, array1, 'b', 'DisplayName', ['DRO1, C.Avg.=', num2str(average1, '%.1f'), 'm'], 'LineWidth', 1.5);
        semilogy(time_vector, array2, 'r', 'DisplayName', ['DRO2, C.Avg.=', num2str(average2, '%.1f'), 'm'], 'LineWidth', 1.5);
        semilogy(time_vector, array3, 'g', 'DisplayName', ['DRO3, C.Avg.=', num2str(average3, '%.1f'), 'm'], 'LineWidth', 1.5);
    end

    set(gca, 'MinorGridLineStyle', 'none');
    xlim([0 xlimit]);
    xlabel('Time (hr)');
    ylabel('Position Error (km)');
    title('Position OD Error in Norm');
    legend('Location', 'best');
    hold off;
    ax = gca;
    set(ax,'YScale','log');
    ytickformat(ax,'%.0e');
end

function plot_vErrorsNorm(array1, array2, array3, dt_p, xlimit, startHr1, startHr2, startHr3)

    % Compute the time vector (scaled by 4.43 to convert to days)
    time_vector = (0:dt_p:(length(array1) - 1) * dt_p) * 4.43*24;
    hold on; grid on;

    n = size(array1);

    if isempty(array2) && isempty(array3)

        start_index = round(startHr1/300*n) + 1;
        last_elements = array1(start_index:n);
        average1 = mean(last_elements)*1000;

        semilogy(time_vector, array1, 'b', 'DisplayName', ['C.Avg.=', num2str(average1, '%.2f'), 'mm/s'], 'LineWidth', 1.5);

    elseif isempty(array3)

        start_index = round(startHr1/300*n) + 1;
        last_elements = array1(start_index:n);
        average1 = mean(last_elements)*1000;

        start_index = round(startHr2/300*n) + 1;
        last_elements = array2(start_index:n);
        average2 = mean(last_elements)*1000;

        semilogy(time_vector, array1, 'b', 'DisplayName', ['DRO1, C.Avg.=', num2str(average1, '%.2f'), 'mm/s'], 'LineWidth', 1.5);
        semilogy(time_vector, array2, 'r', 'DisplayName', ['DRO2, C.Avg.=', num2str(average2, '%.2f'), 'mm/s'], 'LineWidth', 1.5);

    else

        start_index = round(startHr1/300*n) + 1;
        last_elements = array1(start_index:n);
        average1 = mean(last_elements)*1000;

        start_index = round(startHr2/300*n) + 1;
        last_elements = array2(start_index:n);
        average2 = mean(last_elements)*1000;

        start_index = round(startHr3/300*n) + 1;
        last_elements = array3(start_index:n);
        average3 = mean(last_elements)*1000;

        semilogy(time_vector, array1, 'b', 'DisplayName', ['DRO1, C.Avg.=', num2str(average1, '%.2f'), 'mm/s'], 'LineWidth', 1.5);
        semilogy(time_vector, array2, 'r', 'DisplayName', ['DRO2, C.Avg.=', num2str(average2, '%.2f'), 'mm/s'], 'LineWidth', 1.5);
        semilogy(time_vector, array3, 'g', 'DisplayName', ['DRO3, C.Avg.=', num2str(average3, '%.2f'), 'mm/s'], 'LineWidth', 1.5);
    end

    set(gca, 'MinorGridLineStyle', 'none');
    xlim([0 xlimit]);
    xlabel('Time (hr)');
    ylabel('Velocity Error (m/s)')
    title('Velocity OD Error in Norm');
    legend('Location', 'best');
    hold off;
    ax = gca; 
    set(ax,'YScale','log');
    ytickformat(ax,'%.0e');
end

function plot_SiteVisib(visArray, dt_p, xlimit, string)

    % Compute the time vector (scaled by 4.43 to convert to days)
    time_vector = (0:dt_p:(length(visArray) - 1) * dt_p) * 4.43*24;

    plot(time_vector, visArray, 'k', 'LineWidth', 1.5);
    ylim([0 4]);
    xlim([0 xlimit]);
    ylabel('Selected Site #');
    xlabel('Time (hr)');
    title(['Site Selection for ', string, ' Satellite in Ground-Based Measurements']);
    grid on;
end

function plot_ISLVisib(visArray, dt_p, xlimit)

    % Compute the time vector (scaled by 4.43 to convert to days)
    time_vector = (0:dt_p:(length(visArray) - 1) * dt_p) * 4.43*24;

    plot(time_vector, visArray, 'k', 'LineWidth', 1.5);
    ylim([0 max(visArray)+1]);
    xlim([0 xlimit]);
    ylabel('Relay Availability');
    xlabel('Time (hr)');
    title('Number of Available Relays for ELFO Satellite in SST Measurements')
    grid on;
end

function plot_PositionEllipsoid(array1, array2, array3, dt_p, xlimit, string, startHr)

    % Compute the time vector (scaled by 4.43 to convert to days)
    time_vector = (0:dt_p:(length(array1) - 1) * dt_p) * 4.43*24;
    hold on; grid on;

    n = size(array1);
    start_index = round(startHr/300*n) + 1;
    last_elements = array1(start_index:n);
    average = mean(last_elements);
    disp(['Average Position Ellipsoid = ', num2str(average), ' ' string]);

    if isempty(array2) && isempty(array3)
        semilogy(time_vector, array1, 'b', 'LineWidth', 1.5);
    else
        semilogy(time_vector, array1, 'b', 'LineWidth', 1.5,'DisplayName','DRO1');
        semilogy(time_vector, array2, 'g', 'LineWidth', 1.5,'DisplayName','DRO2');
        semilogy(time_vector, array3, 'r', 'LineWidth', 1.5,'DisplayName','DRO3');
        legend('Location','best');
    end

    xlim([0 xlimit]);
    ylabel('Volume (m^3)');
    xlabel('Time (hr)');
    hold off;
    set(gca, 'MinorGridLineStyle', 'none');
    ax = gca;
    set(ax,'YScale','log'); 
    ytickformat(ax,'%.0e');
    title(['1-Sigma Position Covariance Ellipsoid', newline, 'Volume of ', string])
end

function plot_VelocityEllipsoid(array1, array2, array3, dt_p, xlimit, string, startHr)

    % Compute the time vector (scaled by 4.43 to convert to days)
    time_vector = (0:dt_p:(length(array1) - 1) * dt_p) * 4.43*24;
    hold on; grid on;

    n = size(array1);
    start_index = round(startHr/300*n) + 1;
    last_elements = array1(start_index:n);
    average = mean(last_elements);
    disp(['Average Velocity Ellipsoid = ', num2str(average), ' ' string]);

    if isempty(array2) && isempty(array3)
        semilogy(time_vector, array1, 'b', 'LineWidth', 1.5);
    else
        semilogy(time_vector, array1, 'b', 'LineWidth', 1.5,'DisplayName','DRO1');
        semilogy(time_vector, array2, 'g', 'LineWidth', 1.5,'DisplayName','DRO2');
        semilogy(time_vector, array3, 'r', 'LineWidth', 1.5,'DisplayName','DRO3');
        legend('Location','best');
    end

    xlim([0 xlimit]);
    ylabel('Volume (mm^3/s^3)');
    xlabel('Time (hr)');
    hold off;
    set(gca, 'MinorGridLineStyle', 'none'); 
    ax = gca;    
    set(ax,'YScale','log');    
    ytickformat(ax,'%.0e');       
    title(['1-Sigma Velocity Covariance Ellipsoid', newline, 'Volume of ', string])
end

function plot_NEES(NEES, dt_p, N, k, xlimit, string, startHr)

    time_vector = (0:dt_p:(length(NEES) - 1) * dt_p) * 4.43*24;

    n = size(NEES, 1);
    start_index = round(startHr/300*n) + 1;
    last_elements = NEES(start_index:n);
    average = mean(last_elements);

    alpha = 0.05;
    CI_avg_NEES_lower = chi2inv(alpha/2, N*k)/N * ones(length(NEES), 1);
    CI_avg_NEES_upper = chi2inv(1 - alpha/2, N*k)/N * ones(length(NEES), 1);
    CI_avg_NEES_mid   = k * ones(length(NEES), 1);

    hold on; grid on;

    semilogy(time_vector, NEES, 'DisplayName', ['NEES, C.Avg=', num2str(average, '%.2f')], 'LineWidth', 1);
    semilogy(time_vector, CI_avg_NEES_upper, 'DisplayName', 'Upper bound', 'LineWidth', 1, 'LineStyle','--');
    semilogy(time_vector, CI_avg_NEES_lower, 'DisplayName', 'Lower bound', 'LineWidth', 1, 'LineStyle','--');
    semilogy(time_vector, CI_avg_NEES_mid, ...
         'DisplayName', ['DoF=', num2str(k)], ...
         'LineStyle', '--', ...
         'Color', 'k', ...
         'LineWidth', 1);

    xlim([0 xlimit]);
    ylabel('NEES');
    xlabel('Time (hr)');
    set(gca, 'MinorGridLineStyle', 'none');
    ax = gca;
    set(ax,'YScale','log'); 
    ytickformat(ax,'%.0e'); 
    title(['NEES of ', string]);
    legend('Location', 'best');
    hold off;

    CI_avg_NEES_lower = chi2inv(alpha/2, N*k)/N * ones(length(last_elements), 1);
    CI_avg_NEES_upper = chi2inv(1 - alpha/2, N*k)/N * ones(length(last_elements), 1);
    check_NEES_consistency(last_elements, CI_avg_NEES_lower, CI_avg_NEES_upper, string, average)
end

function check_NEES_consistency(NEES_array, lower_bound, upper_bound, string, average)

    within_bounds = (NEES_array >= lower_bound) & (NEES_array <= upper_bound);

    percent_in_bound = sum(within_bounds) / numel(NEES_array) * 100;

    fprintf('%s: NEES Chi-Square Percentage: %.2f%%, ', string, percent_in_bound);
    fprintf('Average%.3f, ', average);
    fprintf('95%% Chi-Square Interval: [%.3f, %.3f], ', lower_bound(1), upper_bound(1));

    outliers = find(~within_bounds);
    fprintf('Number of Epochs out of Interval: %d(Total %d) \n', length(outliers), length(NEES_array));
end

function plot_NIS(NIS, dt_p, N, k, xlimit, string, startHr)

    time_vector = (0:dt_p:(length(NIS) - 1) * dt_p) * 4.43*24;

    n = size(NIS, 1);
    start_index = round(startHr/300*n) + 1;
    last_elements = NIS(start_index:n);
    average = mean(last_elements);

    alpha = 0.05;
    CI_avg_NIS_lower = chi2inv(alpha/2, N*k)/N * ones(length(NIS), 1);
    CI_avg_NIS_upper = chi2inv(1 - alpha/2, N*k)/N * ones(length(NIS), 1);
    CI_avg_NIS_mid   = k * ones(length(NIS), 1); 

    hold on; grid on;

    semilogy(time_vector, NIS, 'DisplayName', ['NIS, C.Avg=', num2str(average, '%.2f')], 'LineWidth', 1);
    semilogy(time_vector, CI_avg_NIS_upper, 'DisplayName', 'Upper bound', 'LineWidth', 1, 'LineStyle','--');
    semilogy(time_vector, CI_avg_NIS_lower, 'DisplayName', 'Lower bound', 'LineWidth', 1, 'LineStyle','--');
    semilogy(time_vector, CI_avg_NIS_mid, ...
         'DisplayName', ['DoF=', num2str(k)], ...
         'LineStyle', '--', ...
         'Color', 'k', ...
         'LineWidth', 1);

    xlim([0 xlimit]);
    ylabel('NIS');
    xlabel('Time (hr)');
    set(gca, 'MinorGridLineStyle', 'none'); 
    ax = gca;                 
    set(ax,'YScale','log');     
    ytickformat(ax,'%.0e');     
    title(['NIS of ', string]);
    legend('Location', 'best');
    hold off;

    CI_avg_NIS_lower = chi2inv(alpha/2, N*k)/N * ones(length(last_elements), 1);
    CI_avg_NIS_upper = chi2inv(1 - alpha/2, N*k)/N * ones(length(last_elements), 1);
    check_NIS_consistency(last_elements, CI_avg_NIS_lower, CI_avg_NIS_upper, string, average)
end

function check_NIS_consistency(NIS_array, lower_bound, upper_bound, string, average)

    within_bounds = (NIS_array >= lower_bound) & (NIS_array <= upper_bound);

    percent_in_bound = sum(within_bounds) / numel(NIS_array) * 100;

    fprintf('%s: NEES Chi-Square Percentage: %.2f%%, ', string, percent_in_bound);
    fprintf('Average%.3f, ', average);
    fprintf('95%% Chi-Square Interval: [%.3f, %.3f], ', lower_bound(1), upper_bound(1));

    outliers = find(~within_bounds);
    fprintf('Number of Epochs out of Interval: %d(Total %d) \n', length(outliers), length(NIS_array));
end

function compare3casesPos(array1, array2, array3, dt_p)

    time_vector = (0:dt_p:(length(array1) - 1) * dt_p) * 4.43*24;
    hold on; grid on;

    array1 = movmean(array1, 10000);
    array2 = movmean(array2, 10000);
    array3 = movmean(array3, 10000);

    semilogy(time_vector, array1, 'DisplayName', '1+1', 'LineWidth', 1.5);
    semilogy(time_vector, array2, 'DisplayName', '2+1', 'LineWidth', 1.5);
    semilogy(time_vector, array3, 'DisplayName', '3+1', 'LineWidth', 1.5);

    set(gca, 'MinorGridLineStyle', 'none'); 
    xlabel('Time (hr)');
    ylabel('Position Error (km)');
    legend('Location', 'best');
    hold off;
    ax = gca; 
    set(ax,'YScale','log'); 
    ytickformat(ax,'%.0e');
end

function compare3casesVel(array1, array2, array3, dt_p)

    time_vector = (0:dt_p:(length(array1) - 1) * dt_p) * 4.43*24;
    hold on; grid on;

    array1 = movmean(array1, 10000);
    array2 = movmean(array2, 10000);
    array3 = movmean(array3, 10000);

    semilogy(time_vector, array1, 'DisplayName', '1+1', 'LineWidth', 1.5);
    semilogy(time_vector, array2, 'DisplayName', '2+1', 'LineWidth', 1.5);
    semilogy(time_vector, array3, 'DisplayName', '3+1', 'LineWidth', 1.5);

    set(gca, 'MinorGridLineStyle', 'none');
    xlabel('Time (hr)');
    ylabel('Velocity Error (m/s)');
    legend('Location', 'best');
    hold off;
    ax = gca; 
    set(ax,'YScale','log');
    ytickformat(ax,'%.0e');
end





