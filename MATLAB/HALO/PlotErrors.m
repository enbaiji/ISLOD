addpath('output');
load('ORBdataSST.mat');

posErrors_DRO    = orb.EstDRO.seq.a.posErrors;
velErrors_DRO    = orb.EstDRO.seq.a.velErrors;
posErrors_ELFO   = orb.EstELFO.seq.a.posErrors;
velErrors_ELFO   = orb.EstELFO.seq.a.velErrors;

BiasErrors_DRO   = orb.EstDRO.seq.a.clcErrors(:, 1);
DriftErrors_DRO  = orb.EstDRO.seq.a.clcErrors(:, 2);
BiasErrors_ELFO  = orb.EstELFO.seq.a.clcErrors(:, 1);
DriftErrors_ELFO = orb.EstELFO.seq.a.clcErrors(:, 2);

c = 299792458; % m/s

%% DRO
figure;
subplot(2,2,1);
semilogy((orb.EstDRO.seq.a.t - orb.EstDRO.seq.a.t(1)) / 3600, posErrors_DRO*1e3, 'LineWidth', 1.5);
grid on;
xlabel("Elapsed Time (hr)");
ylabel("Error (m)");
title("DRO Position Errors");
set(gca, 'MinorGridLineStyle', 'none'); 

subplot(2,2,2);
semilogy((orb.EstDRO.seq.a.t - orb.EstDRO.seq.a.t(1)) / 3600, velErrors_DRO*1e6, 'LineWidth', 1.5);
grid on;
xlabel("Elapsed Time (hr)");
ylabel("Error (mm/s)");
title("DRO Velocity Errors");
set(gca, 'MinorGridLineStyle', 'none'); 

subplot(2,2,3);
semilogy((orb.EstDRO.seq.a.t - orb.EstDRO.seq.a.t(1)) / 3600, BiasErrors_DRO * c, 'LineWidth', 1.5);
grid on;
xlabel("Elapsed Time (hr)");
ylabel("Error (m)");
title("DRO Clock Bias Errors");
set(gca, 'MinorGridLineStyle', 'none'); 

subplot(2,2,4);
semilogy((orb.EstDRO.seq.a.t - orb.EstDRO.seq.a.t(1)) / 3600, DriftErrors_DRO * c * 1e3, 'LineWidth', 1.5);
grid on;
xlabel("Elapsed Time (hr)");
ylabel("Error (mm/s)");
title("DRO Clock Drift Errors");
set(gca, 'MinorGridLineStyle', 'none'); 

%% ELFO
figure;
subplot(2,2,1);
semilogy((orb.EstELFO.seq.a.t - orb.EstELFO.seq.a.t(1)) / 3600, posErrors_ELFO*1e3, 'LineWidth', 1.5);
hold on;
xline(48, 'DisplayName', 'Arc Finish');
hold off;
grid on;
xlabel("Elapsed Time (hr)");
ylabel("Error (m)");
title("ELFO Position Errors");
set(gca, 'MinorGridLineStyle', 'none'); 

subplot(2,2,2);
semilogy((orb.EstELFO.seq.a.t - orb.EstELFO.seq.a.t(1)) / 3600, velErrors_ELFO*1e6, 'LineWidth', 1.5);
hold on;
xline(48, 'DisplayName', 'Arc Finish');
hold off;
grid on;
xlabel("Elapsed Time (hr)");
ylabel("Error (mm/s)");
title("ELFO Velocity Errors");
set(gca, 'MinorGridLineStyle', 'none'); 

subplot(2,2,3);
semilogy((orb.EstELFO.seq.a.t - orb.EstELFO.seq.a.t(1)) / 3600, BiasErrors_ELFO * c, 'LineWidth', 1.5);
hold on;
xline(48, 'DisplayName', 'Arc Finish');
hold off;
grid on;
xlabel("Elapsed Time (hr)");
ylabel("Error (m)");
title("ELFO Clock Bias Errors");
set(gca, 'MinorGridLineStyle', 'none');

subplot(2,2,4);
semilogy((orb.EstELFO.seq.a.t - orb.EstELFO.seq.a.t(1)) / 3600, DriftErrors_ELFO * c * 1e3, 'LineWidth', 1.5);
hold on;
xline(48, 'DisplayName', 'Arc Finish');
hold off;
grid on;
xlabel("Elapsed Time (hr)");
ylabel("Error (mm/s)");
title("ELFO Clock Drift Errors");
set(gca, 'MinorGridLineStyle', 'none');

%%
alpha = 0.05;

%% ACC NEES
DoF = 6 * ones(length(arrayNEES_DRO_acc), 1);
arrayNEES_DRO_acc(arrayNEES_DRO_acc==0) = NaN;                   
arrayNEES_DRO_acc = fillmissing(arrayNEES_DRO_acc,'previous');  
avg = mean(arrayNEES_DRO_acc(round(length(arrayNEES_DRO_acc)/2):end));
lower = chi2inv(alpha/2, N*6)/N * ones(length(arrayNEES_DRO_acc), 1);
upper = chi2inv(1 - alpha/2, N*6)/N * ones(length(arrayNEES_DRO_acc), 1);

figure;
subplot(3,2,1);
semilogy((orb.EstDRO.seq.a.t - orb.EstDRO.seq.a.t(1)) / 3600, arrayNEES_DRO_acc, ...
    'LineWidth', 1.5, 'DisplayName', ['NEES, C.avg=', num2str(avg, '%.2f')]);
hold on;
semilogy((orb.EstDRO.seq.a.t - orb.EstDRO.seq.a.t(1)) / 3600, DoF, ...
    'LineWidth', 1.5, 'DisplayName', 'DoF=6');
semilogy((orb.EstDRO.seq.a.t - orb.EstDRO.seq.a.t(1)) / 3600, lower, ...
    'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Lower Bound');
semilogy((orb.EstDRO.seq.a.t - orb.EstDRO.seq.a.t(1)) / 3600, upper, ...
    'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Upper Bound');
hold off;
legend('Location', 'best');
grid on;
xlabel("Elapsed Time (hr)");
ylabel("NEES");
title("DRO NEES (Orbit)");
set(gca, 'MinorGridLineStyle', 'none');
fprintf('DRO NEES (Orbit)=%.2f\n', avg);

arrayNEES_ELFO_acc(arrayNEES_ELFO_acc==0) = NaN;                   
arrayNEES_ELFO_acc = fillmissing(arrayNEES_ELFO_acc,'previous');  
avg = mean(arrayNEES_ELFO_acc(round(length(arrayNEES_ELFO_acc)/2):end));
lower = chi2inv(alpha/2, N*6)/N * ones(length(arrayNEES_ELFO_acc), 1);
upper = chi2inv(1 - alpha/2, N*6)/N * ones(length(arrayNEES_ELFO_acc), 1);

subplot(3,2,2);
semilogy((orb.EstELFO.seq.a.t-orb.EstELFO.seq.a.t(1))/3600, arrayNEES_ELFO_acc, ...
    'LineWidth', 1.5, 'DisplayName', ['NEES, C.avg=', num2str(avg, '%.2f')]);
hold on;
xline(48, 'DisplayName', 'Arc Finish');
semilogy((orb.EstELFO.seq.a.t - orb.EstELFO.seq.a.t(1)) / 3600, DoF, ...
    'LineWidth', 1.5, 'DisplayName', 'DoF=6');
semilogy((orb.EstELFO.seq.a.t - orb.EstELFO.seq.a.t(1)) / 3600, lower, ...
    'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Lower Bound');
semilogy((orb.EstELFO.seq.a.t - orb.EstELFO.seq.a.t(1)) / 3600, upper, ...
    'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Upper Bound');
hold off;
legend('Location', 'best');
grid on;
xlabel("Elapsed Time (hr)");
ylabel("NEES");
title("ELFO NEES (Orbit)");
set(gca, 'MinorGridLineStyle', 'none'); 
fprintf('ELFO NEES (Orbit)=%.2f\n', avg);

%% CLC NEES
DoF = 2 * ones(length(arrayNEES_DRO_clc), 1);
arrayNEES_DRO_clc(arrayNEES_DRO_clc==0) = NaN;                   
arrayNEES_DRO_clc = fillmissing(arrayNEES_DRO_clc,'previous');  
avg = mean(arrayNEES_DRO_clc(round(length(arrayNEES_DRO_clc)/2):end));
lower = chi2inv(alpha/2, N*2)/N * ones(length(arrayNEES_DRO_clc), 1);
upper = chi2inv(1 - alpha/2, N*2)/N * ones(length(arrayNEES_DRO_clc), 1);

subplot(3,2,3);
semilogy((orb.EstDRO.seq.a.t - orb.EstDRO.seq.a.t(1)) / 3600, arrayNEES_DRO_clc, ...
    'LineWidth', 1.5, 'DisplayName', ['NEES, C.avg=', num2str(avg, '%.2f')]);
hold on;
semilogy((orb.EstDRO.seq.a.t - orb.EstDRO.seq.a.t(1)) / 3600, DoF, ...
    'LineWidth', 1.5, 'DisplayName', 'DoF=2');
semilogy((orb.EstDRO.seq.a.t - orb.EstDRO.seq.a.t(1)) / 3600, lower, ...
    'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Lower Bound');
semilogy((orb.EstDRO.seq.a.t - orb.EstDRO.seq.a.t(1)) / 3600, upper, ...
    'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Upper Bound');
hold off;
legend('Location', 'best');
grid on;
xlabel("Elapsed Time (hr)");
ylabel("NEES");
title("DRO NEES (Clock)");
set(gca, 'MinorGridLineStyle', 'none');
fprintf('DRO NEES (Clock)=%.2f\n', avg);

arrayNEES_ELFO_clc(arrayNEES_ELFO_clc==0) = NaN;                  
arrayNEES_ELFO_clc = fillmissing(arrayNEES_ELFO_clc,'previous');  
avg = mean(arrayNEES_ELFO_clc(round(length(arrayNEES_ELFO_clc)/2):end));
lower = chi2inv(alpha/2, N*2)/N * ones(length(arrayNEES_ELFO_clc), 1);
upper = chi2inv(1 - alpha/2, N*2)/N * ones(length(arrayNEES_ELFO_clc), 1);

subplot(3,2,4);
semilogy((orb.EstELFO.seq.a.t-orb.EstELFO.seq.a.t(1))/3600, arrayNEES_ELFO_clc, ...
    'LineWidth', 1.5, 'DisplayName', ['NEES, C.avg=', num2str(avg, '%.2f')]);
hold on;
xline(48, 'DisplayName', 'Arc Finish');
semilogy((orb.EstELFO.seq.a.t - orb.EstELFO.seq.a.t(1)) / 3600, DoF, ...
    'LineWidth', 1.5, 'DisplayName', 'DoF=2');
semilogy((orb.EstELFO.seq.a.t - orb.EstELFO.seq.a.t(1)) / 3600, lower, ...
    'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Lower Bound');
semilogy((orb.EstELFO.seq.a.t - orb.EstELFO.seq.a.t(1)) / 3600, upper, ...
    'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Upper Bound');
hold off;
legend('Location', 'best');
grid on;
xlabel("Elapsed Time (hr)");
ylabel("NEES");
title("ELFO NEES (Clock)");
set(gca, 'MinorGridLineStyle', 'none');
fprintf('ELFO NEES (Clock)=%.2f\n', avg);

%% NIS
DoF = 4 * ones(length(arrayNIS_DRO), 1);
arrayNIS_DRO(arrayNIS_DRO==0) = NaN;                   
arrayNIS_DRO = fillmissing(arrayNIS_DRO,'previous');  
avg = mean(arrayNIS_DRO(round(length(arrayNIS_DRO)/2):end));
lower = chi2inv(alpha/2, N*4)/N * ones(length(arrayNIS_DRO), 1);
upper = chi2inv(1 - alpha/2, N*4)/N * ones(length(arrayNIS_DRO), 1);

subplot(3,2,5);
semilogy((orb.EstDRO.seq.a.t - orb.EstDRO.seq.a.t(1)) / 3600, arrayNIS_DRO, ...
    'LineWidth', 1.5, 'DisplayName', ['NIS, C.avg=', num2str(avg, '%.2f')]);
hold on;
semilogy((orb.EstDRO.seq.a.t - orb.EstDRO.seq.a.t(1)) / 3600, DoF, ...
    'LineWidth', 1.5, 'DisplayName', 'DoF=4');
semilogy((orb.EstDRO.seq.a.t - orb.EstDRO.seq.a.t(1)) / 3600, lower, ...
    'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Lower Bound');
semilogy((orb.EstDRO.seq.a.t - orb.EstDRO.seq.a.t(1)) / 3600, upper, ...
    'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Upper Bound');
hold off;
legend('Location', 'best');
grid on;
xlabel("Elapsed Time (hr)");
ylabel("NIS");
title("DRO NIS");
set(gca, 'MinorGridLineStyle', 'none');
fprintf('DRO NIS=%.2f\n', avg);

DoF1 = 4 * ones(floor(length(arrayNIS_ELFO)/5*2), 1);
DoF2 = 2 * ones(floor(length(arrayNIS_ELFO)/5*3)+1, 1);
DoF = [DoF1; DoF2];
arrayNIS_ELFO(arrayNIS_ELFO==0) = NaN;                  
arrayNIS_ELFO = fillmissing(arrayNIS_ELFO,'previous');  
avg = mean(arrayNIS_ELFO(round(length(arrayNIS_ELFO)/2):end));
lower1 = chi2inv(alpha/2, N*4)/N * ones(floor(length(arrayNIS_ELFO)/5*2), 1);
lower2 = chi2inv(alpha/2, N*2)/N * ones(floor(length(arrayNIS_ELFO)/5*3)+1, 1);
upper1 = chi2inv(1 - alpha/2, N*4)/N * ones(floor(length(arrayNIS_ELFO)/5*2), 1);
upper2 = chi2inv(1 - alpha/2, N*2)/N * ones(floor(length(arrayNIS_ELFO)/5*3)+1, 1);
lower = [lower1; lower2];
upper = [upper1; upper2];

subplot(3,2,6);
semilogy((orb.EstELFO.seq.a.t-orb.EstELFO.seq.a.t(1))/3600, arrayNIS_ELFO, ...
    'LineWidth', 1.5, 'DisplayName', ['NIS, C.avg=', num2str(avg, '%.2f')]);
hold on;
xline(48, 'DisplayName', 'Arc Finish');
semilogy((orb.EstELFO.seq.a.t - orb.EstELFO.seq.a.t(1)) / 3600, DoF, ...
    'LineWidth', 1.5, 'DisplayName', 'DoF=4 or 2');
semilogy((orb.EstELFO.seq.a.t - orb.EstELFO.seq.a.t(1)) / 3600, lower, ...
    'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Lower Bound');
semilogy((orb.EstELFO.seq.a.t - orb.EstELFO.seq.a.t(1)) / 3600, upper, ...
    'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Upper Bound');
hold off;
legend('Location', 'best');
grid on;
xlabel("Elapsed Time (hr)");
ylabel("NIS");
title("ELFO NIS");
set(gca, 'MinorGridLineStyle', 'none');
fprintf('ELFO NIS=%.2f\n', avg);


