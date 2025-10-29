tic;
clear;

%% ------------ Script Setting -------------
nRelay = 3;                            % Number of relays（1，2，3）
R_GB_factor = 1;                       % Ground measurement noise factor
R_sst_factor = 0.1;                    % ISL measurement nosie factor for filter

%% ------------ Constants -------------
mu = 1.215058560962404E-2;             % Mass ratio Earth‑Moon
LU = 389703;                           % Length unit [km]
TU = 382981;                           % Time unit  [s]
VU = LU / TU;                          % Velocity unit [km/s]
radius_Earth = 6371;                   % Earth radius [km]
radius_Moon  = 1737.1;                 % Moon radius  [km]
T            = 2357418;                % One lunar period [s]
beta_moon = atan(radius_Moon / LU);    % Moon occulation angle
DSN.omega = 7.2921159e-5;              % Earth rotational rate [rad/s]

%% ------------ DSN Station coortinates -------------
% ── Latitude(°), Longitude(°E), Altitude(km) ──
DSN.lat  =  [ 35.2472,  40.4314, -35.4020 ];     % Goldstone, Madrid, Canberra
DSN.lon  = [-116.7933,  -4.2486, 148.9819 ];     % East longitude is positive
DSN.alt  =  [    1.0,     0.9,      0.7   ];     

%% ------------ Initialisation -------------
Period     = 1.3987087584241628E+1;    % days
Period_TU  = Period / 4.43;            % [TU]

hrs = 300;                             % simulation span in hrs
span = (hrs/24) / 4.43;                % simulation span
dt_p = 0.000157 / 5;                   % propagation step (0.2 min)

arc_dro  = 0.0188 * 100 * 100;         % 2 h
arc_elfo = 0.0188 * 100 * 100;         % 2 h

% Monte Carlo setting
kRuns = 100;                           % Times of run
rng(0);                                % Fixed seed

% UKF parameters
alpha = 0.5;
beta = 2;

% For SST
switch nRelay
    case 1
        alpha_sst_p = 1;   beta_sst_p = 2;     kappa_sst_p = 0;
        alpha_sst_m = 1;     beta_sst_m = -10;   kappa_sst_m = 0.2;
        SST_noise_factor = 1;
        dt_sst_try = 280;
    case 2
        alpha_sst_p = 0.5;   beta_sst_p = 2;    kappa_sst_p = 5;
        alpha_sst_m = 0.5;   beta_sst_m = -4;   kappa_sst_m = 6;
        SST_noise_factor = 1;
        dt_sst_try = 330;
    case 3
        alpha_sst_p = 0.5;   beta_sst_p = 2;    kappa_sst_p = 5;
        alpha_sst_m = 0.5;   beta_sst_m = -4;   kappa_sst_m = 6;
        SST_noise_factor = 1;
        dt_sst_try = 360;
end

%% ------------ Base measurement interval -------------
dt_dsn_ref  = 0.000157;     % Ground (1min)
dt_sst_ref = 0.000157;      % SST (1min)

%% ------------ Actual measurement interval -------------
for i = 1:length(dt_sst_try)

dt_sst_min = dt_sst_try(i);

dt_dro   = 0.000157*2;      % Ground (2min)
dt_elfo  = dt_dro;          % Ground (2min)
dt_sst   = 0.000157*dt_sst_min;     % SST (0.5min)

%% ------------ Measurement Noise（Adaptive） ----------------
% Scaling factor
scale_dsn = dt_dsn_ref / dt_dro;        % DSN: Δt_ref / Δt
scale_sst = dt_sst_ref / dt_sst;        % ISL: Δt_ref / Δt

% Ground
R_dsn_base = diag([0.001/LU, 50E-6, 50E-6, 0.1e-6/VU].^2);  % Base covariances (σ²)

% Actual covariances(Δt scaling)
R = R_dsn_base * R_GB_factor;
R(1,1) = R(1,1) * scale_dsn;        % range variance ∝ 1/Δt
R(4,4) = R(4,4) * scale_dsn^3;      % Doppler variance ∝ 1/Δt^3

% SST standard deviation
sigma_range     = 0.001/LU;
sigma_doppler   = 0.1e-6/VU;

% Scaling
range_var     = (sigma_range^2)   * scale_sst;
doppler_var   = (sigma_doppler^2) * scale_sst^3;

R_ISL1 = diag([range_var, doppler_var]);
R_ISL2 = diag([range_var, range_var, doppler_var, doppler_var]);
R_ISL3 = diag([range_var, range_var, range_var, doppler_var, doppler_var, doppler_var]);

R_ISL1_filt = diag([range_var, doppler_var]) * R_sst_factor;
R_ISL2_filt = diag([range_var, range_var, doppler_var, doppler_var]) * R_sst_factor;
R_ISL3_filt = diag([range_var, range_var, range_var, doppler_var, doppler_var, doppler_var]) * R_sst_factor;

%% ------------ Process noise (Adaptive) ----------------
%% Constants
scaleTU2_LU = (TU^2) / LU;   % 2.6608e‑6

%% ---------- DRO ----------
sigma_a_truth_dro   = 1e-12;          % km/s^2
sigma_a_filter_dro  = sigma_a_truth_dro * 1;

aT_nd = sigma_a_truth_dro * scaleTU2_LU;
aF_nd = sigma_a_filter_dro * scaleTU2_LU;

Q_dro_truth = diag((aT_nd).^2 * ones(1,3));
Q_dro_filt  = diag((aF_nd).^2 * ones(1,3));

S_i_dro_truth = [Q_dro_truth*dt_p^3/3, Q_dro_truth*dt_p^2/2;
                 Q_dro_truth*dt_p^2/2, Q_dro_truth*dt_p];
S_i_dro_filt  = [Q_dro_filt*dt_p^3/3,  Q_dro_filt*dt_p^2/2;
                 Q_dro_filt*dt_p^2/2,  Q_dro_filt*dt_p];

%% ---------- ELFO ----------
sigma_a_truth_elfo   = 3e-12;
sigma_a_filter_elfo  = sigma_a_truth_elfo * 1;

aT_nd = sigma_a_truth_elfo * scaleTU2_LU;
aF_nd = sigma_a_filter_elfo * scaleTU2_LU;

Q_elfo_truth = diag((aT_nd).^2 * ones(1,3));
Q_elfo_filt  = diag((aF_nd).^2 * ones(1,3));

S_i_elfo_truth = [Q_elfo_truth*dt_p^3/3, Q_elfo_truth*dt_p^2/2;
                  Q_elfo_truth*dt_p^2/2, Q_elfo_truth*dt_p];
S_i_elfo_filt  = [Q_elfo_filt*dt_p^3/3,  Q_elfo_filt*dt_p^2/2;
                  Q_elfo_filt*dt_p^2/2,  Q_elfo_filt*dt_p];

%% ---------- SST ----------
sigma_a_filter_SST  = sigma_a_truth_elfo * SST_noise_factor;

aF_nd = sigma_a_filter_SST * scaleTU2_LU;

Q_sst_filt  = diag((aF_nd).^2 * ones(1,3));

S_i_sst_filt  = [Q_sst_filt*dt_p^3/3,  Q_sst_filt*dt_p^2/2;
                  Q_sst_filt*dt_p^2/2,  Q_sst_filt*dt_p];

disp(['alpha=', num2str(alpha), ', beta=', num2str(beta), ...
    ', alpha_sst=', num2str(alpha_sst_m), ', beta_sst=', num2str(beta_sst_m), ', kappa_sst=', num2str(kappa_sst_m), ...
    ', dt_sst=', num2str(dt_sst_min), ', SST noise factor=', num2str(SST_noise_factor)]);

%% ------------ Monte Carlo Initialisations ----------------
nLength = ceil(span / dt_p) + 1;

% Error related computations for DRO 1
dro1_errors = zeros(nLength, 6, kRuns);
dro1_rErrorsNorm = zeros(nLength, 1, kRuns); % In km
dro1_vErrorsNorm = zeros(nLength, 1, kRuns); % In m/s
dro1_NEES_M = zeros(nLength, kRuns);
dro1_NIS_M = zeros(nLength, kRuns);

% Error related computations for DRO 2
dro2_errors = zeros(nLength, 6, kRuns);
dro2_rErrorsNorm = zeros(nLength, 1, kRuns); % In km
dro2_vErrorsNorm = zeros(nLength, 1, kRuns); % In m/s
dro2_NEES_M = zeros(nLength, kRuns);
dro2_NIS_M = zeros(nLength, kRuns);

% Error related computations for DRO 3
dro3_errors = zeros(nLength, 6, kRuns);
dro3_rErrorsNorm = zeros(nLength, 1, kRuns); % In km
dro3_vErrorsNorm = zeros(nLength, 1, kRuns); % In m/s
dro3_NEES_M = zeros(nLength, kRuns);
dro3_NIS_M = zeros(nLength, kRuns);

% Error related computations for ELFO
elfo_errors = zeros(nLength, 6, kRuns);
elfo_rErrorsNorm = zeros(nLength, 1, kRuns); % In km
elfo_vErrorsNorm = zeros(nLength, 1, kRuns); % In m/s
elfo_NEES_M = zeros(nLength, kRuns);
elfo_NIS_M = zeros(nLength, kRuns);

% Error related computations for ELFO using ISL
sst_errors = zeros(nLength, 6, kRuns);
sst_rErrorsNorm = zeros(nLength, 1, kRuns); % In km
sst_vErrorsNorm = zeros(nLength, 1, kRuns); % In m/s
sst_NEES_M = zeros(nLength, kRuns);
sst_NIS_M = zeros(nLength, kRuns);

dro1_ellipsoid_r_k = zeros(nLength, kRuns);
dro1_ellipsoid_v_k = zeros(nLength, kRuns);
dro2_ellipsoid_r_k = zeros(nLength, kRuns);
dro2_ellipsoid_v_k = zeros(nLength, kRuns);
dro3_ellipsoid_r_k = zeros(nLength, kRuns);
dro3_ellipsoid_v_k = zeros(nLength, kRuns);
elfo_ellipsoid_r_k = zeros(nLength, kRuns);
elfo_ellipsoid_v_k = zeros(nLength, kRuns);
sst_ellipsoid_r_k = zeros(nLength, kRuns);
sst_ellipsoid_v_k = zeros(nLength, kRuns);

elfo_visArrays = zeros(nLength, 1, kRuns);
dro1_visArrays = zeros(nLength, 1, kRuns);
sst_visArrays = zeros(nLength, 1, kRuns);

elfo_arcArrays = zeros(nLength, 1, kRuns);
dro1_arcArrays = zeros(nLength, 1, kRuns);
sst_arcArrays = zeros(nLength, 1, kRuns);

%% ------------ Main Function ----------------
parfor k = 1:kRuns

    rng(k);

    % Logical mask indicating which DROs participate as relays
    isRelay = false(1,3);
    isRelay(1:nRelay) = true;

    time = 0;
    last_measure_time1 = 0;  last_measure_time2 = 0;  last_measure_time3 = 0;
    last_elfo_measure_time = 0;  last_ISL_measure_time = 0;
    condition_met_count = 0;  n = 1;
    cumulative_duration_elfo = 0;
    cumulative_duration_dro1 = 0;  cumulative_duration_dro2 = 0;  cumulative_duration_dro3 = 0;

    %% STATE‑VECTOR INITIALISATIONS
    % DRO‑1
    xTrue_dro1 = [0.8082345151982595	2.230630555067976e-23	4.300655269891714e-25	...
        -8.433177870650808e-13	0.5164471457797999	-1.393414024105238e-25]';

    % DRO‑2
    xTrue_dro2 = [1.0829328987799418	0.20511422300703588	-2.1728733456390955e-25	...
        0.31976943762396287	-0.2515704044237798	-6.818359799318023e-25]';

    if nRelay == 2 % 180° phase distribution for the 2nd if 2 DROs
        xTrue_dro2 = [1.176169613224601	-0.0006842568544340667	-4.789908432217674e-25	...
            -0.0011832404429543827	-0.4949076145781299	-2.841763929488676e-25]';
    end

    % DRO‑3
    xTrue_dro3 = [1.084810773715486	-0.20361889753888288	-4.817724235727507e-25	...
        -0.3179642049417935	-0.25625067738359253	2.4334825883068064e-25]';

    % ELFO
    r0_elfo = [0.985300655778438 0.00348235528237756 0.00557943029159241];
    v0_elfo = [-1.54413275759651 -0.614032806699692 8.65832686101901e-17];
    xTrue_elfo = [r0_elfo, v0_elfo]';

    % NRHO (not acting as relay)
    r0_nrho = [1.021176128690498 0.00901228052473247 -0.18149355071754605];
    v0_nrho = [-0.0008949603180085014 -0.10140741960410796 -1.67619460939642e-12];
    
    %%  VISIBILITY / RELAY‑PATH CHECKS (modified to honour nRelay)
    % Calculate visibility of DRO from Earth
    rTrue_dro1 = norm([xTrue_dro1(1) + mu, xTrue_dro1(2), xTrue_dro1(3)]);
    yznorm_dro1 = (xTrue_dro1(2)^2 + xTrue_dro1(3)^2)^0.5;
    theta_dro1 = asin(yznorm_dro1 / rTrue_dro1);
    
    % Calculate visibility of ELFO from Earth
    rTrue_elfo = norm([xTrue_elfo(1) + mu, xTrue_elfo(2), xTrue_elfo(3)]);
    yznorm_elfo = (xTrue_elfo(2)^2 + xTrue_elfo(3)^2)^0.5;
    theta_elfo = asin(yznorm_elfo / rTrue_elfo);
    
    % DRO‑1 visibility
    r_dro_moon1 = [1 - mu; 0; 0] - xTrue_dro1(1:3);
    r_dro_elfo1 = xTrue_elfo(1:3) - xTrue_dro1(1:3);
    theta_proj1 = asin(norm(cross(r_dro_moon1, r_dro_elfo1)) /(norm(r_dro_moon1)*norm(r_dro_elfo1)));
    beta1       = atan(radius_Moon/LU/ norm(r_dro_moon1));
    is_visible1 = theta_proj1 > beta1 || norm(r_dro_elfo1)*cos(theta_proj1) < norm(r_dro_moon1);
    
    % DRO‑2 visibility
    r_dro_moon2 = [1 - mu; 0; 0] - xTrue_dro2(1:3);
    r_dro_elfo2 = xTrue_elfo(1:3) - xTrue_dro2(1:3);
    theta_proj2 = asin(norm(cross(r_dro_moon2, r_dro_elfo2)) /(norm(r_dro_moon2)*norm(r_dro_elfo2)));
    beta2       = atan(radius_Moon/LU/ norm(r_dro_moon2));
    is_visible2 = theta_proj2 > beta2 || norm(r_dro_elfo2)*cos(theta_proj2) < norm(r_dro_moon2);
    
    % DRO‑3 visibility
    r_dro_moon3 = [1 - mu; 0; 0] - xTrue_dro3(1:3);
    r_dro_elfo3 = xTrue_elfo(1:3) - xTrue_dro3(1:3);
    theta_proj3 = asin(norm(cross(r_dro_moon3, r_dro_elfo3)) /(norm(r_dro_moon3)*norm(r_dro_elfo3)));
    beta3       = atan(radius_Moon/LU/ norm(r_dro_moon3));
    is_visible3 = theta_proj3 > beta3 || norm(r_dro_elfo3)*cos(theta_proj3) < norm(r_dro_moon3);
    
    % >>> NEW: enforce user‑selected relay configuration <<<
    is_visible1 = isRelay(1) && is_visible1;
    is_visible2 = isRelay(2) && is_visible2;
    is_visible3 = isRelay(3) && is_visible3;
    
    %% ------------ Pre-allocating Spaces for Results----------------
    xDR_dro1 = xTrue_dro1;
    dro1_xDR = zeros(nLength, 6);
    dro1_xDR(1, :) = xDR_dro1';
    dro1_xTrue = zeros(nLength, 6);
    dro1_xTrue(1, :) = xTrue_dro1';
    dro1_xEst = zeros(nLength, 6);
    
    xDR_dro2 = xTrue_dro2;
    dro2_xDR = zeros(nLength, 6);
    dro2_xDR(1, :) = xDR_dro2';
    dro2_xTrue = zeros(nLength, 6);
    dro2_xTrue(1, :) = xTrue_dro2';
    dro2_xEst = zeros(nLength, 6);
    
    xDR_dro3 = xTrue_dro3;
    dro3_xDR = zeros(nLength, 6);
    dro3_xDR(1, :) = xDR_dro3';
    dro3_xTrue = zeros(nLength, 6);
    dro3_xTrue(1, :) = xTrue_dro3';
    dro3_xEst = zeros(nLength, 6);
    
    xDR_elfo = xTrue_elfo;
    elfo_xDR = zeros(nLength, 6);
    elfo_xDR(1, :) = xDR_elfo';
    elfo_xTrue = zeros(nLength, 6);
    elfo_xTrue(1, :) = xTrue_elfo';
    elfo_xEst = zeros(nLength, 6);
    
    sst_xEst = zeros(nLength, 6);

    dro1_visSlice = zeros(nLength, 1);
    elfo_visSlice = zeros(nLength, 1);
    sst_visSlice = zeros(nLength, 1);
    
    dro1_arcSlice = zeros(nLength, 1);
    elfo_arcSlice = zeros(nLength, 1);
    sst_arcSlice = ones(nLength, 1);

    dro1_ellipsoid_r = zeros(nLength, 1);
    dro1_ellipsoid_v = zeros(nLength, 1);
    dro2_ellipsoid_r = zeros(nLength, 1);
    dro2_ellipsoid_v = zeros(nLength, 1);
    dro3_ellipsoid_r = zeros(nLength, 1);
    dro3_ellipsoid_v = zeros(nLength, 1);
    elfo_ellipsoid_r = zeros(nLength, 1);
    elfo_ellipsoid_v = zeros(nLength, 1);
    sst_ellipsoid_r = zeros(nLength, 1);
    sst_ellipsoid_v = zeros(nLength, 1);

    dro1_NEES = zeros(nLength, 1);
    dro2_NEES = zeros(nLength, 1);
    dro3_NEES = zeros(nLength, 1);
    elfo_NEES = zeros(nLength, 1);
    sst_NEES = zeros(nLength, 1);

    dro1_NIS = zeros(nLength, 1);
    dro2_NIS = zeros(nLength, 1);
    dro3_NIS = zeros(nLength, 1);
    elfo_NIS = zeros(nLength, 1);
    sst_NIS = zeros(nLength, 1);
    
    if theta_dro1 > beta_moon || xTrue_dro1(1) < 1-mu
        dro1_visSlice(n) = 1;
        if cumulative_duration_dro1 <= arc_dro
            dro1_arcSlice(n) = 1;
        end
    end
    
    if theta_elfo > beta_moon || xTrue_elfo(1) < 1-mu
        elfo_visSlice(n) = 1;
        if cumulative_duration_elfo <= arc_elfo
            elfo_arcSlice(n) = 1;
        end
    end
    
    % Visibility array
    ISL_visibility = [is_visible1, is_visible2, is_visible3];
    sst_visSlice(n) = sum(ISL_visibility);
    
    %% ------------ Initial Measurements ----------------
    meas_noise = sqrt(diag(R)) .* randn(4,1);

    idxSta = choose_station(xTrue_dro1, time, DSN, 10);
    rGS = dsn_pos_CR3BP(idxSta, time, DSN);
    vGS = cross([0 0 DSN.omega]', rGS*LU)/LU;
    y_dro1 = observation_model_H(xTrue_dro1, rGS, vGS, meas_noise);
    y_coordinates_dro1 = calculate_coordinates(y_dro1, rGS, vGS)';

    idxSta = choose_station(xTrue_dro2, time, DSN, 10);
    rGS = dsn_pos_CR3BP(idxSta, time, DSN);
    vGS = cross([0 0 DSN.omega]', rGS*LU)/LU;
    y_dro2 = observation_model_H(xTrue_dro2, rGS, vGS, meas_noise);
    y_coordinates_dro2 = calculate_coordinates(y_dro2, rGS, vGS)';

    idxSta = choose_station(xTrue_dro3, time, DSN, 10);
    rGS = dsn_pos_CR3BP(idxSta, time, DSN);
    vGS = cross([0 0 DSN.omega]', rGS*LU)/LU;
    y_dro3 = observation_model_H(xTrue_dro3, rGS, vGS, meas_noise);
    y_coordinates_dro3 = calculate_coordinates(y_dro3, rGS, vGS)';

    idxSta = choose_station(xTrue_elfo, time, DSN, 10);
    rGS = dsn_pos_CR3BP(idxSta, time, DSN);
    vGS = cross([0 0 DSN.omega]', rGS*LU)/LU;
    y_elfo = observation_model_H(xTrue_elfo, rGS, vGS, meas_noise);
    y_coordinates_elfo = calculate_coordinates(y_elfo, rGS, vGS)';
    
    dro1_y = y_coordinates_dro1;
    dro2_y = y_coordinates_dro2;
    dro3_y = y_coordinates_dro3;
    elfo_y = y_coordinates_elfo;
    
    %% ------------ Initial UKF Estimation ----------------
    P0 = diag(([1/LU 1/LU 1/LU 1e-3/VU 1e-3/VU 1e-3/VU]).^2);
    initial_position_error = 10/LU;   initial_velocity_error = 1e-3/VU;
    
    xEst_dro1 = xTrue_dro1 + [initial_position_error*randn(3,1); initial_velocity_error*randn(3,1)];
    PEst_dro1 = P0;
    NEES_dro1 = (xTrue_dro1-xEst_dro1)' * PEst_dro1^(-1) * (xTrue_dro1-xEst_dro1);
    NIS_dro1 = 0;

    xEst_dro2 = xTrue_dro2 + [initial_position_error*randn(3,1); initial_velocity_error*randn(3,1)];
    PEst_dro2 = P0;  
    NEES_dro2 = (xTrue_dro2-xEst_dro2)' * PEst_dro2^(-1) * (xTrue_dro2-xEst_dro2);
    NIS_dro2 = 0;

    xEst_dro3 = xTrue_dro3 + [initial_position_error*randn(3,1); initial_velocity_error*randn(3,1)];
    PEst_dro3 = P0;  
    NEES_dro3 = (xTrue_dro3-xEst_dro3)' * PEst_dro3^(-1) * (xTrue_dro3-xEst_dro3);
    NIS_dro3 = 0;

    xEst_elfo = xTrue_elfo + [initial_position_error*randn(3,1); initial_velocity_error*randn(3,1)];
    PEst_elfo = P0;  
    NEES_elfo = (xTrue_elfo-xEst_elfo)' * PEst_elfo^(-1) * (xTrue_elfo-xEst_elfo);
    NIS_elfo = 0;

    xEst_sst = xEst_elfo;
    PEst_sst = P0;
    NEES_sst = (xTrue_elfo-xEst_sst)' * PEst_sst^(-1) * (xTrue_elfo-xEst_sst);
    NIS_sst = 0;

    dro1_xEst(1, :) = xEst_dro1';
    dro2_xEst(1, :) = xEst_dro2';
    dro3_xEst(1, :) = xEst_dro3';
    elfo_xEst(1, :) = xEst_elfo';
    sst_xEst(1, :) = xEst_sst';

    dro1_NEES(1) = NEES_dro1;
    dro2_NEES(1) = NEES_dro2;
    dro3_NEES(1) = NEES_dro3;
    elfo_NEES(1) = NEES_elfo;
    sst_NEES(1) = NEES_sst;

    dro1_NIS(1) = NIS_dro1;
    dro2_NIS(1) = NIS_dro2;
    dro3_NIS(1) = NIS_dro3;
    elfo_NIS(1) = NIS_elfo;
    sst_NIS(1) = NIS_sst;

    dro1_ellipsoid_r(1) = 4*pi/3*sqrt(det(PEst_dro1(1:3,1:3)*(LU*1e3)^2)); % m^3
    dro1_ellipsoid_v(1) = 4*pi/3*sqrt(det(PEst_dro1(4:6,4:6)*(VU*1e6)^2)); % mm^3/s^3
    dro2_ellipsoid_r(1) = 4*pi/3*sqrt(det(PEst_dro2(1:3,1:3)*(LU*1e3)^2));
    dro2_ellipsoid_v(1) = 4*pi/3*sqrt(det(PEst_dro2(4:6,4:6)*(VU*1e6)^2));
    dro3_ellipsoid_r(1) = 4*pi/3*sqrt(det(PEst_dro3(1:3,1:3)*(LU*1e3)^2));
    dro3_ellipsoid_v(1) = 4*pi/3*sqrt(det(PEst_dro3(4:6,4:6)*(VU*1e6)^2));
    elfo_ellipsoid_r(1) = 4*pi/3*sqrt(det(PEst_elfo(1:3,1:3)*(LU*1e3)^2));
    elfo_ellipsoid_v(1) = 4*pi/3*sqrt(det(PEst_elfo(4:6,4:6)*(VU*1e6)^2));
    sst_ellipsoid_r(1) = 4*pi/3*sqrt(det(PEst_sst(1:3,1:3)*(LU*1e3)^2));
    sst_ellipsoid_v(1) = 4*pi/3*sqrt(det(PEst_sst(4:6,4:6)*(VU*1e6)^2));
    
    %% ------------ Main Loop ----------------
    while span >= time
        time = time + dt_p;  % Increment time by propagation step
        % disp(['time=' num2str(time)]);
        n = n + 1;
    
        % Simulate true state propagation
        process_noise_dro = chol(S_i_dro_truth)' * randn(6, 1); % Generate 6x1 noise vector
        xTrue_dro1 = motion_model_F(xTrue_dro1, process_noise_dro, dt_p);
        process_noise_dro = chol(S_i_dro_truth)' * randn(6, 1); % Generate 6x1 noise vector
        xTrue_dro2 = motion_model_F(xTrue_dro2, process_noise_dro, dt_p);
        process_noise_dro = chol(S_i_dro_truth)' * randn(6, 1); % Generate 6x1 noise vector
        xTrue_dro3 = motion_model_F(xTrue_dro3, process_noise_dro, dt_p);
        process_noise_elfo = chol(S_i_elfo_truth)' * randn(6, 1); % Generate 6x1 noise vector
        xTrue_elfo = motion_model_F(xTrue_elfo, process_noise_elfo, dt_p);
    
        % Update dead reckoning state
        xDR_dro1 = motion_model_F(xDR_dro1, [], dt_p);  % No process noise for dead reckoning
        xDR_dro2 = motion_model_F(xDR_dro2, [], dt_p);  % No process noise for dead reckoning
        xDR_dro3 = motion_model_F(xDR_dro3, [], dt_p);  % No process noise for dead reckoning
        xDR_elfo = motion_model_F(xDR_elfo, [], dt_p);  % No process noise for dead reckoning
    
        % Calculate visibility of DRO from Earth
        rTrue_dro1 = norm([xTrue_dro1(1) + mu, xTrue_dro1(2), xTrue_dro1(3)]);
        yznorm_dro1 = (xTrue_dro1(2)^2 + xTrue_dro1(3)^2)^0.5;
        theta_dro1 = asin(yznorm_dro1 / rTrue_dro1);
    
        rTrue_dro2 = norm([xTrue_dro2(1) + mu, xTrue_dro2(2), xTrue_dro2(3)]);
        yznorm_dro2 = (xTrue_dro2(2)^2 + xTrue_dro2(3)^2)^0.5;
        theta_dro2 = asin(yznorm_dro2 / rTrue_dro2);
    
        rTrue_dro3 = norm([xTrue_dro3(1) + mu, xTrue_dro3(2), xTrue_dro3(3)]);
        yznorm_dro3 = (xTrue_dro3(2)^2 + xTrue_dro3(3)^2)^0.5;
        theta_dro3 = asin(yznorm_dro3 / rTrue_dro3);
    
        % Calculate visibility of ELFO from Earth
        rTrue_elfo = norm([xTrue_elfo(1) + mu, xTrue_elfo(2), xTrue_elfo(3)]);
        yznorm_elfo = (xTrue_elfo(2)^2 + xTrue_elfo(3)^2)^0.5;
        theta_elfo = asin(yznorm_elfo / rTrue_elfo);
        
        % DRO 1 Visibility to ELFO
        r_dro_moon1 = [1 - mu; 0; 0] - xTrue_dro1(1:3);
        r_dro_elfo1 = xTrue_elfo(1:3) - xTrue_dro1(1:3);
        theta_proj1 = asin(norm(cross(r_dro_moon1, r_dro_elfo1)) / (norm(r_dro_moon1) * norm(r_dro_elfo1)));
        beta1 = atan(radius_Moon / LU / norm(r_dro_moon1));
        is_visible1 = theta_proj1 > beta1 || norm(r_dro_elfo1) * cos(theta_proj1) < norm(r_dro_moon1);
        
        % DRO 2 Visibility to ELFO
        r_dro_moon2 = [1 - mu; 0; 0] - xTrue_dro2(1:3);
        r_dro_elfo2 = xTrue_elfo(1:3) - xTrue_dro2(1:3);
        theta_proj2 = asin(norm(cross(r_dro_moon2, r_dro_elfo2)) / (norm(r_dro_moon2) * norm(r_dro_elfo2)));
        beta2 = atan(radius_Moon / LU / norm(r_dro_moon2));
        is_visible2 = theta_proj2 > beta2 || norm(r_dro_elfo2) * cos(theta_proj2) < norm(r_dro_moon2);
        
        % DRO 3 Visibility to ELFO
        r_dro_moon3 = [1 - mu; 0; 0] - xTrue_dro3(1:3);
        r_dro_elfo3 = xTrue_elfo(1:3) - xTrue_dro3(1:3);
        theta_proj3 = asin(norm(cross(r_dro_moon3, r_dro_elfo3)) / (norm(r_dro_moon3) * norm(r_dro_elfo3)));
        beta3 = atan(radius_Moon / LU / norm(r_dro_moon3));
        is_visible3 = theta_proj3 > beta3 || norm(r_dro_elfo3) * cos(theta_proj3) < norm(r_dro_moon3);
    
        is_visible1 = isRelay(1) && is_visible1;
        is_visible2 = isRelay(2) && is_visible2;
        is_visible3 = isRelay(3) && is_visible3;
        
        % SST visibility array
        ISL_visibility = [is_visible1, is_visible2, is_visible3];
        sst_visSlice(n) = sum(ISL_visibility);
    
        % Observation arc
        cumulative_duration_elfo = cumulative_duration_elfo + dt_p;
        cumulative_duration_dro1 = cumulative_duration_dro1 + dt_p;
        cumulative_duration_dro2 = cumulative_duration_dro2 + dt_p;
        cumulative_duration_dro3 = cumulative_duration_dro3 + dt_p;
    
        %% ELFO
        idxSta = choose_station(xTrue_elfo, time, DSN, 10);
        if idxSta ~= 0
            elfo_visSlice(n) = idxSta;
            if cumulative_duration_elfo <= arc_elfo
                elfo_arcSlice(n) = 1;
            end
        end
    
        if time - last_elfo_measure_time >= dt_elfo && idxSta ~= 0 && cumulative_duration_elfo <= arc_elfo
            rGS = dsn_pos_CR3BP(idxSta, time, DSN);
            vGS = cross([0 0 DSN.omega]', rGS*LU)/LU; 

            % Simulate observation
            meas_noise = sqrt(diag(R)) .* randn(4, 1);
            y_elfo = observation_model_H(xTrue_elfo, rGS, vGS, meas_noise);
    
            % Perform UKF estimation with measurement update
            [xEst_elfo, PEst_elfo, NIS_elfo] = UKF(xEst_elfo, PEst_elfo, y_elfo, S_i_elfo_filt, R, ...
                dt_p, rGS, vGS, alpha, beta);
            last_elfo_measure_time = time;
        else
            % Perform UKF state propagation without measurement update
            [xEst_elfo, PEst_elfo, ~] = UKF(xEst_elfo, PEst_elfo, [], S_i_elfo_filt, [], ...
                dt_p, [], [], alpha, beta);
        end

        NEES_elfo = (xTrue_elfo-xEst_elfo)' * PEst_elfo^(-1) * (xTrue_elfo-xEst_elfo);
    
        if cumulative_duration_elfo >= 0.2256
            cumulative_duration_elfo = 0;
        end
    
        %% DRO 1
        idxSta = choose_station(xTrue_dro1, time, DSN, 10);
        if idxSta ~= 0
            dro1_visSlice(n) = idxSta;
            if cumulative_duration_dro1 <= arc_dro
                dro1_arcSlice(n) = 1;
            end
        end
    
        if time - last_measure_time1 >= dt_dro && idxSta ~= 0 && cumulative_duration_dro1 <= arc_dro
            rGS = dsn_pos_CR3BP(idxSta, time, DSN);
            vGS = cross([0 0 DSN.omega]', rGS*LU)/LU; 

            % Simulate observation for DRO 1
            meas_noise = sqrt(diag(R)) .* randn(4, 1);
            y_dro1 = observation_model_H(xTrue_dro1, rGS, vGS, meas_noise);
    
            % Perform UKF estimation with measurement update
            [xEst_dro1, PEst_dro1, NIS_dro1] = UKF(xEst_dro1, PEst_dro1, y_dro1, S_i_dro_filt, R, ...
                dt_p, rGS, vGS, alpha, beta);
            last_measure_time1 = time;
        else
            % Perform UKF state propagation without measurement update
            [xEst_dro1, PEst_dro1, ~] = UKF(xEst_dro1, PEst_dro1, [], S_i_dro_filt, [], ...
                dt_p, [], [], alpha, beta);
        end

        NEES_dro1 = (xTrue_dro1-xEst_dro1)' * PEst_dro1^(-1) * (xTrue_dro1-xEst_dro1);
    
        if cumulative_duration_dro1 >= 0.2256 % Length of 1 day in normalised time unit
            cumulative_duration_dro1 = 0;
        end
    
        %% DRO 2
        if nRelay == 2 || nRelay == 3
            idxSta = choose_station(xTrue_dro2, time, DSN, 10);
            if time - last_measure_time2 >= dt_dro && idxSta ~= 0 && cumulative_duration_dro2 <= arc_dro
                rGS = dsn_pos_CR3BP(idxSta, time, DSN);
                vGS = cross([0 0 DSN.omega]', rGS*LU)/LU; 
    
                % Simulate observation for DRO 2
                meas_noise = sqrt(diag(R)) .* randn(4, 1);
                y_dro2 = observation_model_H(xTrue_dro2, rGS, vGS, meas_noise);
        
                % Perform UKF estimation with measurement update
                [xEst_dro2, PEst_dro2, NIS_dro2] = UKF(xEst_dro2, PEst_dro2, y_dro2, S_i_dro_filt, R, ...
                    dt_p, rGS, vGS, alpha, beta);
                last_measure_time2 = time;
            else
                % Perform UKF state propagation without measurement update
                [xEst_dro2, PEst_dro2, ~] = UKF(xEst_dro2, PEst_dro2, [], S_i_dro_filt, [], ...
                    dt_p, [], [], alpha, beta);
            end
        
            NEES_dro2 = (xTrue_dro2-xEst_dro2)' * PEst_dro2^(-1) * (xTrue_dro2-xEst_dro2);
    
            if cumulative_duration_dro2 >= 0.2256
                cumulative_duration_dro2 = 0;
            end
        end

        %% DRO 3
        if nRelay == 3
            idxSta = choose_station(xTrue_dro3, time, DSN, 10);
            if time - last_measure_time3 >= dt_dro && idxSta ~= 0 && cumulative_duration_dro3 <= arc_dro
                rGS = dsn_pos_CR3BP(idxSta, time, DSN);
                vGS = cross([0 0 DSN.omega]', rGS*LU)/LU; 
    
                % Simulate observation for DRO 3
                meas_noise = sqrt(diag(R)) .* randn(4, 1);
                y_dro3 = observation_model_H(xTrue_dro3, rGS, vGS, meas_noise);
        
                % Perform UKF estimation with measurement update
                [xEst_dro3, PEst_dro3, NIS_dro3] = UKF(xEst_dro3, PEst_dro3, y_dro3, S_i_dro_filt, R, ...
                    dt_p, rGS, vGS, alpha, beta);
                last_measure_time3 = time;
            else
                % Perform UKF state propagation without measurement update
                [xEst_dro3, PEst_dro3, ~] = UKF(xEst_dro3, PEst_dro3, [], S_i_dro_filt, [], ...
                    dt_p, [], [], alpha, beta);
            end
    
            NEES_dro3 = (xTrue_dro3-xEst_dro3)' * PEst_dro3^(-1) * (xTrue_dro3-xEst_dro3);
        
            if cumulative_duration_dro3 >= 0.2256
                cumulative_duration_dro3 = 0;
            end
        end
    
        %% ELFO of SST
        if time - last_ISL_measure_time >= dt_sst && sum(ISL_visibility) == 3
    
            % Simulate ISL observation
            ISL_noise = sqrt(diag(R_ISL3)) .* randn(6, 1);
            delta = ISL_model(xTrue_elfo, xTrue_dro1, xTrue_dro2, xTrue_dro3, ISL_noise);
    
            % Perform UKF estimation with measurement update
            [xEst_sst, PEst_sst, NIS_sst] = UKF_ISL(xEst_sst, PEst_sst, delta, S_i_sst_filt, R_ISL3_filt, dt_p, ...
                xEst_dro1, PEst_dro1, xEst_dro2, PEst_dro2, xEst_dro3, PEst_dro3, alpha_sst_p, beta_sst_p, kappa_sst_p, alpha_sst_m, beta_sst_m, kappa_sst_m);
            last_ISL_measure_time = time;
    
        elseif time - last_ISL_measure_time >= dt_sst && sum(ISL_visibility) == 2
    
            invisible_idx = find(ISL_visibility == 0);

            if invisible_idx == 1
                xEst_dro_1 = xEst_dro2;
                PEst_dro_1 = PEst_dro2;
                xTrue_dro_1 = xTrue_dro2;
                xEst_dro_2 = xEst_dro3;
                PEst_dro_2 = PEst_dro3;
                xTrue_dro_2 = xTrue_dro3;
    
            elseif invisible_idx == 2
                xEst_dro_1 = xEst_dro1;
                PEst_dro_1 = PEst_dro1;
                xTrue_dro_1 = xTrue_dro1;
                xEst_dro_2 = xEst_dro3;
                PEst_dro_2 = PEst_dro3;
                xTrue_dro_2 = xTrue_dro3;
    
            elseif invisible_idx == 3
                xEst_dro_1 = xEst_dro1;
                PEst_dro_1 = PEst_dro1;
                xTrue_dro_1 = xTrue_dro1;
                xEst_dro_2 = xEst_dro2;
                PEst_dro_2 = PEst_dro2;
                xTrue_dro_2 = xTrue_dro2;
            end
    
            % Simulate ISL observation
            ISL_noise = sqrt(diag(R_ISL2)) .* randn(4, 1);
            delta = ISL_model(xTrue_elfo, xTrue_dro_1, xTrue_dro_2, [], ISL_noise);
    
            % Perform UKF estimation with measurement update
            [xEst_sst, PEst_sst, NIS_sst] = UKF_ISL(xEst_sst, PEst_sst, delta, S_i_sst_filt, R_ISL2_filt, dt_p, ...
                xEst_dro_1, PEst_dro_1, xEst_dro_2, PEst_dro_2, [], [], alpha_sst_p, beta_sst_p, kappa_sst_p, alpha_sst_m, beta_sst_m, kappa_sst_m);
            last_ISL_measure_time = time;
    
        elseif time - last_ISL_measure_time >= dt_sst && sum(ISL_visibility) == 1
    
            visible_idx = find(ISL_visibility, 1);
    
            if visible_idx == 1
                xEst_dro = xEst_dro1;
                PEst_dro = PEst_dro1;
                xTrue_dro = xTrue_dro1;
            elseif visible_idx == 2
                xEst_dro = xEst_dro2;
                PEst_dro = PEst_dro2;
                xTrue_dro = xTrue_dro2;
            else
                xEst_dro = xEst_dro3;
                PEst_dro = PEst_dro3;
                xTrue_dro = xTrue_dro3;
            end
    
            % Simulate ISL observation
            ISL_noise = sqrt(diag(R_ISL1)) .* randn(2, 1);
            delta = ISL_model(xTrue_elfo, xTrue_dro, [], [], ISL_noise);
    
            % Perform UKF estimation with measurement update
            [xEst_sst, PEst_sst, NIS_sst] = UKF_ISL(xEst_sst, PEst_sst, delta, S_i_sst_filt, R_ISL1_filt, dt_p, ...
                xEst_dro, PEst_dro, [], [], [], [], alpha_sst_p, beta_sst_p, kappa_sst_p, alpha_sst_m, beta_sst_m, kappa_sst_m);
            last_ISL_measure_time = time;

        else
            % Perform UKF state propagation without measurement update
            [xEst_sst, PEst_sst, ~] = UKF(xEst_sst, PEst_sst, [], S_i_elfo_filt, [], dt_p, [], [], alpha, beta);
        end

        NEES_sst = (xTrue_elfo-xEst_sst)' * PEst_sst^(-1) * (xTrue_elfo-xEst_sst);
    
        %%
        % Update history arrays
        dro1_xDR(n, :) = xDR_dro1';
        dro1_xTrue(n, :) = xTrue_dro1';
        dro1_xEst(n, :) = xEst_dro1';
        dro1_ellipsoid_r(n) = 4*pi/3*sqrt(det(PEst_dro1(1:3,1:3)*(LU*1e3)^2));
        dro1_ellipsoid_v(n) = 4*pi/3*sqrt(det(PEst_dro1(4:6,4:6)*(VU*1e6)^2));
        dro1_NEES(n) = NEES_dro1;
        dro1_NIS(n) = NIS_dro1;

        dro2_xDR(n, :) = xDR_dro2';
        dro2_xTrue(n, :) = xTrue_dro2';
        dro2_xEst(n, :) = xEst_dro2';
        dro2_ellipsoid_r(n) = 4*pi/3*sqrt(det(PEst_dro2(1:3,1:3)*(LU*1e3)^2));
        dro2_ellipsoid_v(n) = 4*pi/3*sqrt(det(PEst_dro2(4:6,4:6)*(VU*1e6)^2));
        dro2_NEES(n) = NEES_dro2;
        dro2_NIS(n) = NIS_dro2;

        dro3_xDR(n, :) = xDR_dro3';
        dro3_xTrue(n, :) = xTrue_dro3';
        dro3_xEst(n, :) = xEst_dro3';
        dro3_ellipsoid_r(n) = 4*pi/3*sqrt(det(PEst_dro3(1:3,1:3)*(LU*1e3)^2));
        dro3_ellipsoid_v(n) = 4*pi/3*sqrt(det(PEst_dro3(4:6,4:6)*(VU*1e6)^2));
        dro3_NEES(n) = NEES_dro3;
        dro3_NIS(n) = NIS_dro3;

        elfo_xDR(n, :) = xDR_elfo';
        elfo_xTrue(n, :) = xTrue_elfo';
        elfo_xEst(n, :) = xEst_elfo';
        elfo_ellipsoid_r(n) = 4*pi/3*sqrt(det(PEst_elfo(1:3,1:3)*(LU*1e3)^2));
        elfo_ellipsoid_v(n) = 4*pi/3*sqrt(det(PEst_elfo(4:6,4:6)*(VU*1e6)^2));
        elfo_NEES(n) = NEES_elfo;
        elfo_NIS(n) = NIS_elfo;

        sst_xEst(n, :) = xEst_sst';
        sst_ellipsoid_r(n) = 4*pi/3*sqrt(det(PEst_sst(1:3,1:3)*(LU*1e3)^2));
        sst_ellipsoid_v(n) = 4*pi/3*sqrt(det(PEst_sst(4:6,4:6)*(VU*1e6)^2));
        sst_NEES(n) = NEES_sst;
        sst_NIS(n) = NIS_sst;

    end
    
    %% ------------ Plotting Calculations ----------------
    % Error related computations for DRO 1
    dro1_errors(:,:,k) = dro1_xTrue - dro1_xEst;
    dro1_rErrorsNorm(:,:,k) = vecnorm((dro1_xTrue(:,1:3) - dro1_xEst(:,1:3)) * LU, 2, 2); % km
    dro1_vErrorsNorm(:,:,k) = vecnorm((dro1_xTrue(:,4:6) - dro1_xEst(:,4:6)) * VU * 1e3, 2, 2); % m/s
    dro1_NEES_M(:,k) = dro1_NEES;
    dro1_NIS_M(:,k) = dro1_NIS;
    
    % Error related computations for DRO 2
    dro2_errors(:,:,k) = dro2_xTrue - dro2_xEst;
    dro2_rErrorsNorm(:,:,k) = vecnorm((dro2_xTrue(:,1:3) - dro2_xEst(:,1:3)) * LU, 2, 2);
    dro2_vErrorsNorm(:,:,k) = vecnorm((dro2_xTrue(:,4:6) - dro2_xEst(:,4:6)) * VU * 1e3, 2, 2);
    dro2_NEES_M(:,k) = dro2_NEES;
    dro2_NIS_M(:,k) = dro2_NIS;
    
    % Error related computations for DRO 3
    dro3_errors(:,:,k) = dro3_xTrue - dro3_xEst;
    dro3_rErrorsNorm(:,:,k) = vecnorm((dro3_xTrue(:,1:3) - dro3_xEst(:,1:3)) * LU, 2, 2);
    dro3_vErrorsNorm(:,:,k) = vecnorm((dro3_xTrue(:,4:6) - dro3_xEst(:,4:6)) * VU * 1e3, 2, 2);
    dro3_NEES_M(:,k) = dro3_NEES;
    dro3_NIS_M(:,k) = dro3_NIS;
    
    % Error related computations for ELFO
    elfo_errors(:,:,k) = elfo_xTrue - elfo_xEst;
    elfo_rErrorsNorm(:,:,k) = vecnorm((elfo_xTrue(:,1:3) - elfo_xEst(:,1:3)) * LU, 2, 2);
    elfo_vErrorsNorm(:,:,k) = vecnorm((elfo_xTrue(:,4:6) - elfo_xEst(:,4:6)) * VU * 1e3, 2, 2);
    elfo_NEES_M(:,k) = elfo_NEES;
    elfo_NIS_M(:,k) = elfo_NIS;
    
    % Error related computations for ELFO using ISL
    sst_errors(:,:,k) = elfo_xTrue - sst_xEst;
    sst_rErrorsNorm(:,:,k) = vecnorm((elfo_xTrue(:,1:3) - sst_xEst(:,1:3)) * LU, 2, 2);
    sst_vErrorsNorm(:,:,k) = vecnorm((elfo_xTrue(:,4:6) - sst_xEst(:,4:6)) * VU * 1e3, 2, 2);
    sst_NEES_M(:,k) = sst_NEES;
    sst_NIS_M(:,k) = sst_NIS;

    dro1_ellipsoid_r_k(:,k) = dro1_ellipsoid_r;
    dro1_ellipsoid_v_k(:,k) = dro1_ellipsoid_v;
    dro2_ellipsoid_r_k(:,k) = dro2_ellipsoid_r;
    dro2_ellipsoid_v_k(:,k) = dro2_ellipsoid_v;
    dro3_ellipsoid_r_k(:,k) = dro3_ellipsoid_r;
    dro3_ellipsoid_v_k(:,k) = dro3_ellipsoid_v;
    elfo_ellipsoid_r_k(:,k) = elfo_ellipsoid_r;
    elfo_ellipsoid_v_k(:,k) = elfo_ellipsoid_v;
    sst_ellipsoid_r_k(:,k) = sst_ellipsoid_r;
    sst_ellipsoid_v_k(:,k) = sst_ellipsoid_v;

    elfo_visArrays(:,1,k) = elfo_visSlice;
    dro1_visArrays(:,1,k) = dro1_visSlice;
    sst_visArrays(:,1,k) = sst_visSlice;

    elfo_arcArrays(:,1,k) = elfo_arcSlice;
    dro1_arcArrays(:,1,k) = dro1_arcSlice;
    sst_arcArrays(:,1,k) = sst_arcSlice;
end

%% ------------ Storage ----------------
% Get the mean by Monte Carlo method
dro1.errors = sqrt(mean(dro1_errors.^2, 3));
dro1.rErrorsNorm = sqrt(mean(dro1_rErrorsNorm.^2, 3)); % In km
dro1.vErrorsNorm = sqrt(mean(dro1_vErrorsNorm.^2, 3)); % In m/s
% Get the mean by Monte Carlo method
dro2.errors = sqrt(mean(dro2_errors.^2, 3));
dro2.rErrorsNorm = sqrt(mean(dro2_rErrorsNorm.^2, 3)); % In km
dro2.vErrorsNorm = sqrt(mean(dro2_vErrorsNorm.^2, 3)); % In m/s
% Get the mean by Monte Carlo method
dro3.errors = sqrt(mean(dro3_errors.^2, 3));
dro3.rErrorsNorm = sqrt(mean(dro3_rErrorsNorm.^2, 3)); % In km
dro3.vErrorsNorm = sqrt(mean(dro3_vErrorsNorm.^2, 3)); % In m/s
% Get the mean by Monte Carlo method
elfo.errors = sqrt(mean(elfo_errors.^2, 3));
elfo.rErrorsNorm = sqrt(mean(elfo_rErrorsNorm.^2, 3)); % In km
elfo.vErrorsNorm = sqrt(mean(elfo_vErrorsNorm.^2, 3)); % In m/s
% Get the mean by Monte Carlo method
elfo.errors_ISL = sqrt(mean(sst_errors.^2, 3));
elfo.rErrorsNorm_ISL = sqrt(mean(sst_rErrorsNorm.^2, 3)); % In km
elfo.vErrorsNorm_ISL = sqrt(mean(sst_vErrorsNorm.^2, 3)); % In m/s

% Ellipsoid volumnes
dro1.Ellipsoid_r = sqrt(mean(dro1_ellipsoid_r_k.^2, 2));
dro1.Ellipsoid_v = sqrt(mean(dro1_ellipsoid_v_k.^2, 2));
dro2.Ellipsoid_r = sqrt(mean(dro2_ellipsoid_r_k.^2, 2));
dro2.Ellipsoid_v = sqrt(mean(dro2_ellipsoid_v_k.^2, 2));
dro3.Ellipsoid_r = sqrt(mean(dro3_ellipsoid_r_k.^2, 2));
dro3.Ellipsoid_v = sqrt(mean(dro3_ellipsoid_v_k.^2, 2));
elfo.Ellipsoid_r = sqrt(mean(elfo_ellipsoid_r_k.^2, 2));
elfo.Ellipsoid_v = sqrt(mean(elfo_ellipsoid_v_k.^2, 2));
elfo.ISL_Ellipsoid_r = sqrt(mean(sst_ellipsoid_r_k.^2, 2));
elfo.ISL_Ellipsoid_v = sqrt(mean(sst_ellipsoid_v_k.^2, 2));

% NEES
dro1.NEES = mean(dro1_NEES_M, 2);
dro2.NEES = mean(dro2_NEES_M, 2);
dro3.NEES = mean(dro3_NEES_M, 2);
elfo.NEES = mean(elfo_NEES_M, 2);
elfo.NEES_sst = mean(sst_NEES_M, 2);

% NIS
dro1.NIS = mean(dro1_NIS_M, 2);
dro2.NIS = mean(dro2_NIS_M, 2);
dro3.NIS = mean(dro3_NIS_M, 2);
elfo.NIS = mean(elfo_NIS_M, 2);
elfo.NIS_sst = mean(sst_NIS_M, 2);

dro1.visArrays = dro1_visArrays;
dro1.measArrays = dro1_arcArrays;
elfo.visArrays = elfo_visArrays;
elfo.measArrays = elfo_arcArrays;
elfo.ISLvisArrays = sst_visArrays;
elfo.ISLmeasArrays = sst_arcArrays;

if nRelay == 1
    save('output/dro1(1)','dro1');
    save('output/elfo(1)','elfo');
elseif nRelay == 2
    save('output/dro1(2)','dro1');
    save('output/dro2(2)','dro2');
    save('output/elfo(2)','elfo');
elseif nRelay == 3
    save('output/dro1(3)','dro1');
    save('output/dro2(3)','dro2');
    save('output/dro3(3)','dro3');
    save('output/elfo(3)','elfo');
end

n = size(elfo.NEES_sst, 1);
last_fifth = ceil(n / 3 * 2);
start_index = n - last_fifth + 1;
last_elements = elfo.NEES_sst(start_index:n);
avgNEES(i) = mean(last_elements);

n = size(elfo.NIS_sst, 1);
last_fifth = ceil(n / 3 * 2);
start_index = n - last_fifth + 1;
last_elements = elfo.NIS_sst(start_index:n);
avgNIS(i) = mean(last_elements);
end

toc;

function [xEst, PEst, NIS] = UKF(xEst, PEst, y, Pv, Pn, dt_p, rGS, vGS, alpha, beta)

    switch isempty(y)

        case true

            % Dimensions
            L_x = numel(xEst); % Number of dimensions of states
            L_v = size(Pv, 1); % Number of dimensions of process noise
            L_a = L_x + L_v; % Dimension of augmented state
        
            % UKF parameters
            kappa = 3-L_a; % Secondary scaling parameter
        
            % Scaling parameters
            lambda = alpha^2 * (L_a + kappa) - L_a; % Scaling factor
            gamma = sqrt(L_a + lambda); % Sigma point spread
        
            % Weight factors
            Wm = [lambda / (L_a + lambda); repmat(1 / (2 * (L_a + lambda)), 2 * L_a, 1)];
            Wc = Wm;
            Wc(1) = Wc(1) + (1 - alpha^2 + beta);
        
            %%%%%%%%%%%%%%%%%%%%%%% Calculate sigma points %%%%%%%%%%%%%%%%%%%%%%%%
            % Augment state and covariance
            xEst_a = [xEst; zeros(L_v, 1)]; % Augmented state vector
            P_a = blkdiag(PEst, Pv); % Augmented covariance matrix
        
            % Generate sigma points
            sigma_points = zeros(L_a, 2 * L_a + 1);
            sigma_points(:, 1) = xEst_a;
        
            sqrtP_a = chol(P_a, 'lower'); % Cholesky decomposition
            for i = 1:L_a
                sigma_points(:, i + 1) = xEst_a + gamma * sqrtP_a(:, i);
                sigma_points(:, i + L_a + 1) = xEst_a - gamma * sqrtP_a(:, i);
            end
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Time update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Extract state and process noise sigma points
            x_sigma = sigma_points(1:L_x, :); % State sigma points
            v_sigma = sigma_points(L_x + 1:L_x + L_v, :); % Process noise sigma points
        
            % Time update in UKF
            sigma_points_pred = zeros(L_x, 2 * L_a + 1);
            for i = 1:2 * L_a + 1
                % Use v_sigma directly as process noise
                sigma_points_pred(:, i) = motion_model_F(x_sigma(:, i), v_sigma(:, i), dt_p);
            end
        
            % Predict state mean
            xPred = sigma_points_pred * Wm;
        
            % Predict state covariance
            PPred = zeros(L_x, L_x);
            for i = 1:2 * L_a + 1
                diff = sigma_points_pred(:, i) - xPred;
                PPred = PPred + Wc(i) * (diff * diff');
            end
        
            % Update the state and covariance for the next iteration (before measurement)
            xEst = xPred;
            PEst = PPred;
            NIS = [];
        
        case false

            %%%%%%%%%%%%%%%%%%%%%%%%%%% Measurement update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Dimensions
            L_x = numel(xEst); % Number of dimensions of states
            L_v = size(Pv, 1); % Number of dimensions of process noise
            L_n = size(Pn, 1); % Number of dimensions of measurement noise
            L_a = L_x + L_v + L_n; % Dimension of augmented state
        
            % UKF parameters
            kappa = 3-L_a; % Secondary scaling parameter
        
            % Scaling parameters
            lambda = alpha^2 * (L_a + kappa) - L_a; % Scaling factor
            gamma = sqrt(L_a + lambda); % Sigma point spread
        
            % Weight factors
            Wm = [lambda / (L_a + lambda); repmat(1 / (2 * (L_a + lambda)), 2 * L_a, 1)];
            Wc = Wm;
            Wc(1) = Wc(1) + (1 - alpha^2 + beta);
        
            %%%%%%%%%%%%%%%%%%%%%%% Calculate sigma points %%%%%%%%%%%%%%%%%%%%%%%%
            % Augment state and covariance
            xEst_a = [xEst; zeros(L_v, 1); zeros(L_n, 1)]; % Augmented state vector
            P_a = blkdiag(PEst, Pv, Pn); % Augmented covariance matrix
        
            % Generate sigma points
            sigma_points = zeros(L_a, 2 * L_a + 1);
            sigma_points(:, 1) = xEst_a;
        
            sqrtP_a = chol(P_a, 'lower'); % Cholesky decomposition
            for i = 1:L_a
                sigma_points(:, i + 1) = xEst_a + gamma * sqrtP_a(:, i);
                sigma_points(:, i + L_a + 1) = xEst_a - gamma * sqrtP_a(:, i);
            end
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Time update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Extract state and process noise sigma points
            x_sigma = sigma_points(1:L_x, :); % State sigma points
            v_sigma = sigma_points(L_x + 1:L_x + L_v, :); % Process noise sigma points
        
            % Time update in UKF
            sigma_points_pred = zeros(L_x, 2 * L_a + 1);
            for i = 1:2 * L_a + 1
                sigma_points_pred(:, i) = motion_model_F(x_sigma(:, i), v_sigma(:, i), dt_p);
            end
        
            % Predict state mean
            xPred = sigma_points_pred * Wm;
        
            % Predict state covariance
            PPred = zeros(L_x, L_x);
            for i = 1:2 * L_a + 1
                diff = sigma_points_pred(:, i) - xPred;
                PPred = PPred + Wc(i) * (diff * diff');
            end
        
            % Extract measurement noise sigma points
            n_sigma = sigma_points(L_x + L_v + 1:end, :);
    
            % Predict measurement sigma points
            sigma_points_meas = zeros(L_n, 2 * L_a + 1); % L_n = 4 (number of measurements)
            for i = 1:2 * L_a + 1
                % Pass noise vector
                sigma_points_meas(:, i) = observation_model_H(sigma_points_pred(:, i), ...
                    rGS, vGS, n_sigma(:, i)); 
            end
    
            % Predict measurement mean
            yPred = sigma_points_meas * Wm;
    
            % Predict measurement covariance
            Pyy = zeros(L_n, L_n);
            for i = 1:2 * L_a + 1
                diff = sigma_points_meas(:, i) - yPred;
                Pyy = Pyy + Wc(i) * (diff * diff');
            end
    
            % Cross covariance between state and measurement
            Pxy = zeros(L_x, L_n);
            for i = 1:2 * L_a + 1
                diffX = sigma_points_pred(:, i) - xPred;
                diffY = sigma_points_meas(:, i) - yPred;
                Pxy = Pxy + Wc(i) * (diffX * diffY');
            end
    
            % Kalman gain
            K = Pxy / Pyy;
        
            xEst = xPred + K * (y - yPred); % Update state mean
            PEst = PPred - K * Pyy * K'; % Update state covariance
            NIS = (y - yPred)' * Pyy^(-1) * (y - yPred);
    end
end

function [xEst, PEst, NIS] = UKF_ISL( ...
        xEst, PEst, delta, Pv, Pn, dt_p, ...
        xEst_dro1, PEst_dro1, xEst_dro2, PEst_dro2, xEst_dro3, PEst_dro3, ...
        alpha, beta, kappa, alpha_sst, beta_sst, kappa_sst)

%% ------------------------- Dimensions ----------------------------------
Lx = numel(xEst);      % state
Lv = size(Pv,1);       % process noise
Ln = size(Pn,1);       % meas noise
La = Lx + Lv + Ln;     % augmented

%% ------------------------- Prediction-UT parameters ---------------------
lambda_p = alpha^2 * (La + kappa) - La;
gamma_p  = sqrt(La + lambda_p);

Wm_p = [lambda_p /(La+lambda_p);  repmat(1/(2*(La+lambda_p)),2*La,1)];
Wc_p = Wm_p;  Wc_p(1) = Wc_p(1) + (1-alpha^2+beta);

%% ------------------------- Augmented sigma points ----------------------
xAug = [xEst; zeros(Lv,1); zeros(Ln,1)];
PAug = blkdiag(PEst, Pv, Pn);

S = chol(PAug,'lower');
sig = zeros(La,2*La+1);
sig(:,1) = xAug;
for i = 1:La
    sig(:,i+1)     = xAug + gamma_p*S(:,i);
    sig(:,i+La+1)  = xAug - gamma_p*S(:,i);
end

%% ------------------------- Time update ---------------------------------
x_sig = sig(1:Lx,:);     % state part
v_sig = sig(Lx+1:Lx+Lv,:);

xPred_sig = zeros(Lx,2*La+1);
for i = 1:2*La+1
    xPred_sig(:,i) = motion_model_F(x_sig(:,i), v_sig(:,i), dt_p);
end
xPred = xPred_sig * Wm_p;

PPred = zeros(Lx);
for i = 1:2*La+1
    d = xPred_sig(:,i) - xPred;
    PPred = PPred + Wc_p(i)*(d*d.');
end

%% ------------------------- Measurement update --------------------------
switch size(delta,1)

case 2
    z_sig = zeros(Ln,2*La+1);
    n_sig = sig(Lx+Lv+1:end,:);
    for i = 1:2*La+1
        z_sig(:,i) = ISL_model(xPred_sig(:,i), xEst_dro1, [], [], n_sig(:,i));
    end

    H  = build_H_SST(xPred, xEst_dro1);
    R_eff = H * PEst_dro1 * H.';

case 4 
    z_sig = zeros(Ln,2*La+1);
    n_sig = sig(Lx+Lv+1:end,:);
    for i = 1:2*La+1
        z_sig(:,i) = ISL_model(xPred_sig(:,i), xEst_dro1, xEst_dro2, [], n_sig(:,i));
    end

    H1 = build_H_SST(xPred, xEst_dro1);
    H2 = build_H_SST(xPred, xEst_dro2);
    H  = [H1; H2];                          

    R_eff = zeros(4);
    R_eff([1 3],[1 3]) = H1*PEst_dro1*H1.';
    R_eff([2 4],[2 4]) = H2*PEst_dro2*H2.';

case 6  
    z_sig = zeros(Ln,2*La+1);
    n_sig = sig(Lx+Lv+1:end,:);
    for i = 1:2*La+1
        z_sig(:,i) = ISL_model(xPred_sig(:,i), xEst_dro1, xEst_dro2, xEst_dro3, n_sig(:,i));
    end

    H1 = build_H_SST(xPred, xEst_dro1);
    H2 = build_H_SST(xPred, xEst_dro2);
    H3 = build_H_SST(xPred, xEst_dro3);
    H  = [H1; H2; H3];                          

    R_eff = zeros(6);
    R_eff([1 4],[1 4]) = H1*PEst_dro1*H1.';
    R_eff([2 5],[2 5]) = H2*PEst_dro2*H2.';
    R_eff([3 6],[3 6]) = H3*PEst_dro3*H3.';

otherwise
    error('Unsupported measurement dimension.');
end

%% -------- Measurement-UT parameters (fixed α_sst,β_sst,κ_sst) ----------
lambda_m = alpha_sst^2 * (La + kappa_sst) - La;
gamma_m  = sqrt(La + lambda_m); 
Wm_m = [lambda_m/(La+lambda_m); repmat(1/(2*(La+lambda_m)),2*La,1)];
Wc_m = Wm_m;  Wc_m(1) = Wc_m(1) + (1-alpha_sst^2+beta_sst);

%% -------- Measurement statistics --------------------------------------
zPred = z_sig * Wm_m;

Pyy = zeros(Ln);           Pxy = zeros(Lx,Ln);
for i = 1:2*La+1
    dz = z_sig(:,i) - zPred;
    dx = xPred_sig(:,i) - xPred;
    Pyy = Pyy + Wc_m(i)*(dz*dz.');
    Pxy = Pxy + Wc_m(i)*(dx*dz.');
end
Pyy = Pyy + R_eff;

%% -------- Kalman update ------------------------------------------------
K    = Pxy / Pyy;
xEst = xPred + K*(delta - zPred);
PEst = PPred - K*Pyy*K.';
NIS  = (delta - zPred).' / Pyy * (delta - zPred);

end

function x = motion_model_F(x, process_noise, dt)

    r = x(1:3);
    v = x(4:6);
    mu = 1.215058560962404E-2;

    % Define the ODE function for the system dynamics
    function dx = dynamics(r, v)
        r1 = sqrt((r(1) + mu)^2 + r(2)^2 + r(3)^2);
        r2 = sqrt((r(1) - 1 + mu)^2 + r(2)^2 + r(3)^2);
        a = [r(1) + 2 * v(2) - (1 - mu) * (r(1) + mu) / r1^3 - mu * (r(1) - 1 + mu) / r2^3;
             r(2) - 2 * v(1) - (1 - mu) * r(2) / r1^3 - mu * r(2) / r2^3;
             -(1 - mu) * r(3) / r1^3 - mu * r(3) / r2^3];
        dx = [v; a];
    end

    % Perform Runge-Kutta 4th order integration
    k1 = dynamics(r, v);
    k2 = dynamics(r + 0.5 * dt * k1(1:3), v + 0.5 * dt * k1(4:6));
    k3 = dynamics(r + 0.5 * dt * k2(1:3), v + 0.5 * dt * k2(4:6));
    k4 = dynamics(r + dt * k3(1:3), v + dt * k3(4:6));
    x = x + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4);

    % Add process noise if provided
    if ~isempty(process_noise)
        x = x + process_noise;
    end
end

function z = observation_model_H(xSat, rGS, vGS, observe_noise)

    rho = xSat(1:3) - rGS;        
    range = norm(rho);

    elev = asin(rho(3)/range);
    azim = atan2(rho(2),rho(1));

    range_rate = dot(rho , (xSat(4:6) - vGS)) / range;

    z = [range; elev; azim; range_rate] + observe_noise;
end

function Z_xyz = calculate_coordinates(Z, rGS, vGS)

    rho      = Z(1);
    elev   = Z(2);
    az     = Z(3);
    rho_dot  = Z(4);

    u_los = [ cos(elev)*cos(az);
              cos(elev)*sin(az);
              sin(elev)           ];

    Z_xyz(1:3,1) = rGS(:) + rho * u_los;

    Z_xyz(4:6,1) = vGS(:) + rho_dot * u_los;
end

function [delta] = ISL_model(target, source1, source2, source3, noise)

    if isempty(source2) && isempty (source3)

        delta(1, :) = norm(target(1:3) - source1(1:3)) + noise(1); % Range
    
        dr  = target(1:3) - source1(1:3);
        dv  = target(4:6) - source1(4:6);
        rho_dot = dot(dr, dv)/norm(dr);
        delta(2, :) = rho_dot + noise(2); % Speed

    elseif isempty(source3)

        % Compute range and range rate for two sources relative to the target
        delta(1, :) = norm(target(1:3) - source1(1:3)) + noise(1); % Range 1
        delta(2, :) = norm(target(1:3) - source2(1:3)) + noise(2); % Range 2
    
        dr  = target(1:3) - source1(1:3);
        dv  = target(4:6) - source1(4:6);
        rho_dot = dot(dr, dv)/norm(dr);
        delta(3, :) = rho_dot + noise(3); % Speed value
    
        dr  = target(1:3) - source2(1:3);
        dv  = target(4:6) - source2(4:6);
        rho_dot = dot(dr, dv)/norm(dr);
        delta(4, :) = rho_dot + noise(4); % Speed value

    else

        % Range and range rate data from 3 source spacecrafts to the target
        delta(1, :) = norm(target(1:3) - source1(1:3)) + noise(1); % Range 1
        delta(2, :) = norm(target(1:3) - source2(1:3)) + noise(2); % Range 2
        delta(3, :) = norm(target(1:3) - source3(1:3)) + noise(3); % Range 3
    
        dr  = target(1:3) - source1(1:3);
        dv  = target(4:6) - source1(4:6);
        rho_dot = dot(dr, dv)/norm(dr);
        delta(4, :) = rho_dot + noise(4);
    
        dr  = target(1:3) - source2(1:3);
        dv  = target(4:6) - source2(4:6);
        rho_dot = dot(dr, dv)/norm(dr);
        delta(5, :) = rho_dot + noise(5);
        
        dr  = target(1:3) - source3(1:3);
        dv  = target(4:6) - source3(4:6);
        rho_dot = dot(dr, dv)/norm(dr);
        delta(6, :) = rho_dot + noise(6);

    end
end 

function rCR3BP = dsn_pos_CR3BP(iStation,t,DSN)

    LU = 389703;
    TU = 382981; 
    mu = 1.215058560962404E-2;

    a  = 6378.137;        
    e2 = 0.00669437999014; 
    lat = deg2rad(DSN.lat(iStation));
    lon = deg2rad(DSN.lon(iStation));
    h   = DSN.alt(iStation);

    N = a / sqrt(1 - e2*sin(lat)^2);
    xE = (N+h)*cos(lat)*cos(lon);
    yE = (N+h)*cos(lat)*sin(lon);
    zE = (N*(1-e2)+h)*sin(lat);
    rECEF0 = [xE; yE; zE];      

    theta_E  = DSN.omega * t * TU;
    Rz_eci   = [cos(-theta_E) -sin(-theta_E) 0;
                sin(-theta_E)  cos(-theta_E) 0;
                     0               0       1];
    rECEF = Rz_eci * rECEF0;   

    rCR3BP = [-mu;0;0] + rECEF/LU;
end

function idx = choose_station(xSat, time, DSN, elMask)

    mu = 1.215058560962404E-2;

    elMask = deg2rad(elMask);     
    elev   = -Inf(1,3);           

    for k = 1:3
        rGS  = dsn_pos_CR3BP(k, time, DSN);
        rho  = xSat(1:3) - rGS;   
        range = norm(rho);

        zHat = -(rGS + [mu;0;0])/norm(rGS + [mu;0;0]);
        elev(k) = asin( dot(rho/range, zHat) );
    end

    visible = elev >= elMask;

    if any(visible)
        [~, idx] = max(elev .* visible);
    else
        idx = 0;  
    end
end

function H = build_H_SST(x_nav, x_ref)

    dr   = x_nav(1:3) - x_ref(1:3);
    dv   = x_nav(4:6) - x_ref(4:6);
    rho  = norm(dr);
    rho_d= dot(dr,dv)/rho;
    u_r  = dr.'/rho;

    H11 = -u_r;
    H12 = zeros(1,3);
    H21 = -(dv.'/rho) + (rho_d/rho^2)*u_r;
    H22 = -u_r;

    H = [H11 H12;
         H21 H22];
end
