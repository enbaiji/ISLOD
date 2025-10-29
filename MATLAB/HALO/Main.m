% Author: Enbai Ji @ UNSW.
% Email: enbaiji@outlook.com.
% All rights reserved.

%%%%%%%%%%%%%%%%%%% Initialisation of Environment/Variables %%%%%%%%%%%%%%%%%%%
tic;
clear;
ts = cputime;
format long;
rng(0);
N = 1;

addpath("mice/lib");
addpath("mice/src/mice");
addpath('prop');
addpath('input');
addpath('output');

load('ORBdataTrue1.mat');
orbTrue1 = orb;
clear orb;

load('ORBdataTrue2.mat');
orbTrue2 = orb;
clear orb;

% Low fidelity model for estimation
orb        = struct();
orb.model  = "Estimation";
orb.type   = "halo";
orbEs1 = LoadSequential(orb);

orb         = struct();
orb.model   = "Estimation";
orb.type    = "elfo2";
orbEst2 = LoadSequential(orb);

%% General Settings
options      = odeset('RelTol',1e-6,'AbsTol',1e-9); % The option for the ODE integrater
tStep        = 60;                                  % Minimum time step (default: 60s)
measIntvBase = 60;                                  % Baseline masurement interval (default 60s)
dtGB     = tStep * 10;                          % Selected GB measurement update interval
dtSST    = tStep * 10;                          % Selected SST measurement update interval
iniGBArc     = 172800;                            % Initial ground arc (default: 2d)

%% Process Noise
% Acceleration PSD [km^2/s^3]
accNoise_DRO  = 4e-16;           
accNoise_ELFO = 4e-16;         

% Process noise of clock error (PSD)
q_DRO  = orbTrue1.noise.clcNoise;        % (clcNoise = 1e-28)
q_ELFO = orbTrue2.noise.clcNoise;       % (clcNoise = 1e-28)
q_SST  = orbTrue2.noise.clcNoise;       % (clcNoise = 1e-28)

% Process noise
Q_DRO  = diag([accNoise_DRO, accNoise_DRO, accNoise_DRO]);  
Q_ELFO = diag([accNoise_ELFO, accNoise_ELFO, accNoise_ELFO]);  

Q_DRO_SST  = Q_DRO;  
Q_ELFO_SST = Q_ELFO;  

%% Measurement Noise
% Baseline measurement noise
R_GB = diag([ 1e-3, ...         % 1 m  
              5e-6, ...         % 1 arcsec
              5e-6, ...         % 1 arcsec
              5e-8].^2);        % 0.05 mm/s (Δf ≈ 0.1 Hz)

% R_SST1 = 0.001^2;                                       % [m]
R_SST2 = diag([1e-6*5, 0.01e-6*5].^2);                    % [mm, mm/s]
% R_SST4 = diag([0.001, 50E-6, 50E-6, 0.1e-6].^2);        

% Measurement noise factor for time integration (k=dt*/dt)
k_GB  = measIntvBase / dtGB;
k_SST = measIntvBase / dtSST;

% Integrated measurement noise
R_GB(4,4)   = R_GB(4,4)    * k_GB^2;
R_SST2(2,2) = R_SST2(2,2)  * k_SST^2;
% R_SST4(4,4) = R_SST4(4,4)  * k_SST^2;

%% Initial States and Covariances
% Initial standard deviation of clock error
Pc0_DRO  = diag([1e-9, 1e-12]);
Pc0_ELFO = diag([1e-9, 1e-12]);

% Augmented initial states and covariance matrices
xEst0_DRO  = [orbEs1.sat.X0iner' + [1*randn(3,1); 1e-3*randn(3,1)]; ...
    orbTrue1.sat.b0' + [1e-9*randn(1,1); 1e-12*randn(1,1)]];

PEst0_DRO  = blkdiag(eye(3), eye(3) * 1e-3, Pc0_DRO).^2;

xEst0_ELFO = [orbEst2.sat.X0iner' + [1*randn(3,1); 1e-3*randn(3,1)]; ...
    orbTrue2.sat.b0' + [1e-9*randn(1,1); 1e-12*randn(1,1)]];

PEst0_ELFO = blkdiag(eye(3), eye(3) * 1e-3, Pc0_ELFO).^2;

%% UKF Parameters
% DRO GB
estParaDRO   = [0.25, 2, 0];

% ELFO GB
estParaELFO  = [0.1, 2, 0];

% SST
predParaDRO  = [0.25, 2, 0];
predParaELFO = [0.1, 2, 0];
estParaSST   = [0.1, 2, 0];

%%
%%%%%%%%%%%%%%%%%%% Main Function %%%%%%%%%%%%%%%%%%%

span      = orbEs1.seq.a.span;
t0        = orbEs1.seq.Time;
tEnd      = t0 + span;

% Pre-allocate space for ephemeris array
nStep      = floor(span/tStep) + 1;

MarrayT             = zeros(nStep, 1, N);
MarrayX_DRO         = zeros(nStep, 8, N);
MarrayX_ELFO        = zeros(nStep, 8, N);
MarrayNEES_DRO_acc  = zeros(nStep, 1, N);
MarrayNEES_ELFO_acc = zeros(nStep, 1, N);
MarrayNEES_DRO_clc  = zeros(nStep, 1, N);
MarrayNEES_ELFO_clc = zeros(nStep, 1, N);
MarrayNIS_DRO       = zeros(nStep, 1, N);
MarrayNIS_ELFO      = zeros(nStep, 1, N);

for n = 1:N

    % SPICE initialisations
    cspice_kclear;
    metakernelcheck;
    cspice_furnsh('metakernel.tm');

    rng(n);

    arrayT             = zeros(nStep, 1);
    arrayX_DRO         = zeros(nStep, 8);
    arrayX_ELFO        = zeros(nStep, 8);
    arrayNEES_DRO_acc  = zeros(nStep, 1);
    arrayNEES_ELFO_acc = zeros(nStep, 1);
    arrayNEES_DRO_clc  = zeros(nStep, 1);
    arrayNEES_ELFO_clc = zeros(nStep, 1);
    arrayNIS_DRO       = zeros(nStep, 1);
    arrayNIS_ELFO      = zeros(nStep, 1);

    % Store the initial ephemeris
    xEst_DRO  = xEst0_DRO;
    xEst_ELFO = xEst0_ELFO;
    PEst_DRO  = PEst0_DRO;
    PEst_ELFO = PEst0_ELFO;

    arrayT(1)         = t0;
    arrayX_DRO(1, :)  = xEst_DRO';
    arrayX_ELFO(1, :) = xEst_ELFO';

    t        = t0;
    % tLastSST = t0 + iniGBArc;       % SST starts from the end of GB arc
    P12 = [];
    
    % Main loop
    while t < tEnd

            if t < (t0 + iniGBArc)
                dt = min(dtGB,  tEnd - t);
            else
                dt = min(dtSST, tEnd - t);
            end
             
            tIndexMeas = floor((t+dt - t0)/tStep)+1; % Measurement epoch
    
            % Extract true states
            xTrue_DRO  = orbTrue1.seq.a.XJ2000(tIndexMeas, :)';
            xTrue_ELFO = orbTrue2.seq.a.XJ2000(tIndexMeas, :)';
            
            % if within the initial GB arc, use ground measurements
            if t < (t0 + iniGBArc)

                fprintf('In GB arc:\n');

                % Propagate for interpolations
                [tTemp_DRO, Xtemp_DRO] = ode113(@(t,x) prophpopClock(t,x,orbEs1), ...
                    t : tStep : t+dt, xEst_DRO, options);
                [~, Xtemp_ELFO] = ode113(@(t,x) prophpopClock(t,x,orbEst2), ...
                    t : tStep : t+dt, xEst_ELFO, options);
        
                % Interpolations    
                timeInterpo = tTemp_DRO(2:end-1);
                rowsInterpo = numel(timeInterpo);
                arrayT     (tIndexMeas-rowsInterpo : tIndexMeas-1)    = tTemp_DRO (2:end-1);
                arrayX_DRO (tIndexMeas-rowsInterpo : tIndexMeas-1, :) = Xtemp_DRO (2:end-1, :);
                arrayX_ELFO(tIndexMeas-rowsInterpo : tIndexMeas-1, :) = Xtemp_ELFO(2:end-1, :);

                % DRO process noise
                Pc_DRO = [(dt^3)/3, (dt^2)/2; (dt^2)/2, dt] * q_DRO;
                P_DRO  = blkdiag([Q_DRO*dt^3/3, Q_DRO*dt^2/2; Q_DRO*dt^2/2, Q_DRO*dt], Pc_DRO);
       
                % ELFO process noise
                Pc_ELFO = [(dt^3)/3, (dt^2)/2; (dt^2)/2, dt] * q_ELFO;
                P_ELFO  = blkdiag([Q_ELFO*dt^3/3, Q_ELFO*dt^2/2; Q_ELFO*dt^2/2, Q_ELFO*dt], Pc_ELFO);

                zMeas_DRO = MeasurementModel(xTrue_DRO, t+dt, 'Earth', []);
                noiseVec  = chol(R_GB) * randn(length(zMeas_DRO), 1);
                zMeas_DRO = zMeas_DRO + noiseVec;

                % DRO ground measurement update
                [xEst_DRO, PEst_DRO, NIS_DRO] = UKF( ...
                    xEst_DRO, PEst_DRO, ...    
                    P_DRO, R_GB, ...
                    t, dt, orbEs1, ...
                    zMeas_DRO, 'Earth', [], [], ... 
                    estParaDRO, options);

                zMeas_ELFO = MeasurementModel(xTrue_ELFO, t+dt, 'Earth', []);
                noiseVec   = chol(R_GB) * randn(length(zMeas_ELFO), 1);
                zMeas_ELFO = zMeas_ELFO + noiseVec;

                % ELFO ground measurement update
                [xEst_ELFO, PEst_ELFO, NIS_ELFO] = UKF( ...
                    xEst_ELFO, PEst_ELFO, ...
                    P_ELFO, R_GB, ...
                    t, dt, orbEst2, ...
                    zMeas_ELFO, 'Earth', [], [], ...
                    estParaELFO, options);
    
            % if exceed the initial GB arc but not exceed the SST interval
            else

                fprintf('In SST arc\n');
                
                % Propagate for interpolations
                [tTemp_DRO, Xtemp_DRO] = ode113(@(t,x) prophpopClock(t,x,orbEs1), ...
                    t : tStep : t+dt, xEst_DRO, options);
                [~, Xtemp_ELFO] = ode113(@(t,x) prophpopClock(t,x,orbEst2), ...
                    t : tStep : t+dt, xEst_ELFO, options);
        
                % Interpolations    
                timeInterpo = tTemp_DRO(2:end-1);
                rowsInterpo = numel(timeInterpo);
                arrayT     (tIndexMeas-rowsInterpo : tIndexMeas-1)    = tTemp_DRO (2:end-1);
                arrayX_DRO (tIndexMeas-rowsInterpo : tIndexMeas-1, :) = Xtemp_DRO (2:end-1, :);
                arrayX_ELFO(tIndexMeas-rowsInterpo : tIndexMeas-1, :) = Xtemp_ELFO(2:end-1, :);

                Pc_DRO = [(dt^3)/3, (dt^2)/2; (dt^2)/2, dt] * q_DRO;
                P_DRO  = blkdiag([Q_DRO_SST*dt^3/3, Q_DRO_SST*dt^2/2; Q_DRO_SST*dt^2/2, Q_DRO_SST*dt], Pc_DRO);

                Pc_ELFO = [(dt^3)/3, (dt^2)/2; (dt^2)/2, dt] * q_ELFO;
                P_ELFO  = blkdiag([Q_ELFO_SST*dt^3/3, Q_ELFO_SST*dt^2/2; Q_ELFO_SST*dt^2/2, Q_ELFO_SST*dt], Pc_ELFO);
    
                yMeas_SST = MeasurementModel(xTrue_ELFO, t+dt, 'SST2', xTrue_DRO);
                noiseVec  = chol(R_SST2) * randn(length(yMeas_SST), 1);
                yMeas_SST = yMeas_SST + noiseVec;
    
                [xEst_ELFO, PEst_ELFO, xEst_DRO, PEst_DRO, P12, NIS] = JUKF( ...
                    xEst_ELFO, PEst_ELFO, xEst_DRO, PEst_DRO, P12, ...    
                    P_ELFO, P_DRO, R_SST2, ...
                    t, dt, orbEst2, orbEs1, ...
                    yMeas_SST, "SST2", ...   
                    predParaELFO, predParaDRO, estParaSST, options);

                NIS_DRO = NIS;
                NIS_ELFO = NIS;
                
            end
  
            %%
            % Store estimations to ephemeris
            arrayT     (tIndexMeas)    = t+dt; % Measurement epoch
            arrayX_DRO (tIndexMeas, :) = xEst_DRO';
            arrayX_ELFO(tIndexMeas, :) = xEst_ELFO';
    
            diffxDRO  = xTrue_DRO - xEst_DRO;
            diffxELFO = xTrue_ELFO - xEst_ELFO;
    
            NEES_DRO_acc  = diffxDRO(1:6)' * (PEst_DRO(1:6,1:6) \ diffxDRO(1:6));
            NEES_ELFO_acc = diffxELFO(1:6)' * (PEst_ELFO(1:6,1:6) \ diffxELFO(1:6));
            NEES_DRO_clc  = diffxDRO(7:8)' * (PEst_DRO(7:8,7:8) \ diffxDRO(7:8));
            NEES_ELFO_clc = diffxELFO(7:8)' * (PEst_ELFO(7:8,7:8) \ diffxELFO(7:8));
    
            arrayNEES_DRO_acc(tIndexMeas)  = NEES_DRO_acc;
            arrayNEES_ELFO_acc(tIndexMeas) = NEES_ELFO_acc;
            arrayNEES_DRO_clc(tIndexMeas)  = NEES_DRO_clc;
            arrayNEES_ELFO_clc(tIndexMeas) = NEES_ELFO_clc;

            arrayNIS_DRO(tIndexMeas)  = NIS_DRO;
            arrayNIS_ELFO(tIndexMeas) = NIS_ELFO;
    
            % Update time
            t = t+dt;
    
            fprintf(['tIndex= %.0f, t(hr)=%.1f, DRO Error=%.3f, ELFO Error=%.3f, ' ...
                'DRO NEES(orb)=%.2f, ELFO NEES(orb)=%.2f, DRO NEES(clc)=%.2f, ELFO NEES(clc)=%.2f, ELFO_NIS=%.2f\n'], ...
                tIndexMeas, (t-t0)/3600, norm(xEst_DRO(1:3)-xTrue_DRO(1:3)), norm(xEst_ELFO(1:3)-xTrue_ELFO(1:3)), ...
                NEES_DRO_acc, NEES_ELFO_acc, NEES_DRO_clc, NEES_ELFO_clc, NIS_ELFO);
    
    end

    MarrayT(:,:,n)      = arrayT;
    MarrayX_DRO(:,:,n)  = arrayX_DRO;
    MarrayX_ELFO(:,:,n) = arrayX_ELFO;

    MarrayNEES_DRO_acc(:,:,n)  = arrayNEES_DRO_acc;
    MarrayNEES_ELFO_acc(:,:,n) = arrayNEES_ELFO_acc;
    MarrayNEES_DRO_clc(:,:,n)  = arrayNEES_DRO_clc;
    MarrayNEES_ELFO_clc(:,:,n) = arrayNEES_ELFO_clc;

    MarrayNIS_DRO(:,:,n)  = arrayNIS_DRO;
    MarrayNIS_ELFO(:,:,n) = arrayNIS_ELFO;

end

%% Save and Output
for i = 1:N

    diff_DRO = orbTrue1.seq.a.XJ2000 - MarrayX_DRO(:,:,i);

    posErrors_DRO(:,:,i) = vecnorm(diff_DRO(:,1:3), 2, 2);
    velErrors_DRO(:,:,i) = vecnorm(diff_DRO(:,4:6), 2, 2);
    clcErrors_DRO(:,:,i) = diff_DRO(:,7:8);

    diff_ELFO = orbTrue2.seq.a.XJ2000 - MarrayX_ELFO(:,:,i);

    posErrors_ELFO(:,:,i) = vecnorm(diff_ELFO(:,1:3), 2, 2);
    velErrors_ELFO(:,:,i) = vecnorm(diff_ELFO(:,4:6), 2, 2);
    clcErrors_ELFO(:,:,i) = diff_ELFO(:,7:8);

end

posRMSE_DRO = sqrt(mean(posErrors_DRO.^2, 3));
velRMSE_DRO = sqrt(mean(velErrors_DRO.^2, 3));
clcRMSE_DRO = sqrt(mean(clcErrors_DRO.^2, 3));

posRMSE_ELFO = sqrt(mean(posErrors_ELFO.^2, 3));
velRMSE_ELFO = sqrt(mean(velErrors_ELFO.^2, 3));
clcRMSE_ELFO = sqrt(mean(clcErrors_ELFO.^2, 3));

arrayNEES_DRO_acc = mean(MarrayNEES_DRO_acc, 3);
arrayNEES_ELFO_acc = mean(MarrayNEES_ELFO_acc, 3);
arrayNEES_DRO_clc = mean(MarrayNEES_DRO_clc, 3);
arrayNEES_ELFO_clc = mean(MarrayNEES_ELFO_clc, 3);
arrayNIS_DRO = mean(MarrayNIS_DRO, 3);
arrayNIS_ELFO = mean(MarrayNIS_ELFO, 3);

orbEs1.seq.a.t         = MarrayT(:,:,1);
orbEs1.seq.a.posErrors = posRMSE_DRO;
orbEs1.seq.a.velErrors = velRMSE_DRO;
orbEs1.seq.a.clcErrors = clcRMSE_DRO;

orbEst2.seq.a.t         = MarrayT(:,:,1);
orbEst2.seq.a.posErrors = posRMSE_ELFO;
orbEst2.seq.a.velErrors = velRMSE_ELFO;
orbEst2.seq.a.clcErrors = clcRMSE_ELFO;

orb.EstDRO   = orbEs1;
orb.TrueDRO  = orbTrue1;
orb.EstELFO  = orbEst2;
orb.TrueELFO = orbTrue2;

save('output/ORBdataSST','orb');

%% Close Path
rmpath('prop');
rmpath('input');
rmpath("mice/lib");
rmpath("mice/src/mice");

toc;

function [xEst, PEst, NIS] = UKF( ...
    xEst, PEst, ...   
    P, R, ...           
    time, dt, orb, ...        
    yMeas, measType, xSource, PSource, ... 
    para, options)

    alpha = para(1); beta = para(2); kappa = para(3);
    
    % Dimensions
    Lx = length(xEst);                          % Augmented state dimension (=8)
    Lv = size(P, 1);                            % Process noise dimension (=8)
    Ln = size(R, 1);                            % Measurement dimension (=4)
    La = Lx + Lv + Ln;                          % Total dimension (=20)
    
    % Augmented state & augmented covariance
    xAug = [xEst; zeros(Lv,1); zeros(Ln,1)];    % [x; 0; 0]
    P_a  = blkdiag(PEst, P, R);                 % [Px 0  0
                                                %  0  Q  0
                                                %  0  0  R]
    % Weights
    lambda = alpha^2*(La + kappa) - La;
    gamma  = sqrt(La + lambda);
    
    % Generate augmented sigma points
    sqrtPa  = chol_spd(P_a, 1e-12);      % La × La
    xSigAug = [ xAug , ...
                xAug + gamma*sqrtPa , ...
                xAug - gamma*sqrtPa ];   % La × (2·La+1)
    
    % Disassemble sigma points
    xSig = xSigAug(1:Lx, :);                  % State
    vSig = xSigAug(Lx+1:Lx+Lv, :);            % Process noise
    nSig = xSigAug(Lx+Lv+1:Lx+Lv+Ln, :);      % Measurement noise
    
    % Calculate predictions for sigma points
    xSigPred = zeros(Lx, 2*La+1);
    for i = 1:(2*La+1)
        % Pure Propagation using the low-fidelity dynamic model
        [~, Xtemp] = ode113(@(t,X) prophpopClock(t, X, orb), [time (time+dt)], xSig(:, i), options);
        % Add process noise
        xSigPred(:, i) = Xtemp(end,:)' + vSig(:,i);
    end
    
    % UKF weights
    N  = 2*La + 1;                    
    c0 = lambda / (La + lambda);      
    c  = 0.5   / (La + lambda);      
    
    Wm          = [c0 ;  c*ones(N-1,1)];   % N×1
    Wc          =  Wm;                   
    Wc(1)       =  Wc(1) + (1 - alpha^2 + beta);

    % Mean of prediction
    xPred = xSigPred * Wm;
    
    % Covariance of prediction
    diffX = xSigPred - xPred;        % n × N  
    PPred = (diffX .* Wc.') * diffX.';   % n × n 
    
    % If there is no measurement, return the prediction
    if isempty(yMeas)
        xEst = xPred;
        PEst = PPred;
        PEst = (PEst + PEst.')/2;
        PEst = nearestSPD(PEst);
        PEst = (PEst + PEst.')/2 + 1e-12*eye(size(PEst));  
        NIS = NaN;
        return;
    end
    
    % Measurement update
    switch measType

        case 'Earth' % Ground measurement

            % 7.1 Calculate measurement results for sigma points
            ySig = zeros(Ln, 2*La+1);
            for i = 1:(2*La+1)
                % Pure predicted measurement
                y_det_i = MeasurementModel(xSigPred(:,i), time+dt, measType, xSource);
                % Add measurement noise
                ySig(:,i) = y_det_i + nSig(:,i);
            end
            
            % 7.2 Mean of predicted measurement
            yPred = ySig * Wm;
        
            % Angle wrapping for residuals
            angIdx = [2, 3]; % [elev, azim]
            
            diffX = xSigPred - xPred;      % n × N
            diffZ = ySig       - yPred;    % m × N
            
            for j = angIdx
                diffZ(j,:) = atan2(sin(diffZ(j,:)), cos(diffZ(j,:)));
            end
            
            innov = yMeas - yPred;
            for j = angIdx
                innov(j) = atan2(sin(innov(j)), cos(innov(j)));
            end
            
            Pxy = (diffX .* Wc.') * diffZ.';   % n × m
            Pyy = (diffZ .* Wc.') * diffZ.';   % m × m
            Pyy = (Pyy + Pyy.')/2 + 1e-12*eye(size(Pyy));
            
            K   = Pxy / Pyy;
            xEst = xPred + K*innov;
            PEst = PPred - K*Pyy*K';
            PEst = (PEst + PEst.')/2;
            PEst = nearestSPD(PEst);
            PEst = (PEst + PEst.')/2 + 1e-12*eye(size(PEst)); 
            NIS  = innov' * (Pyy \ innov);


        case 'SST2' % Relay SST measurement

            % 7.1 Calculate measurement results for sigma points
            ySig = zeros(Ln, 2*La+1);
            for i = 1:(2*La+1)
                % Pure predicted measurement
                y_det_i = MeasurementModel(xSigPred(:,i), time+dt, measType, xSource);
                % Add measurement noise
                ySig(:,i) = y_det_i + nSig(:,i);
            end
    
            % 7.2 Mean of predicted measurement
            yPred = ySig * Wm;
        
            % 7.3 Measurement covariance Pzz and cross-covariance Pxz
            diffX = xSigPred - xPred;          % n × N
            diffZ = ySig       - yPred;        % m × N

            Pxy = (diffX .* Wc.') * diffZ.';   % n × m
            Pyy = (diffZ .* Wc.') * diffZ.';   % m × m

            H    = BuildJacobianMatrix(xPred, xSource);
            Reff = H * PSource * H.';

            Pyy  = Pyy + Reff;
            Pyy  = (Pyy + Pyy.')/2 + 1e-15*eye(size(Pyy));

            % Measurement update and return
            K    = Pxy / Pyy;
            xEst = xPred + K*(yMeas - yPred);
            PEst = PPred - K*Pyy*K';
            PEst = (PEst + PEst.')/2;
            PEst = nearestSPD(PEst);
            PEst = (PEst + PEst.')/2 + 1e-12*eye(size(PEst));
            NIS  = (yMeas - yPred)' * (Pyy \ (yMeas - yPred));

    end
end

function [xEst1, PEst1, xEst2, PEst2, P12_out, NIS] = JUKF( ...
    xEst1, PEst1, xEst2, PEst2, P12_in, ...          
    P1, P2, R, time, dt, orb1, orb2, ...             
    yMeas, measType, ...                             
    para1, para2, paraEst, options)                 

    % UKF Parameters (1 is ELFO, 2 is DRO)
    [alpha1,beta1,kappa1] = deal(para1(1), para1(2), para1(3));
    [alpha2,beta2,kappa2] = deal(para2(1), para2(2), para2(3));
    [alphaE,betaE,kappaE] = deal(paraEst(1), paraEst(2), paraEst(3));
    
    % Dimensions
    Lx1 = length(xEst1);   Lx2 = length(xEst2);
    Lv1 = size(P1,1);      Lv2 = size(P2,1);
    
    %% ================= ELFO Prediction =======================
    La1 = Lx1 + Lv1;
    [gamma1,Wm1,Wc1] = CalculateWeights(alpha1, beta1, kappa1, La1);
    
    xAug1 = [xEst1; zeros(Lv1,1)];
    Paug1 = blkdiag(PEst1,P1);
    S1    = chol_spd(Paug1, 1e-12);
    X1    = [xAug1, xAug1 + gamma1*S1, xAug1 - gamma1*S1];
    
    X1p = zeros(Lx1,2*La1+1);
    for i = 1:2*La1+1
        [~,Xt] = ode113(@(t,x) prophpopClock(t,x,orb1), ...
                        [time,time+dt], X1(1:Lx1,i), options);
        X1p(:,i) = Xt(end,:)' + X1(Lx1+1:end,i);
    end
    xPred1 = X1p*Wm1;
    dX1    = X1p - xPred1;
    PPred1 = (dX1.*Wc1')*dX1';
    
    %% ================= DRO Prediction =======================
    La2 = Lx2 + Lv2;
    [gamma2,Wm2,Wc2] = CalculateWeights(alpha2, beta2, kappa2, La2);
    
    xAug2 = [xEst2; zeros(Lv2,1)];
    Paug2 = blkdiag(PEst2,P2);
    S2 = chol_spd(Paug2, 1e-12);
    X2    = [xAug2, xAug2 + gamma2*S2, xAug2 - gamma2*S2];
    
    X2p = zeros(Lx2,2*La2+1);
    for i = 1:2*La2+1
        [~,Xt] = ode113(@(t,x) prophpopClock(t,x,orb2), ...
                        [time,time+dt], X2(1:Lx2,i), options);
        X2p(:,i) = Xt(end,:)' + X2(Lx2+1:end,i);
    end
    xPred2 = X2p*Wm2;
    dX2    = X2p - xPred2;
    PPred2 = (dX2.*Wc2')*dX2';
    
    %% =============== %%% =====================
    if isempty(P12_in), P12_prev = zeros(Lx1,Lx2); else, P12_prev = P12_in; end
    
    intStep = 10;                            
    Nstep = max(1, round(dt/intStep));
    
    Phi1 = eye(Lx1);   xInt1 = xEst1;
    for k = 1:Nstep
        A1 = numJac_dim(@(x) prophpopClock(time+(k-1)*intStep, x, orb1), xInt1);
        Phi1 = (eye(Lx1)+A1*intStep)*Phi1;
        [~,xt] = ode113(@(t,x) prophpopClock(t,x,orb1),[0 intStep],xInt1,options);
        xInt1 = xt(end,:)';
    end
    
    Phi2 = eye(Lx2);   xInt2 = xEst2;
    for k = 1:Nstep
        A2 = numJac_dim(@(x) prophpopClock(time+(k-1)*intStep, x, orb2), xInt2);
        Phi2 = (eye(Lx2)+A2*intStep)*Phi2;
        [~,xt] = ode113(@(t,x) prophpopClock(t,x,orb2),[0 intStep],xInt2,options);
        xInt2 = xt(end,:)';
    end
    
    P12_pred = Phi1 * P12_prev * Phi2';
    
    % Joint Predictions
    xPred = [xPred1; xPred2];
    PPred = [PPred1, P12_pred; P12_pred', PPred2];
    PPred = (PPred + PPred.')/2;

    S1 = chol_spd(PPred1, 1e-12);
    S2 = chol_spd(PPred2, 1e-12);

    M  = S1 \ (P12_pred / S2');
    smax = svds(M, 1);
    limit = 0.98;
    if smax >= limit
        P12_pred = P12_pred * (limit / smax * 0.999);
        PPred = [PPred1, P12_pred; P12_pred', PPred2];
        PPred = (PPred + PPred.')/2;
    end
    
    PPred = nearestSPD(PPred);
    
    %% ================= Joint Measurement Update ==================
    m  = size(R,1);
    Ls = Lx1 + Lx2;          
    La = Ls + m;              
    [gammaE, WmE, WcE] = CalculateWeights(alphaE, betaE, kappaE, La);
    
    xAug = [xPred; zeros(m,1)];
    Paug = blkdiag(PPred, R);   
    
    Saug = chol_spd(Paug, 1e-12);
    X = [xAug, xAug + gammaE*Saug, xAug - gammaE*Saug];
    
    Xstate = X(1:Ls, :);
    Nmeas  = X(Ls+1:end, :);
    
    ySig = zeros(m, 2*La+1);
    for i = 1:2*La+1
        x1_i = Xstate(1:Lx1, i);
        x2_i = Xstate(Lx1+1:end, i);
        ySig(:, i) = MeasurementModel(x1_i, time+dt, measType, x2_i) + Nmeas(:, i);
    end
    
    yPred = ySig * WmE;
    dX    = Xstate - xPred;
    dZ    = ySig   - yPred;
    
    Pxy = (dX .* WcE') * dZ.'; 
    Pyy = (dZ .* WcE') * dZ.';   
    Pyy = Pyy + 1e-12*eye(size(Pyy));
    
    opts.SYM = true; opts.POSDEF = true;
    K    = linsolve(Pyy.', Pxy.', opts).';
    innov = yMeas - yPred;
    
    xEst = xPred + K*innov;
    PEst = PPred - K*Pyy*K';
    PEst = (PEst + PEst.')/2;
    PEst = nearestSPD(PEst);
    PEst = (PEst + PEst.')/2 + 1e-12*eye(size(PEst));  
    
    xEst1  = xEst(1:Lx1);
    xEst2  = xEst(Lx1+1:end);
    PEst1  = PEst(1:Lx1, 1:Lx1);
    PEst2  = PEst(Lx1+1:end, Lx1+1:end);
    P12_out= PEst(1:Lx1, Lx1+1:end);
    
    NIS = innov' * (linsolve(Pyy, innov, opts));


end

function A = numJac(f, x0)
    n = length(x0);
    A = zeros(n);
    for j = 1:n
        h = max(1e-6, 1e-6*abs(x0(j)));       % Adaptive step size
        dx = zeros(n,1); dx(j) = h;
        A(:,j) = (f(x0+dx) - f(x0-dx)) / (2*h);
    end
end

function J = numJac_dim(f, x0, opts)

    x0 = x0(:);
    n  = numel(x0);
    y0 = f(x0);
    m  = numel(y0);
    J  = zeros(m, n);

    if ~exist('opts','var') || isempty(opts), opts = struct(); end
    if ~isfield(opts,'abs_step') || isempty(opts.abs_step)
        abs_step = [1e-6*ones(3,1);  ... % pos [km]
                    1e-9*ones(3,1);  ... % vel [km/s]
                    1e-11;          ...  % clock bias [s]
                    1e-14];               % clock drift [s/s]

        if numel(abs_step) ~= n
            abs_step = 1e-8*ones(n,1);
        end
    else
        abs_step = opts.abs_step(:);
    end

    if ~isfield(opts,'rel_step') || isempty(opts.rel_step)
        rel_step = 1e-6 * ones(n,1);     
    else
        rel_step = opts.rel_step(:);
    end
    if ~isfield(opts,'h_min') || isempty(opts.h_min)
        h_min = abs_step;               
    else
        h_min = opts.h_min(:);
    end

    h = max(abs_step, rel_step .* max(abs(x0), 1)); 
    h = max(h, h_min);                              

    use_complex = isfield(opts,'use_complex') && opts.use_complex;

    if use_complex

        for j = 1:n
            xh = x0;
            xh(j) = xh(j) + 1i*h(j);
            y = f(xh);
            J(:,j) = imag(y) / h(j);
        end
    else

        for j = 1:n
            xh_plus  = x0;  xh_plus(j)  = xh_plus(j)  + h(j);
            xh_minus = x0;  xh_minus(j) = xh_minus(j) - h(j);
            y_plus  = f(xh_plus);
            y_minus = f(xh_minus);
            J(:,j) = (y_plus - y_minus) / (2*h(j));
        end
    end
end

function z = MeasurementModel(xTarget, et, measType, xSource) 

    c = 299792.458; % km/s

    switch measType

        % Ground measurement
        case 'Earth'
            % Calculate the measurement [range; elev; azim; rangeRate]
            % by spacecraft's state [rx ry rz vx vy vz]

            % Obtain the coordinates of the moon relative to the Earth
            [x_Moon2Earth, ~] = cspice_spkezr('MOON', et, 'J2000', 'LT+S', 'EARTH');
            % The coordinates of the satellite in the Earth J2000 system
            x_sat2Earth = xTarget(1:6) + x_Moon2Earth;
        
            rx = x_sat2Earth(1); ry = x_sat2Earth(2); rz = x_sat2Earth(3);
            vx = x_sat2Earth(4); vy = x_sat2Earth(5); vz = x_sat2Earth(6);
            r_vec = [rx; ry; rz];
            v_vec = [vx; vy; vz];
        
            r = norm(r_vec) + c * xTarget(7); % Radial range
            elev = atan2(rz, sqrt(rx^2 + ry^2));
            elev = wrapPi(elev);
            azim = atan2(ry, rx);
            azim = wrapPi(azim);
            rdot = dot(r_vec, v_vec)/r + c * xTarget(8);  % Radial rangerate
        
            z = [r; elev; azim; rdot];

        % Range only SST
        case 'SST1'
            % Measurement elements (radial range only)
            z = norm(xTarget(1:3) - xSource(1:3)) + c * (xTarget(7) - xSource(7));
        
        % Range and range-rate SST
        case 'SST2'
            relativePos = xTarget(1:3) - xSource(1:3);
            relativeVel = xTarget(4:6) - xSource(4:6);
            rho = norm(relativePos);
            rho = max(rho, 1e-15);  
            vLOS = dot(relativePos, relativeVel) / rho;
            z = [rho; vLOS] + [c*(xTarget(7)-xSource(7)); c*(xTarget(8)-xSource(8))];

        otherwise
            z = [];
    end

end

function isVisible = IsVisible(target, et, sourceType, xSource)

    % Both the input of target and origin (if it is a satellite) are in Moon
    % centered J2000 reference system

    % Get the Moon’s position in Earth J2000
    [xMoonFromEarth, ~] = cspice_spkezr('MOON', et, 'J2000', 'LT+S', 'EARTH');

    % Decide how to define the observer
    if strcmpi(sourceType, 'Earth') % Origin is Earth
        
        x_MoonFromSource = xMoonFromEarth(1:3);
        x_TargetFromSource  = target(1:3) + xMoonFromEarth(1:3);

    elseif strcmpi(sourceType, 'Satellite') % Origin is another satellite
        
        x_MoonFromSource = -xSource(1:3); 
        x_TargetFromSource  = target(1:3) - xSource(1:3);

    else
        error('originType must be either "Earth" or "Satellite".');
    end

    % Compute distances
    rMoon = norm(x_MoonFromSource);
    rSat  = norm(x_TargetFromSource);

    % "Project" the target onto the sphere at the Moon's distance 
    if rSat > 0
        x_targetProjected = (rMoon / rSat) .* x_TargetFromSource;
    else
        % Degenerate case: if rSat=0, target and origin are the same point!
        % Let's just set them so no division by zero occurs:
        x_targetProjected = zeros(3,1);
    end

    % Apparent offset from Moon center 
    rRelative = norm(x_targetProjected - x_MoonFromSource);

    if rMoon > 0
        theta = asin(rRelative / rMoon);        % Target's apparent angle from Moon center
        beta  = asin(1737.4 / rMoon);           % Moon's angular radius
    else
        % Another degenerate case: if rMoon=0, the observer is at the Moon's center
        theta = 0;
        beta  = pi;  
                     
    end

    % Angle between satellite and Moon center, as seen by the origin
    cos_alpha = dot(x_TargetFromSource, x_MoonFromSource) / (rSat * rMoon + eps);
    cos_alpha = max(-1, min(1, cos_alpha)); 
    alpha     = acos(cos_alpha);

    % Logical judgements
    outsideLimb  = (theta > beta);
    inFrontMoon  = (rSat < rMoon) && (alpha < pi/2);

    isVisible = outsideLimb || inFrontMoon;
end

function H = BuildJacobianMatrix(xT,xS)
    
    c = 299792.458;

    dr = xT(1:3) - xS(1:3);
    dv = xT(4:6) - xS(4:6);
    rho = norm(dr);
    rho = max(rho, 1e-15);   
    rho_d = dot(dr, dv)/rho;
    u_r = (dr.'/rho);
    
    H = [ -u_r,                      zeros(1,3),  c, 0; 
          -(dv.'/rho)+u_r*rho_d/rho, -u_r,        0, c ];
end

function [gamma, Wm, Wc] = CalculateWeights(alpha, beta, kappa, La)

    lambda = alpha^2*(La + kappa) - La;
    gamma  = sqrt(La + lambda);

    N  = 2*La + 1;                     % Total sigma number
    c0 = lambda / (La + lambda);       % First weight
    c  = 0.5   / (La + lambda);        % Other weight

    Wm    = [c0 ;  c*ones(N-1,1)];   % N×1
    Wc    =  Wm;                     
    Wc(1) =  Wc(1) + (1 - alpha^2 + beta);
end
      
function a = wrapPi(a)
    a = atan2(sin(a), cos(a));
end

function Ahat = nearestSPD(A)
    A  = (A + A')/2;
    [V,D] = eig(A);
    D     = diag(D);
    D(D < 0) = 0;
    Ahat  = V*diag(D)*V';
    Ahat  = (Ahat + Ahat')/2;
    Ahat  = Ahat + 1e-15*eye(size(Ahat));
end

function S = chol_spd(A, base_eps)
    if nargin < 2, base_eps = 1e-12; end
    A = (A + A.')/2;           
    [S,p] = chol(A,'lower');
    if p==0, return; end            

    eps_now = base_eps;
    [S,p] = chol(A + eps_now*eye(size(A)),'lower');
    if p==0, return; end

    Ahat = nearestSPD(A);
    Ahat = (Ahat + Ahat.')/2;
    eps_now = max(1e-12, eps_now);
    for k = 1:3         
        [S,p] = chol(Ahat + eps_now*eye(size(Ahat)),'lower');
        if p==0, return; end
        eps_now = eps_now*10;
    end
    error('chol_spd:Failed','Failed to make SPD');
end

