% HALO, High-precision Analyser of Lunar Orbits, is here built for the propagation of a lunar orbiter.
% The motion is integrated under 
% - the Keplerian attraction force
% - the perturbations due to the asymmetry of the lunar gravity field (up to 350*350 harmonics)
% - the perturbations due to the asymmetry of the Earth gravity field (up to 100*100 harmonics)
% - the third body attractions of the Sun and Jupiter
% - the solar radiation pressure
% - the General relativity
% - the Earth's albedo
%     
% The integration is performed in the "J2000" inertial frame, and the output is the state and time vectors.
    
clear;
ts = cputime;

%% SPICE LIBRARIES
addpath("mice/lib");
addpath("mice/src/mice");
addpath('prop');
addpath('input');

cspice_kclear; 
metakernelcheck;
cspice_furnsh('metakernel.tm');

orbitList = ["halo", "elfo2"];
% options = odeset('RelTol',1e-7,'AbsTol',1e-14); % For low noise and high precision
options = odeset('RelTol',1e-6,'AbsTol',1e-9); 
tStep = 60;
noiseStep = 60;
rng(0);

for n = 1:length(orbitList)

    orbit = orbitList(n);
    
    %% Initialization
    orb = struct();
    orb.model = "Truth";
    orb.type = orbit;

    orb = LoadSequential(orb);
    Time = orb.seq.Time;
    SeqNames = fieldnames(orb.seq);

    orb.noise.clcNoise = 1e-28; % "q" (previously 1e-14)
    % orb.noise.clcNoise = 1e-99; % "q" (previously 1e-14)

    b0  = 0; % Initial clock bias
    db0 = 0; % Initial clock drift
    
    x6   = orb.sat.X0iner;
    x8 = [x6, b0, db0];
    orb.sat.X0iner = x8;
    orb.sat.b0 = [b0, db0];

    %%  Running the sequence
    for i = 2:length(SeqNames)
        Seq = orb.seq.(SeqNames{i});
        
        if Seq.type == "Propag"
            span = Seq.span;
            Tspan = Time:tStep:Time+span;

            N = floor(span / noiseStep);
            NoiseClock = zeros(N,1);
        
            % Construct the table of clock noise to be drawn from
            for j = 1:N
                NoiseClock(j,:) = sqrt(orb.noise.clcNoise / noiseStep) * randn(1,1);
            end
        
            orb.noise.dt  = noiseStep;
            orb.noise.clc = NoiseClock;

            [orb.seq.(SeqNames{i}).t, orb.seq.(SeqNames{i}).XJ2000]...
                = ode113(@prophpopNoisedClock, Tspan, orb.sat.X0iner, options, orb);
        else
            error('Unrecognizable Sequence Type');
        end

        orb.sat.X0iner = orb.seq.(SeqNames{i}).XJ2000(end,:);
        Time = Time+span;
        orb.seq.Time = Time;
    end

    save(['output/ORBdataTrue',num2str(n)],'orb');

end

disp("Run Time: "+num2str(cputime - ts)+"s")

%% Close Path
rmpath('prop');
rmpath('input');
rmpath("mice/lib");
rmpath("mice/src/mice");


