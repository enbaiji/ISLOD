function [orb] = LoadState(Orbit,orb)
%LoadState This function gathers the states of some selected orbits of
%interest. 
%       The Orbits called "RefSpacecraft" are the orbit of interest used in the paper.
%       The NRHO and DRO14 are CR3BP orbits (close to the Capstone and Orion one) fitted to ephemeris model.
%       ELFO is an Elliptical Lunar Frozen Orbit.
%       Two examples/templates for Keplerian coordinates, Kepl and Kepl2, were also added.
%       The user can add new ones or modify the "Perso" Orbit in last position.

    if Orbit == "NRHO"
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)
        
        X = -14007.49388066819; % X-Coordinate
        Y = 37926.40572749509;    % Y-Coordinate
        Z = -59472.47612911355;    % Z-Coordinate
        VX = 0.03815140145811324;    % X-Velocity Coordinate
        VY = 0.05511836319581054;    % Y-Velocity Coordinate
        VZ  = 0.028138675143097314;    % Z-Velocity Coordinate

        Time = '2024 Jan 01 00:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]

    elseif Orbit == "DRO14"
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)
        
        X = 63650.58084169839; % X-Coordinate
        Y = -24686.81610635431;    % Y-Coordinate
        Z = -15448.941380173477;    % Z-Coordinate
        VX = -0.14532575232216205;    % X-Velocity Coordinate
        VY = -0.290549113768819;    % Y-Velocity Coordinate
        VZ  = -0.1502956531877957;    % Z-Velocity Coordinate

        Time = '2024 Jan 01 00:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]

    elseif Orbit == "ELFO"
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)

        X = 12.931; % X-Coordinate
        Y = 235.558;    % Y-Coordinate
        Z = 2388.38;    % Z-Coordinate
        VX = 1.80519;    % X-Velocity Coordinate
        VY = -0.0991;    % Y-Velocity Coordinate
        VZ = -3.83e-07;    % Z-Velocity Coordinate

        Time = '2024 Jan 01 00:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]

    elseif Orbit == "ELFO2"
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)

        X = 396.8442912666942; % X-Coordinate
        Y = -2468.0076512364612;    % Y-Coordinate
        Z = 1132.7859839924058;    % Z-Coordinate
        VX = 1.704357257513478;    % X-Velocity Coordinate
        VY = -0.05316874438678221;    % Y-Velocity Coordinate
        VZ  = -0.08851852089230178;    % Z-Velocity Coordinate

        Time = '2024 Jan 01 00:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]

    elseif Orbit == "L1Halo"
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)

        X = 58641.62585597829; % X-Coordinate
        Y = -28334.082269736246;    % Y-Coordinate
        Z = -3417.061795447473;    % Z-Coordinate
        VX = 0.010324360067413718;    % X-Velocity Coordinate
        VY = 0.0055308532041151715;    % Y-Velocity Coordinate
        VZ  = 0.00259459252901334;    % Z-Velocity Coordinate

        Time = '2024 Jan 01 00:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]

    elseif Orbit == "L2HaloS"
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)

        X = -64083.40215841029; % X-Coordinate
        Y = 44807.07340175446;    % Y-Coordinate
        Z = -23052.134640715994;    % Z-Coordinate
        VX = 0.025796655339742514;    % X-Velocity Coordinate
        VY = 0.029131152894391166;    % Y-Velocity Coordinate
        VZ  = 0.014673293613276073;    % Z-Velocity Coordinate

        Time = '2024 Jan 01 00:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]

    elseif Orbit == "RefClementine"
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)

        X = -2098.4745703999997; % X-Coordinate
        Y = 1298.5243538999998;    % Y-Coordinate
        Z = -3494.2586871999997;    % Z-Coordinate
        VX = 0.8948704814799997;    % X-Velocity Coordinate
        VY = 0.17693000307000004;    % Y-Velocity Coordinate
        VZ = -0.15601097269;    % Z-Velocity Coordinate

        Time = '1994 Apr 15 15:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]

    elseif Orbit == "RefCapstone"
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)

        X = -16983.14075642353; % X-Coordinate
        Y = 21213.5584242304;    % Y-Coordinate
        Z = -58035.6304537942;    % Z-Coordinate
        VX = -0.04457645856905286;    % X-Velocity Coordinate
        VY = -0.05137482728591616;    % Y-Velocity Coordinate
        VZ = 0.1416421540429468;    % Z-Velocity Coordinate

        Time = '2022 Nov 25 00:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]
    
    elseif Orbit == "RefOrion"
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)

        X = 41570.98496247831; % X-Coordinate
        Y = -50015.47166299814;    % Y-Coordinate
        Z = -28760.46002667214;    % Z-Coordinate
        VX = -0.2465381485425481;    % X-Velocity Coordinate
        VY = -0.182296424457602;    % Y-Velocity Coordinate
        VZ = -0.07642694222326024;    % Z-Velocity Coordinate

        Time = '2022 Nov 29 16:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]

    elseif Orbit == "RefLRO"
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)

        X = 820.1514640535461; % X-Coordinate
        Y = 1091.2367697664858;    % Y-Coordinate
        Z = 1262.7894519284164;    % Z-Coordinate
        VX = 0.22458861495091112;    % X-Velocity Coordinate
        VY = 1.1177228435015545;    % Y-Velocity Coordinate
        VZ = -1.1391038535297986;    % Z-Velocity Coordinate
        
        Time = '2020 Feb 01 00:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]
    
    elseif Orbit == "Kepl"
        RefS = 'MOON_ME';   % Reference System (J2000/MOON_ME)
        CoordT = 'Keplerian'; %Coordinates Type (Cartesian/Keplerian)

        SMA = 1738+100;    % semi-major axis             [km]
        ECC = 0;           % eccentricity
        INC = 90;          % inclination                 [°]
        LAN = 0;           % longitude of ascending node [°]
        AOP = 45;          % argument of pericenter      [°]
        MA = 0;            % true anomaly                [°]
        
        Time = '2020 Feb 01 00:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]
    
    elseif Orbit == "Kepl2"
        RefS = 'MOON_ME';   % Reference System (J2000/MOON_ME)
        CoordT = 'Keplerian2'; %Coordinates Type: Similar as Keplerian but 
                               % with a different definition of the ellipse geometry

        PER = 1738+100;    % perilune                    [km]
        APO = 1738+1000;   % apolune                     [km]
        INC = 90;          % inclination                 [°]
        LAN = 0;           % longitude of ascending node [°]
        AOP = 45;          % argument of pericenter      [°]
        MA = 0;            % true anomaly                [°]
        
        Time = '2020 Feb 01 00:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]
    elseif Orbit == "Perso"
        RefS = 'J2000';   % Reference System (J2000/MOON_ME)
        CoordT = 'Cartesian'; %Coordinates Type (Cartesian/Keplerian)

        X = -175.001861; % X-Coordinate
        Y = -1175.17121;    % Y-Coordinate
        Z = 1410.69401;    % Z-Coordinate
        VX = 0.791300523;    % X-Velocity Coordinate
        VY = 1.03420255;    % Y-Velocity Coordinate
        VZ = 0.971046194;    % Z-Velocity Coordinate

        Time = '2022 Nov 26 00:00:00.000';  % start time [yyyy month dd hh:mm:ss.---]
    end

%% Initialization
% position and velocity
    % =========================================================================% 
    orb.frame.initstate = RefS;
    et = cspice_str2et(Time);
    if strcmp(CoordT,'Cartesian')
        orb.sat.X0 = [X,Y,Z,VX,VY,VZ];
        if ~strcmp(orb.frame.initstate,{'ICRF','J2000'})
            Rgeog_iner = cspice_sxform(orb.frame.initstate,orb.frame.to,et);
            orb.sat.X0iner = Rgeog_iner*orb.sat.X0';
        else
            orb.sat.X0iner = orb.sat.X0;
        end
    else
        if strcmp(CoordT,'Keplerian')
            orb.sat.keplstate(1,1:2) = {'SMA',SMA};            % semi-major axis         [km]
            orb.sat.keplstate(2,1:2) = {'ECC',ECC};            % eccentricity
        else
            orb.sat.keplstate(1,1:2) = {'SMA',(APO+PER)/2};            % perilune                [km]
            orb.sat.keplstate(2,1:2) = {'ECC',1-PER/orb.sat.keplstate{1,2}};            % apolune                 [km]
        end
        orb.sat.keplstate(3,1:2) = {'INC',INC*(pi/180)};   % inclination                 [rad]
        orb.sat.keplstate(4,1:2) = {'LAN',LAN*(pi/180)};   % longitude of ascending node [rad]
        orb.sat.keplstate(5,1:2) = {'AOP',AOP*(pi/180)};   % argument of pericenter      [rad]
        orb.sat.keplstate(6,1:2) = {'MA',MA*(pi/180)};     % true anomaly                [rad]
        orb.sat.X0 = (cspice_conics([ orb.sat.keplstate{1,2}*(1-orb.sat.keplstate{2,2}),orb.sat.keplstate{2,2},orb.sat.keplstate{3,2},...
                                    orb.sat.keplstate{4,2},orb.sat.keplstate{5,2},orb.sat.keplstate{6,2},...
                                    et,orb.centralPlanet.GM]',et))';
        if ~strcmp(orb.frame.initstate,{'ICRF','J2000'})
            Rgeog_iner = cspice_pxform(orb.frame.initstate,orb.frame.to,et);
            orb.sat.X0iner(1:3,1) = Rgeog_iner*orb.sat.X0(1,1:3)';
            orb.sat.X0iner(4:6,1) = Rgeog_iner*orb.sat.X0(1,4:6)';
        else
            orb.sat.X0iner = orb.sat.X0;
        end
    end
    orb.sat.t0 = et;
    % =========================================================================% 

end

