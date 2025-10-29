%This code has been deeply optimized to maximize computation speed where
%possible

tic
syms r f X Y Z XF YF ZF xh yh zh
syms lambda_s phi_s h_s

phi_values = linspace(pi/2, -pi/2, 10);
lambda_values = linspace(-pi, pi, 20);
[Lambda, Phi] = meshgrid(lambda_values, phi_values);

a=6541400; % Semi-major axis
e=0.6; % Eccentricity
inc=56.2*pi/180; % Inclination
omega=90*pi/180; % Argument of periapsis
t0=0;
miu = 6.6743*10^(-11)*7.3477*10^22; % Standard gravitational parameter
Ele = 5*pi/180; % Elevation to analysis coverage

t_values = linspace(0, 3600*24*10, 1000);
Omega_values = [0, pi]; % RAAN
t_DOP = 1;

% Preallocate spaces
GDOP = zeros(length(phi_values), length(lambda_values), length(t_values));
PDOP = zeros(length(phi_values), length(lambda_values), length(t_values));
HDOP = zeros(length(phi_values), length(lambda_values), length(t_values));
AA = zeros(4, 4, numel(Phi));

% Outer loop: iterate over t values
parfor i = 1:length(t_values)
    t = t_values(i);

    % % Initialize cell arrays to store H and Sg matrices for different Omega values
    H_matrices = cell(length(Omega_values), 1);
    sat_in_view_matrices = cell(length(Omega_values), 1);

    M0_values = zeros(1, 4); % Mean anomaly

    sat_in_view_matrix = zeros(length(M0_values), length(M0_values));

    % Middle loop: iterate over Omega values
    for o = 1:length(Omega_values)
        Omega = Omega_values(o);

        % Initialize matrices to store H values for different M0 values,
        % along with different coordinates 
        H_matrix = zeros(length(M0_values), 4, numel(Phi));

        % Reassign different values to M0 based on Omega
        if Omega == 0
            M0_values = [0, pi/2, pi, pi/2*3];
        elseif Omega == pi
            M0_values = [pi/4, pi/4*3, pi/4*5, pi/4*7];
            % M0_values = [0, pi/2, pi, pi/2*3];
        else
            error('Invalid value of Omega');
        end

        sat_in_view_matrix = zeros(length(M0_values), length(M0_values), length(phi_values));

        % Inner loop: iterate over M0 values
        for m = 1:length(M0_values)
            M0 = M0_values(m);

            E0 = (miu/a^3)^0.5*(t-t0)+M0;
            options = optimoptions('fsolve', 'Display', 'none');
            eqns_numeric = @(E) (miu/a^3)^0.5*(t-t0)+M0-(E-e*sin(E));
            E_numeric = fsolve(eqns_numeric, E0, options);
            
            % Relationships between eccentricity anomaly E, satellite-to-centre distance r and true anomaly f
            eqns2 = [r*cos(f) == a*(cos(E_numeric)-e), r*sin(f) == a*(1-e^2)^0.5*sin(E_numeric)]; 
            S = solve(eqns2,[r f]);
            r1 = S.r(1);
            f1 = S.f(1);
        
            u = f1+omega; % Angle between satellite and RAAN at specific time t
            x = vpa(r1*cos(u)); % Satellite coordinates in orbital coodinate system
            y = vpa(r1*sin(u));
            z = 0;
            
            % Orbital to Inertial
            A = Rz(Omega).'*Rx(inc).';
            I = A*[x;y;z];
            
            % Inertial to Fixed
            W = Ephemeris(t, t0);
            B = Rz(W); % Transformation matrix
            F = B*I;

            % Fixed to Horizontal
            % Preallocating spaces
            xh = zeros(size(Phi));
            yh = zeros(size(Phi));
            zh = zeros(size(Phi));
            R = zeros(size(Phi));
            
            % Vectorized computation for coordinates of each satellite 
            % about each longtitude and latitude of observer (continueing "fixed to horizontal")
            for n = 1:numel(Phi)
                phi = Phi(n);
                lambda = Lambda(n);
                
                % Transformation matrix for each phi and lambda
                C = Ry(pi/2-phi)*Rz(lambda);
                
                % Transform and adjust by Moon's mean radius
                P = C*F - [0; 0; 1737400];
            
                % Store results
                xh(n) = P(1);
                yh(n) = P(2);
                zh(n) = P(3);
                R(n) = (xh(n).^2+yh(n).^2+zh(n).^2)^0.5;
                H_matrix(m, :, n) = [xh(n)/R(n), yh(n)/R(n), zh(n)/R(n), 1];

                Ele_sat = atan(zh(n)/(xh(n)^2+yh(n)^2)^0.5);
                Ele_sat_deg = Ele_sat*180/pi;

                if Ele_sat > Ele
                    sat_in_view_matrix(m, m, n) = 1;
                end
            end
        end

        % Store matrices for this "o"th loop in cell arrays
        H_matrices{o} = H_matrix;
        sat_in_view_matrices{o} = sat_in_view_matrix;
    end

    % Convert the cell matrix to a numeric matrix
    sat_in_view_cell_array = {sat_in_view_matrices{1}, zeros(4, 4, n); zeros(4, 4, n), sat_in_view_matrices{2}};
    sat_in_view = cell2mat(sat_in_view_cell_array);
    % 3D version of matrix calculation H=(cell2mat(H_matrices).'*sat_in_view).'
    % The code "cell2mat(H_matrices)" concatenate "o" H matrices (which are 3D matrices, "o" is the number of "o" loops)
    % into one large 3D matrix, along the direction of column. Hence the H matrix becomes 8*4*n.
    % Mutiplied by satellite visibility matrices, some rows in it will be replaced to zeros.
    H = pagetranspose(pagemtimes(pagetranspose(cell2mat(H_matrices)), sat_in_view));
    
    % Create a logical matrix and sum along columns
    nonZeroElementsPerRow = sum(sat_in_view ~= 0, 2);
    % Check for non-zero counts
    nonZeroRows = nonZeroElementsPerRow > 0;
    % Sum the result to get the number of non-zero rows
    numberOfSatInView = reshape(sum(nonZeroRows), 1, n);

    AA = pageinv((pagemtimes(pagetranspose(H), H))); % Calculate A matrix from H matrix

    % If the number of visible satellites < 4, the DOP value is meaningless
    % and set to NaN
    for j = 1:length(numberOfSatInView)
        if numberOfSatInView(1, j) < 4
            AA(:, :, j) = NaN(size(AA(:, :, j)));
        end
    end
    
    % Calculate DOPs for this "i"th loop and store them in 3d matrices
    GDOP_array = (AA(1,1,:)+AA(2,2,:)+AA(3,3,:)+AA(4,4,:)).^0.5; % This results to 1*1*n matrices
    GDOP(:, :, i) = reshape(GDOP_array, [length(phi_values), length(lambda_values)]); % Flat the 3rd dimension, time loops become the new 3rd dimension
    PDOP_array = (AA(1,1,:)+AA(2,2,:)+AA(3,3,:)).^0.5;
    PDOP(:, :, i) = reshape(PDOP_array, [length(phi_values), length(lambda_values)]);
    HDOP_array = (AA(1,1,:)+AA(2,2,:)).^0.5;
    HDOP(:, :, i) = reshape(HDOP_array, [length(phi_values), length(lambda_values)]);

end

GDOPvsTime = GDOP(length(phi_values), length(lambda_values)/2, :);
GDOPvsTimeArray = squeeze(GDOPvsTime);
PDOPvsTime = PDOP(length(phi_values), length(lambda_values)/2, :);
PDOPvsTimeArray = squeeze(PDOPvsTime);
HDOPvsTime = HDOP(length(phi_values), length(lambda_values)/2, :);
HDOPvsTimeArray = squeeze(HDOPvsTime);

% This is to avoid oversized values affecting color display
GDOP_t = GDOP(:,:,t_DOP);
for k=1:numel(GDOP_t)
    if GDOP_t(k) > 20
        GDOP_t(k) = NaN;
    end
end

PDOP_t = PDOP(:,:,t_DOP);
for k=1:numel(PDOP_t)
    if PDOP_t(k) > 20
        PDOP_t(k) = NaN;
    end
end

HDOP_t = HDOP(:,:,t_DOP);
for k=1:numel(HDOP_t)
    if HDOP_t(k) > 20
        HDOP_t(k) = NaN;
    end
end

GDOP_max = max(GDOP, [], 3);
for k=1:numel(GDOP_max)
    if GDOP_max(k) > 20
        GDOP_max(k) = NaN;
    end
end

PDOP_max = max(PDOP, [], 3);
for k=1:numel(PDOP_max)
    if PDOP_max(k) > 20
        PDOP_max(k) = NaN;
    end
end

HDOP_max = max(HDOP, [], 3);
for k=1:numel(HDOP_max)
    if HDOP_max(k) > 20
        HDOP_max(k) = NaN;
    end
end

toc

%%
%----------------GDOP Plot--------------------------------------------------
figure;
h = imagesc(lambda_values*180/pi, phi_values*180/pi, GDOP_t);
axis image on;
set(h, 'AlphaData', ~isnan(GDOP_t)); % If a value in "GDOP_t" is "NaN", its color will be transparent
set(gca, 'YDir', 'normal'); % Adjust the Y-axis direction
set(gca,'color',[1 1 1])
xlabel('Longitude [°]');
ylabel('Latitude [°] ');
title(['GDOP at t=', num2str(t_values(t_DOP)), ' [s]']);
colorbar;
grid on;

figure;
h = imagesc(lambda_values*180/pi, phi_values*180/pi, GDOP_max);
axis image on;
set(h, 'AlphaData', ~isnan(GDOP_max)); % If a value in "GDOP_max" is "NaN", its color will be transparent
set(gca, 'YDir', 'normal'); % Adjust the Y-axis direction
set(gca,'color',[1 1 1])
xlabel('Longitude [°]');
ylabel('Latitude [°] ');
title(['Maximum GDOP in a Period of ', num2str(max(t_values)/3600), ' Hours']);
colorbar;
grid on;
%------------------------PDOP Plot----------------------------------------------
figure;
h = imagesc(lambda_values*180/pi, phi_values*180/pi, PDOP_t);
axis image on;
set(h, 'AlphaData', ~isnan(PDOP_t)); % If a value in "PDOP_t" is "NaN", its color will be transparent
set(gca, 'YDir', 'normal'); % Adjust the Y-axis direction
set(gca,'color',[1 1 1])
xlabel('Longitude [°]');
ylabel('Latitude [°] ');
title(['PDOP at t=', num2str(t_values(t_DOP)), ' [s]']);
colorbar;
grid on;

figure;
h = imagesc(lambda_values*180/pi, phi_values*180/pi, PDOP_max);
axis image on;
set(h, 'AlphaData', ~isnan(PDOP_max)); % If a value in "PDOP_max" is "NaN", its color will be transparent
set(gca, 'YDir', 'normal'); % Adjust the Y-axis direction
set(gca,'color',[1 1 1])
xlabel('Longitude [°]');
ylabel('Latitude [°] ');
title(['Maximum PDOP in a Period of ', num2str(max(t_values)/3600), ' Hours']);
colorbar;
grid on;
%---------------------------HDOP Plot------------------------------------------
figure;
h = imagesc(lambda_values*180/pi, phi_values*180/pi, HDOP_t);
axis image on;
set(h, 'AlphaData', ~isnan(HDOP_t)); % If a value in "HDOP_t" is "NaN", its color will be transparent
set(gca, 'YDir', 'normal'); % Adjust the Y-axis direction
set(gca,'color',[1 1 1])
xlabel('Longitude [°]');
ylabel('Latitude [°] ');
title(['HDOP at t=', num2str(t_values(t_DOP)), ' [s]']);
colorbar;
grid on;

figure;
h = imagesc(lambda_values*180/pi, phi_values*180/pi, HDOP_max);
axis image on;
set(h, 'AlphaData', ~isnan(HDOP_max)); % If a value in "HDOP_max" is "NaN", its color will be transparent
set(gca, 'YDir', 'normal'); % Adjust the Y-axis direction
set(gca,'color',[1 1 1])
xlabel('Longitude [°]');
ylabel('Latitude [°] ');
title(['Maximum HDOP in a Period of ', num2str(max(t_values)/3600), ' Hours']);
colorbar;
grid on;
%----------------------------Sphere plot--------------------------------
% P = GDOP_max;
% figure;
% [xx3,yy3,zz3,cc] = sphere3d(P,-pi,pi,-pi/2,pi/2,1,2,'off','spline',.001);
% surf(xx3,yy3,zz3,(cc-1)/0.001); 
% axis off;
% grid off;
% set(gca,'DataAspectRatio',[1 1 1])    
% set(gca,'XDir','normal','YDir','normal'); 
% colorbar;
%-----------------------DOP vs Time Plot----------------------------------------------
figure;
plot(t_values/3600, GDOPvsTimeArray, t_values/3600, PDOPvsTimeArray, t_values/3600, HDOPvsTimeArray);
xlabel('Time [Hr]');
ylabel('DOP');
title('DOPs vs. Time');
legend('GDOP','PDOP','HDOP');
grid on;
%-----------------------------------------------------------------------------------

function [W] = Ephemeris(t, t0)
    % The ephemeris position of the prime meridian (at day d from J2000.0)
    % d=0 is equivalant to 12:00:00, January 1, 2000, or Julian date 2451545.0
    d = t0-2451545+t/3600/24;
    E1 = (125.045-0.0529921*d)*pi/180;
    E2 = (250.089-0.1059842*d)*pi/180;
    E3 = (260.008+13.0120009*d)*pi/180;
    E4 = (176.625+13.3407154*d)*pi/180;
    E5 = (357.529+0.9856003*d)*pi/180;
    E6 = (311.589+26.4057084*d)*pi/180;
    E7 = (134.963+13.0649930*d)*pi/180;
    E8 = (276.617+0.3287146*d)*pi/180;
    E9 = (34.226+1.7484877*d)*pi/180;
    E10 = (15.134-0.1589763*d)*pi/180;
    E11 = (119.743+0.0036096*d)*pi/180;
    E12 = (239.961+0.1643573*d)*pi/180;
    E13 = (25.053+12.9590088*d)*pi/180;
    W = (38.3213 ...
        +13.17635815*d ...
        -1.4*10^-12*d^2 ...
        +3.5610*sin(E1) ...
        +0.1208*sin(E2) ...
        -0.0642*sin(E3) ...
        +0.0158*sin(E4) ...
        +0.0252*sin(E5) ...
        -0.0066*sin(E6) ...
        -0.0047*sin(E7) ...
        -0.0046*sin(E8) ...
        +0.0028*sin(E9) ...
        +0.0052*sin(E10) ...
        +0.0040*sin(E11) ...
        +0.0019*sin(E12) ...
        -0.0044*sin(E13)) ...
        *pi/180;
end

function [rx] = Rx(x)
    rx = [1 0 0;0 cos(x) sin(x);0 -sin(x) cos(x)];
end
function [ry] = Ry(x)
    ry = [cos(x) 0 -sin(x);0 1 0;sin(x) 0 cos(x)];
end
function [rz] = Rz(x)
    rz = [cos(x) sin(x) 0;-sin(x) cos(x) 0;0 0 1];
end
