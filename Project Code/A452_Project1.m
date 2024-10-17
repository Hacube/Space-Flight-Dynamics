%% A452 Project 1 Rendezvous
% Jacob Amezquita
% Ryan Pettigrew

clc; close all; clear;

%Constants
RE = 6378;      %km, radius Earth
mu = 398600;    %km^3/s^2, Mu earth

% S/C A
ecc_A = 0.0002046; 
inc_A = 0.0921;     % degs
rp_A = 35777 + RE;  % km 
ra_A = 35794 + RE;  % km 
raan_A = 79.9111;   % degs 
omega_A = 98.1606;  % degs 
T_A = 1.00272835;   % rev/day
theta_A = 210.4111; % degs, mean anomaly

a_A = (ra_A+rp_A)/2;            % km
h_A = sqrt(a_A*mu*(1-ecc_A^2)); % km^2/s

% S/C B
ecc_B = 0.0002046; 
inc_B = 0.0921;    % degs
rp_B = ra_A;       % km 
ra_B = ra_A;       % km 
raan_B = 79.9111;  % degs 
omega_B = 98.1606; % degs 
T_B = 1.00272835;  % rev/day
theta_B = 210.577; % degs, mean anomaly starting 100 km ahead

a_B = (ra_B+rp_B)/2;            % km
h_B = sqrt(a_B*mu*(1-ecc_B^2)); % km^2/s

% mean motion
n = sqrt(mu/a_A^3);  % rad/s

% target parameters
[~,~,rECI_A,vECI_A] = coes2rv(ecc_A,h_A,inc_A,raan_A,omega_A,theta_A,mu); % km & km/s, target r&v
a_A = h_A^2 / (mu * (1-ecc_A^2)); % km, target semi-major axis
P_A = 2*pi*a_A^(3/2) / sqrt(mu);  % s, target period
hECI_A0 = cross(rECI_A,vECI_A);   % km2/s

% chaser parameters
[~,~,rECI_B,vECI_B] = coes2rv(ecc_B,h_B,inc_B,raan_B,omega_B,theta_B,mu); % km & km/s, chaser r&v

% Step forward (10 periods)
dt = 1;            % s, time step
TOL = 10e-8;       % tolerance for UV
countMax = 1000;   % max iterations for UV
tf = 10*P_A;       % s, final time from start
num = ceil(tf/dt)+1; % number of time steps
t = 0;             % pre-allocate time variable
rho = zeros(num,5);  % pre-allocate matrix for time and rho

hops = zeros(num,7);  % pre-allocate matrix for time and hop

% Mission Target Waypoints
waypt0 = hECI_A0;
waypt1 = hECI_A0 - [0; 100; 0];
% waypt2
t_hop = 2;
t_Fb= 5;
% [time, x, y,z];
wayPoint = [0; waypt0;
            10; 1; 1; 1];

delta_y = 20000;

for i = 1:dt:num
    
    % if i <  wayPoint(1)
    %     % inital position / orbit
    % elseif  i < waypoint(2)
    %     % first maneuver 
    % end

    hECI_A = cross(rECI_A,vECI_A); % km2/s
    aECI_A = (-mu/norm(rECI_A)^3).*rECI_A; % km/s2
    aECI_B = (-mu/norm(rECI_B)^3).*rECI_B; % km/s2

    % 5-term acceleration
    rhoECI = rECI_B - rECI_A; % km, relative positon of chaser
    Omega = hECI_A/(norm(rECI_A)^2); % absolute angular velocity of moving frame
    drhoECI = vECI_B - vECI_A - cross(Omega,rhoECI); % km/s
    dOmega = -2*dot(vECI_A,rECI_A).*Omega ./ norm(rECI_A)^2; % absolute angular acceleration of moving frame
    ddrhoECI = aECI_B - aECI_A - 2.*cross(Omega,drhoECI) - cross(dOmega,rhoECI) + cross(Omega,cross(Omega,rhoECI)); % km/s2
    
    % DCM (ECI to LVLH)
    % ihat = rECI_A / norm(rECI_A);
    % khat = hECI_A / norm(hECI_A);
    % jhat = cross(khat,ihat);
    % Q = [ihat';jhat';khat']; % direction cosine matrix
    Q = ECI2LVLH_DCM(rECI_A, hECI_A); % direction cosine matrix
    
    % rho in LVLH
    rhoLVLH = Q*rhoECI;
    drhoLVLH = Q*drhoECI;
    ddrhoLVLH = Q*ddrhoECI;
    
    % store
    rho(i,1) = t; % s, time step
    rho(i,2:4) = rhoLVLH'; % km, relative position vector
    rho(i,5) = norm(rhoLVLH); % km, relative position magnitude
    
    % re-allocate values
    vECI_A0 = vECI_A;
    
    rECI_A0 = rECI_A;
    vECI_B0 = vECI_B;
    rECI_B0 = rECI_B;
    
    % next step:
    [rECI_A,vECI_A] = UV_rv(dt, vECI_A0, rECI_A0, mu, TOL, countMax);
    [rECI_B,vECI_B] = UV_rv(dt, vECI_B0, rECI_B0, mu, TOL, countMax);
    
    % % hop maneuver
    % [r_final, v_final] = hop(rECI_A, vECI_A, n, t_hop, delta_y);
    % 
    % % store hops
    % hops(i,1) = t; % s, time step
    % hops(i,2:4) = r_final';
    % hops(i,5:7) = v_final';

    t = t + dt;
end

for i = 1:dt:num
    % hop maneuver
    [r_final, v_final] = hop(rECI_A, vECI_A, n, t_hop, delta_y);
    
    % store hops
    hops(i,1) = t; % s, time step
    hops(i,2:4) = r_final';
    hops(i,5:7) = v_final';

    t = t + dt;
end

%% Functions
function [rPeri,vPeri,rECI,vECI] = coes2rv(ecc,h,inc,raan,aop,theta,mu)
    % ecc = eccentricity
    % r = radius, km
    % inc = inclination, degrees
    % raan = RAAN, degrees
    % aop = argument of perigee, degrees
    % theta = true anomaly, degrees
    % Put angles in radians for calculations:
    
    aop = deg2rad(aop);
    inc = deg2rad(inc);
    raan = deg2rad(raan);
    theta = deg2rad(theta);
    
    rPeri = h^2/(mu*(1+ecc*cos(theta)))*[cos(theta);sin(theta);0]; % km, radius vector
    vPeri = (mu/h)*[-sin(theta);ecc+cos(theta);0]; % km/s, velocity vector
    
    R1 = principleRotations('z',aop);
    R2 = principleRotations('x',inc);
    R3 = principleRotations('z',raan);
    
    C_ECI_peri = R1*R2*R3;
    rECI = C_ECI_peri'*rPeri; % km, radius in ECI
    vECI = C_ECI_peri'*vPeri; % km/s, velocity in ECI
end

function [C] = principleRotations(axis,angle)
    % axis as a string 'x','y', or 'z'
    % angle in radians
    
    if axis == 'x'
        C = [1,0,0;0,cos(angle),sin(angle);0,-sin(angle),cos(angle)];
    elseif axis == 'y'
        C = [cos(angle),0,-sin(angle);0,1,0;sin(angle),0,cos(angle)];
    elseif axis == 'z'
        C = [cos(angle),sin(angle),0;-sin(angle),cos(angle),0;0,0,1];
    else
        error("Invalid axis")
    end
end

function [smiley, C, S, alpha] = KeplerUV(dt, v0, r0, mu, TOL, countMax)
    % Find Kepler's Universal Variables
    % dt = time step
    % v0 = velocity vector (km/s)
    % r0 = position vector (km)
    
    r = norm(r0);
    v = norm(v0);
    vr = dot(r0,v0)/r;
    alpha = (2/r)-(v^2)/mu;
    a = 1/alpha;
    smiley = (sqrt(mu)*dt)/abs(a);
    error = 1;
    count = 1;
    
    while (error > TOL && count < countMax)
        [C, S] = stumpff(alpha*smiley^2);
        f = (((r*vr)/sqrt(mu))*smiley^2*C)+((1-(alpha*r))*smiley^3*S)-(sqrt(mu)*dt)+(r*smiley);
        fP = (((r*vr)/sqrt(mu))*smiley*(1-(alpha*smiley^2*S)))+((1-(alpha*r))*(smiley^2*C))+r;
        ratio = f/fP;
        smiley = smiley - ratio;
        error = ratio;
        count = count+1;
    end
end

function [C, S] = stumpff(z)
    if z>0
        C = (1-cos(sqrt(z)))/z;
        S = (sqrt(z)-sin(sqrt(z)))/(sqrt(z))^3;
    elseif z<0
        C = (cosh(sqrt(-z))-1)/(-z);
        S = (sinh(sqrt(-z))-(sqrt(-z)))/(sqrt(-z))^3;
    else % z=0
        S = 1/6;
        C = 1/2;
    end
end

function [rvect,vvect] = UV_rv(dt, v0, r0, mu, TOL, countMax)
    % Finds new r and v vectors for time step using Universal Variable
    
    [smiley, C, S, alpha] = KeplerUV(dt, v0, r0, mu, TOL, countMax);
    
    % Lagrange coefficients
    rn0 = norm(r0);
    f = 1-(smiley^2/rn0)*C;
    g = dt-(1/sqrt(mu))*smiley^3*S;
    
    % a) position
    rvect = (f*r0)+(g*v0);
    rn = norm(rvect);
    fP = (sqrt(mu)/(rn*rn0))*((alpha*(smiley^3)*S)-(smiley));
    gP = 1-((smiley^2/rn)*C);
    
    % b) velocity
    vvect = (fP*r0)+(gP*v0);
end

function [dydt] = linearized_EOM(time,y,mu)
    % Solves linearized EOM for rendezvous and spacecraft separation
    % Use in ODE45
    % y(1:3) - deputy position (km)
    % y(4:6) - deputy velocity (km)
    % y(7:9) - chief position (km)
    % y(10:12) - chief velocity (km)
    
    R = norm(y(7:9)); % km
    h_vect = cross(y(7:9),y(10:12)); % km2/s
    h = norm(h_vect); % km2/s
    rvect = y(7:9); % km, deputy position vector
    vvect = y(10:12); % km/s, deputy velocity vector
    
    % deputy
    dydt(1) = y(4);
    dydt(2) = y(5);
    dydt(3) = y(6);
    dydt(4) = ((2*mu/R^3)+(h^2/R^4))*y(1)-2*(dot(vvect,rvect))*(h/R^4)*y(2)+2*(h/R^2)*y(5);
    dydt(5) = ((-mu/R^3)+(h^2/R^4))*y(2)+2*(dot(vvect,rvect))*(h/R^4)*y(1)-2*(h/R^2)*y(4);
    dydt(6) = -(mu/R^3)*y(3);
    
    % chief
    dydt(7) = y(10);
    dydt(8) = y(11);
    dydt(9) = y(12);
    dydt(10) = -mu*y(7)/R^3;
    dydt(11) = -mu*y(8)/R^3;
    dydt(12) = -mu*y(9)/R^3;
    dydt = dydt';
end

function [r,v] = CW_EOM(r0,v0,n,t)
    % Solves CW equations
    % Assumes no external forces
    % Uses Laplace transform
    % All vectors are column vectors
    % r0 = initial relative position vector (km) [x0;y0;z0]
    % v0 = initial relative velocity vector (km/s) [vx0;vy0;vz0]
    % n = mean motion of target (1/s)
    % t = time from intial (s)
    % r = relative position vector after time t (km) [x;y;z]
    % v = relative velocity vector after time t (km/s) [vx;vy;vz]
    
    x = 4*r0(1) + 2*v0(2)/n + v0(1)/n*sin(n*t) - (3*r0(1) + 2*v0(2)/n)*cos(n*t);
    y = 2*v0(1)/n*cos(n*t) + (6*r0(1) + 4*v0(2)/n)*sin(n*t) - (6*n*r0(1) + 3*v0(2))*t - 2*v0(1)/n + r0(2);
    z = r0(3)*cos(n*t) + v0(3)/n*sin(n*t);
    r = [x;y;z]; % km
    
    vx = v0(1)*cos(n*t) + (3*n*r0(1) + 2*v0(2))*sin(n*t);
    vy = (6*n*r0(1) + 4*v0(2))*cos(n*t) - 2*v0(1)*sin(n*t) - (6*n*r0(1) + 3*v0(2));
    vz = -r0(3)*n*sin(n*t) + v0(3)*cos(n*t);
    v = [vx;vy;vz]; % km/s
end

function [r_f, v_f] = hop(r0, v0, n, t_hop, delta_y)
    % HOP maneuver in CW frame

    % Inputs:
    %   r0: initial pos vector [x; y; z] km
    %   v0: initial vel vector [vx; vy; vz] km/s
    %   n: mean motion of target orbit rad/s
    %   t_hop: Time of the hop maneuver - seconds
    %   delta_y: Desired change in y-position - km
    
    % Outputs:
    %   r_f: final pos after hop [x; y; z] km
    %   v_f: final vel after hop [vx; vy; vz] km/s

    % required dv
    delta_vy = (delta_y - 2*v0(2)/n*(cos(n*t_hop) - 1)) / (4/n*sin(n*t_hop/2)^2);
    v0_new = v0 + [0; delta_vy; 0]; % new inital velocity

    % call CW eqs to propagate
    [r_f, v_f] = CW_EOM(r0, v0_new, n, t_hop); %returns final [pos, vel]
end

function [r_f, v_f] = FB_orbit(r0,v0,n,t_FB,delta_y)
    % Football maneuver in CW frame

    % t_FB: time for football maneuver
    % n: mean motion
    % delta_y: the relative position
    % r0: intial target position
    % v0: initial target velocity
    
    % Required dv
    delta_vx = ((delta_y - r0(2))*n)/(2*(cos(n*t_FB) - 1));
    v0_new = v0 + [delta_vx; 0; 0]; % new inital velocity

    % Call CW eqs to propagate
    [r_f, v_f] = CW_EOM(r0, v0_new, n, t_FB); %returns final [pos, vel]

end

function Q = ECI2LVLH_DCM(rECI_A, hECI_A)
    % DCM (ECI to LVLH)
    ihat = rECI_A / norm(rECI_A);
    khat = hECI_A / norm(hECI_A);
    jhat = cross(khat,ihat);
    Q = [ihat';jhat';khat']; % direction cosine matrix
end