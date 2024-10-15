%% A452 Project 1 Rendezvous
% Jacob Amezquita
% Ryan Pettigrew

%Constants
RE = 6378;       %km, radius Earth
muE = 398600;    %km^3/s^2, Mu earth
mE = 5.9722e24;  %kg, mass Earth

% S/C A
ecc_A = 0.0002046; 
inc_A = 0.0921;    % degs
rp_A = 35777 + RE; % km 
ra_A = 35794 + RE; % km 
raan_A = 79.9111;  % degs 
omega_A = 98.1606; % degs 
T_A = 1.00272835;  % rev/day
ma_A = 210.4111;   % degs, mean anomaly

a_A = (ra_A+rp_A)/2;
h_A = sqrt(a_A*muE*(1-ecc_A^2));

% S/C B
ecc_B = 0.0002046; 
inc_B = 0.0921;    % degs
rp_B = 35777 + RE;      % km 
ra_B = 35794 + RE;      % km 
raan_B = 79.9111;  % degs 
omega_B = 98.1606; % degs 
T_B = 1.00272835;  % rev/day
ma_B = 210.4111;   % degs, mean anomaly

a_B = (ra_B+rp_B)/2;
h_B = sqrt(a_B*muE*(1-ecc_B^2));


% target parameters
[~,~,rECI_A,vECI_A] = coes2rv(ecc_A,h_A,inc_A,raan_A,omega_A,theta_A,mu); % km & km/s, target r&v
a_A = h_A^2 / (mu * (1-ecc_A^2)); % km, target semi-major axis
P_A = 2*pi*a_A^(3/2) / sqrt(mu); % s, target period
hECI_A0 = cross(rECI_A,vECI_A); % km2/s

% chaser parameters
[~,~,rECI_B,vECI_B] = coes2rv(ecc_B,h_B,inc_B,raan_B,omega_B,theta_B,mu); % km & km/s, chaser r&v


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

% Stumpff Function
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




