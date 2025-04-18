% Branden Stahl
% AE 313 Spr 25
clear;
clc;
close all;

% Constants
mu = 398600;
Re = 6378;

% Gauss Refinement Threshold
Epsilon = power(10, -8);

% Determining which values are being used:
textbookUse = NaN;
while textbookUse ~= 0 && textbookUse ~= 1
    in1 = lower(input("Use values from textbook? (y/n t/f 1/0 true/false yes/no)\n", "s"));

    if in1 == "true" || in1 == "yes" || in1 == "y" || in1 == "1" || in1 == "t"
        textbookUse = 1;
        fprintf("Using textbook options...\n")
    elseif in1 == "false" || in1 == "no" || in1 == "n" || in1 == "0" || in1 == "f"
        textbookUse = 0;
        fprintf("Using given options...\n")
    else
        fprintf("You have not entered a valid answer. Try again.\n")
    end
end

if textbookUse == 1
    % Use textbook
    phi = deg2rad(40);
    alt = 1;

    t = [0, 118.1, 237.58];
    alphai = deg2rad([43.537, 54.42, 64.318]);
    deltai = deg2rad([-8.7833, -12.074, -15.105]);
    thetai = deg2rad([44.506, 45, 45.499]);
else
    % Use given
    phi = deg2rad(60);
    alt = 0.5;

    t = 60*[0, 5, 10];
    alphai = deg2rad([157.783, 159.221, 160.526]);
    deltai = deg2rad([24.2403, 27.2993, 29.8982]);
    thetai = deg2rad([150, 151.253, 152.507]);
end
Rmag = Re + alt;

% Solve for LOS Unit Vectors (Yes I misspelled Rho, I use Rho later)
ro1 = sphericalize(deltai(1), alphai(1));
ro2 = sphericalize(deltai(2), alphai(2));
ro3 = sphericalize(deltai(3), alphai(3));

% Solving for Observer ECL position
R1 = Rmag*sphericalize(phi, thetai(1));
R2 = Rmag*sphericalize(phi, thetai(2));
R3 = Rmag*sphericalize(phi, thetai(3));
function ro = sphericalize(a, b)
    ro = [cos(a)*cos(b); cos(a)*sin(b); sin(a)];
end

% Solving for time intervals
T1 = t(1)-t(2);
T3 = t(3)-t(2);
T = T3-T1;

% Separation between iterative and non iterative
%% NON-ITERATIVE
% Solving for 3 position vectors
R = [R1, R2, R3];
ro = [ro1, ro2, ro3];
M = ro\R;

% Solve for A, B, E
A = (M(2,1)*(T3/T)) - M(2,2) - (M(2,3)*(T1/T));
B = M(2,1)*(T3/(6*T))*(T^2-T3^2) - M(2,3)*(T1/(6*T))*(T^2-T1^2);
E = dot(R2, ro2);

% Solve for a, b, c as it relates to 8th Degree Polynomial
a = -((A^2) + 2*A*E + (norm(R2)^2));
b = -2*mu*B*(A+E);
c = -(mu^2)*(B^2);

% Solve for initial r2 by roots
r2mag_roots = roots([1, 0, a, 0, 0, b, 0, 0, c]);
r2mag = r2mag_roots(imag(r2mag_roots) == 0 & real(r2mag_roots) > 0);

% Scaling constants
c1 = (T3/T)*(1+(mu/(6*r2mag^3))*(T^2-T3^2));
c3 = -(T1/T)*(1+(mu/(6*r2mag^3))*(T^2-T1^2));

% Black magic matrix operations (from slides)
RPM = M*[-c1; 1; -c3];
p1 = RPM(1)/c1;
p2 = RPM(2)/-1;
p3 = RPM(3)/c3;

% Solving assumed r1, r2, r3
r1 = R1 + p1*ro1;
r2 = R2 + p2*ro2;
r3 = R3 + p3*ro3;

% Solving lagrange
f1 = 1 - (mu/(2*norm(r2)^3))*T1^2;
f3 = 1 - (mu/(2*norm(r2)^3))*T3^2;
g1 = T1 - (mu/(6*norm(r2)^3))*T1^2;
g3 = T3 - (mu/(6*norm(r2)^3))*T3^2;

% Option one, solving for v2
num = (f1*r3) - (f3*r1);
den = (f1*g3) - (f3*g1);
v2_1 = num/den;
NI_O1_OE = RV2OE(r2, v2_1);

% Option two, using Gibbs
N = norm(r1)*cross(r2, r3) + norm(r2)*cross(r3, r1) + norm(r3)*cross(r1, r2);
D = cross(r1, r2) + cross(r2, r3) + cross(r3, r1);
S = r1*(norm(r2)-norm(r3)) + r2*(norm(r3)-norm(r1)) + r3*(norm(r1)-norm(r2));
v2_2 = (cross(D, r2/norm(r2)) + S)*sqrt(mu/dot(N,D));
NI_O2_OE = RV2OE(r2, v2_2);

%% ITERATIVE
% Clear all vars except project based to avoid using incorrect vars
clearvars -except phi alphai deltai t T T1 T3 Rmag R1 R2 R3 ro1 ro2 ro3 mu Re thetai v2_1 v2_2 NI_O1_OE NI_O2_OE Epsilon

% Getting cross products
p1 = cross(ro2, ro3);
p2 = cross(ro1, ro3);
p3 = cross(ro1, ro2);

% Solving for scalar quantities
D0 = dot(ro1, p1);
D = [dot(R1, p1), dot(R1, p2), dot(R1, p3); dot(R2, p1), dot(R2, p2), dot(R2, p3); dot(R3, p1), dot(R3, p2), dot(R3, p3)];

% Solving A & B
A = (1/D0)*(-D(1,2)*(T3/T) + D(2,2) + D(3,2)*(T1/T));
B = (1/(6*D0))*(D(1,2)*((T3^2)-(T^2))*(T3/T) + D(3,2)*((T^2)-(T1^2))*(T1/T));

% Solving for E
E = dot(R2, ro2);

% Solving for a, b, c
a = -((A^2) + 2*A*E + (norm(R2)^2));
b = -2*mu*B*(A+E);
c = -(mu^2)*(B^2);

% 8th Degree Polynomial
% x^8 + a*x^6 + b*x^3 + c
r2mag_roots = roots([1 0 a 0 0 b 0 0 c]);
r2mag = r2mag_roots(imag(r2mag_roots) == 0 & real(r2mag_roots) > 0);

% Get Rho Magnitudes
rho1 = (1/D0)*((6*(D(3,1)*T1/T3 + D(2,1)*T/T3)*r2mag^3 + mu*D(3,1)*(T^2 - T1^2)*T1/T3)/(6*r2mag^3 + mu*(T^2 - T3^2)) - D(1,1));
rho2 = A + mu*B/r2mag^3;
rho3 = 1/D0*((6*(D(1,3)*T3/T1 - D(2,3)*T/T1)*r2mag^3 + mu*D(1,3)*(T^2 - T3^2)*T3/T1)/(6*r2mag^3 + mu*(T^2 - T1^2)) - D(3,3));

% Get new position vectors
r1 = R1 + rho1*ro1;
r2 = R2 + rho2*ro2;
r3 = R3 + rho3*ro3;

% Solving lagrange
f1 = 1 - (mu/(2*norm(r2)^3))*T1^2;
f3 = 1 - (mu/(2*norm(r2)^3))*T3^2;
g1 = T1 - (mu/(6*norm(r2)^3))*T1^2;
g3 = T3 - (mu/(6*norm(r2)^3))*T3^2;

% Solving for v2
num = (f1*r3) - (f3*r1);
den = (f1*g3) - (f3*g1);
v2 = num/den;

% Open Iterative Loop
i = 0;
conv = 1/Epsilon; % Guarantees it's always higher initially when Epsilon < 1
while conv > Epsilon && i < 100
    i = i + 1; % Increase iteration index

    % Algorithm 5.6; Iterative Refinement
    r2mag = norm(r2);
    v2mag = norm(v2);

    % Semi-Major Axis Reciprocal
    alpha = 2/r2mag - (v2mag^2)/mu;
    
    % Radial Velocity 2
    vr2 = dot(v2, r2/r2mag);
    
    % Solve the Kepler equation (Woah this took too long)
    uni1 = Kepler(T1, r2mag, vr2, alpha, Epsilon);
    uni3 = Kepler(T3, r2mag, vr2, alpha, Epsilon);
    
    % Lagrange Redefinition
    f1 = (1 - ((uni1^2)/r2mag)*Cf(alpha*uni1^2) + f1)/2;
    f3 = (1 - ((uni3^2)/r2mag)*Cf(alpha*uni3^2) + f3)/2;
    g1 = (T1 - (1/sqrt(mu))*(uni1^3)*Sf(alpha*uni1^2) + g1)/2;
    g3 = (T3 - (1/sqrt(mu))*(uni3^3)*Sf(alpha*uni3^2) + g3)/2;
    
    % Scaling Factors
    c1 = g3/(f1*g3 - f3*g1);
    c3 = -g1/(f1*g3 - f3*g1);
    
    % Save old values to compare convergence
    rho1old = rho1;
    rho2old = rho2;
    rho3old = rho3;
    
    % Set new rho values
    rho1 = (1/D0)*(-D(1,1) + D(2,1)/c1 - c3*D(3,1)/c1);
    rho2 = (1/D0)*(-c1*D(1,2) + D(2,2) - c3*D(3,2));
    rho3 = (1/D0)*((-c1/c3)*D(1,3) + D(2,3)/c3 - D(3,3));

    % Positions
    r1 = R1 + rho1*ro1;
    r2 = R2 + rho2*ro2;
    r3 = R3 + rho3*ro3;

    % Get new v2
    num = (f1*r3) - (f3*r1);
    den = (f1*g3) - (f3*g1);
    v2 = num/den;

    conv = (abs(rho1old-rho1)+abs(rho2old-rho2)+abs(rho3old-rho3))/3; % Convergence
end
fprintf("Finished solving Iterative Gauss Method after %d iterations.\n", i);

% Get output OE
I_OE = RV2OE(r2, v2);

% Solve the Universal Kepler Equation for the Universal Variables
function uni = Kepler(dT, r, vr, alpha, Epsilon)
    mu = 398600;
    X0 = sqrt(mu)*abs(alpha)*dT; % Initial Guess
    rat = 1/Epsilon; % Create an initial value to enter loop

    i = 0;
    while abs(rat) > Epsilon && i < 10
        i = i + 1;

        z = alpha*power(X0,2);
        
        fun = (r*vr/sqrt(mu))*power(X0,2)*Cf(z) + (1-alpha*r)*power(X0,3)*Sf(z) + r*X0 - sqrt(mu)*dT;
        funprime = (r*vr/sqrt(mu))*X0*(1-z*Sf(z)) + (1-alpha*r)*power(X0, 2)*Cf(z) + r;

        rat = fun/funprime;

        if abs(rat) > Epsilon
            X0 = X0 - rat;
        end
    end

    uni = X0;
end
function VAL = Sf(z) % Textbook S(z) function
    if z > 0
        VAL = (sqrt(z)-sin(sqrt(z)))/power(sqrt(z), 3);
    elseif z < 0
        VAL = (sinh(sqrt(-z))-sqrt(-z))/power(sqrt(z), 3);
    else
        VAL = 1/6;
    end
end
function VAL = Cf(z) % Textbook C(z) function
    if z > 0
        VAL = (1-cos(sqrt(z)))/z;
    elseif z < 0
        VAL = (cosh(sqrt(-z))-1)/(-z);
    else
        VAL = 1/2;
    end
end
function OE = RV2OE(R, V)
    mu = 398600;

    % Solve for generic elements
    r = norm(R);

    H = cross(R, V);
    h = norm(H);

    E = cross(V, H)/mu - R/r;
    e = norm(E);
    
    N = cross([0, 0, 1], H);
    n = norm(N);

    i = acos(dot(H, [0 0 1])/h);
    i = rad2deg(i);

    % Solve for angles and ambiguity
    theta = acos(dot(R, E)/(e*r));
    if dot(R, V) < 0
        theta = 2*pi-theta;
    end
    theta = rad2deg(theta);
    RAAN = acos(N(1)/n);
    if N(2) < 0
        RAAN = 2*pi - RAAN;
    end
    RAAN = rad2deg(RAAN);
    AP = acos(dot(N, E)/(n*e));
    if E(3) < 0
        AP = 2*pi - AP;
    end
    AP = rad2deg(AP);
    
    OE = [h, e, i, RAAN, AP, theta];
end

% Printing, intentionally ending code at 313
fprintf("\nNon-Iterative Gauss:\n   h = %.4f\n   e = %.4f\n   inclination = %.4f\n   RAAN = %.4f\n   AP = %.4f\n   true anomaly = %.4f\n", NI_O1_OE); fprintf("\nGibbs:\n   h = %.4f\n   e = %.4f\n   inclination = %.4f\n   RAAN = %.4f\n   AP = %.4f\n   true anomaly = %.4f\n", NI_O2_OE); fprintf("\nIterative Gauss:\n   h = %.4f\n   e = %.4f\n   inclination = %.4f\n   RAAN = %.4f\n   AP = %.4f\n   true anomaly = %.4f\n", I_OE);