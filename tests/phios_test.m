%% phios_test.m

% Test input parameters
k = [0; 1; 2];          % Wavenumbers [rad/m]
D = 10;                  % Flexural rigidity [Nm]
DEL = [1000; 500];       % Density contrasts [kg/m^3] (surface, subsurface)
g = 9.81;                % Gravity [m/s^2]

expected_phi = 1 + D * k.^4 / (g * DEL(1));
phi = phios(k,D,DEL,g);

tolerance = 1e-12;
if all(abs(phi - expected_phi) < tolerance)
    disp('phios_test passed: numerical values correct.')
else
    error('phios_test failed: numerical values incorrect.')
end

