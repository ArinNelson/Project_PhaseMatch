function [omega,chi,f,nu] = tideinfo(tref)

    % Frequencies being used IN RADIANS PER SECOND
    omega_m2 = 1.405189  * 1e-4;	% 0.5175 days per cycle, or 12.4200 hrs
    omega_s2 = 1.454441  * 1e-4;    % 0.5274 days per cycle, or 12.6576 hrs
    omega_n2 = 1.378797  * 1e-4;    % 0.5000 days per cycle, or 12.0000 hrs
    omega_k1 = 0.7292117 * 1e-4;    % 0.9973 days per cycle, or 23.9352 hrs
    omega_o1 = 0.6759774 * 1e-4;    % 1.0758 days per cycle, or 25.8192 hrs
    omega    = [omega_m2 omega_s2 omega_n2 omega_k1 omega_o1];
    omega    = omega .* 86400;      % Units of RADIANS PER DAY
    
    % Compute the ephemerides for the 5 constituents (given in degrees)
    [h0,s0,p0,N0,T0,t0]=ray_arguments(year(tref),day(tref));
    chi_m2 = 2*h0-2*s0;
    chi_s2 = 0;
    chi_n2 = 2*h0-3*s0+p0;
    chi_k1 = h0+90;
    chi_o1 = h0-2*s0-90;
    chi    = [chi_m2 chi_s2 chi_n2 chi_k1 chi_o1];

    % Compute nodal factors f for the 5 constituents:
    f_m2 = 1.000-0.037*cosd(N0);
    f_s2 = 1;
    f_n2 = f_m2;
    f_k1 = 1.006+0.115*cosd(N0);
    f_o1 = 1.009+0.187*cosd(N0);
    f    = [f_m2 f_s2 f_n2 f_k1 f_o1];

    % Compute nodal factors u for the 5 constituents
    nu_m2 = -2.1*sind(N0);
    nu_s2 = 0;
    nu_n2 = nu_m2;
    nu_k1 = -8.9*sind(N0);
    nu_o1 = 10.8*sind(N0);
    nu    = [nu_m2 nu_s2 nu_n2 nu_k1 nu_o1];

end