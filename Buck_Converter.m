close all; clear; clc,;

% Reqiured parameters to design a power converter

Vin = 12;                   % Input Voltage
Vout = 6;                   % Output Voltage
Iout = 2;                   % Output Current
i_rip_per = 10;             % ripple current in percent
v_rip_per = 0.2;            % voltage ripple in percent
fsw = 100e3;                % switching frequency
eff_est = 92;               % estimated efficiency
Rds_HS = 0.01;              % Rds on of High side Mosfet
Rds_LS = 0.01;              % Rds on of Low side Mosfet
switching_loss_est = 10;    % Estimated switching loss

% Structure:
Params = struct('Vin', Vin, ...
                'Vout', Vout, ...
                'Iout', Iout, ...
                'i_rip_per', i_rip_per, ...
                'v_rip_per', v_rip_per, ...
                'fsw', fsw, ...
                'eff_est', eff_est, ...
                'Rds_HS', Rds_HS, ...
                'Rds_LS', Rds_LS, ...
                'switching_loss_est', switching_loss_est);

% Clear everything except the Params structure
clearvars -except Params

%%
[Params.L, Params.rL, Params.C, Params.rC, Params.Rout, Params.Duty] = Power_Components_Cal(Params);


%% Controller Design

% Design Goals
Params.fc = Params.fsw / 10;            % Target Crossover Frequency (10 kHz)
Params.Vref = Params.Duty;             % Controller Reference Voltage
[plant, designs] = model_buck_plant(Params);

%%

[Comp_III, Comp_II, loops] = design_compensators(plant, designs, Params);

%%

components = calculate_type3_components(designs, Params, Comp_III, plant);


%%

% Generate Bode plots, Step responses, and Stability Margins
analyze_results(plant, Comp_III, Comp_II, loops, designs, Params);


%%

function [L, rL, C, rC, Rout, Duty] = Power_Components_Cal(p)

eff = p.eff_est / 100;
Ploss = (p.Vout * p.Iout) * ((1 - eff)/eff);
Ploss_HS = p.Iout^2 * p.Rds_HS * (p.Vout / p.Vin);
Ploss_LS = p.Iout^2 * p.Rds_LS * (1 - (p.Vout / p.Vin));
Ploss = Ploss - Ploss_HS - Ploss_LS;
Ploss = Ploss - (Ploss * (p.switching_loss_est/100));
Duty = p.Vout / (p.Vin * eff);
i_ripple = p.Iout * (p.i_rip_per / 100);
v_ripple = p.Vout * (p.v_rip_per / 100);

L = (p.Vout * (1 - Duty))/(i_ripple * p.fsw);
C = (1 - Duty)/(8 * L * (p.v_rip_per / 100) * p.fsw^2);
rC = v_ripple / i_ripple;
Rout = p.Vout / p.Iout;
rL = abs((Ploss / p.Iout^2));
end



function [plant, designs] = model_buck_plant(p)
    % A. Critical Frequencies
    designs.w_sw  = 2 * pi * p.fsw;
    designs.w_c   = 2 * pi * p.fc;
    designs.f0    = 1 / (2 * pi * sqrt(p.L * p.C));     % Resonance (Hz)
    designs.fz_esr= 1 / (2 * pi * p.rC * p.C);          % ESR Zero (Hz)
    
    % Radians
    designs.w0    = 2 * pi * designs.f0;
    designs.w_esr = 2 * pi * designs.fz_esr;
    if isinf(designs.w_esr), designs.w_esr = designs.w_sw; end

    % B. State-Space Matrices (Avg Model with Parasitics)
    R_prime = p.Rout + p.rC;
    
    A = [ -((p.Rout * p.rC)/(p.L * R_prime) + p.rL/p.L),  -p.Rout/(p.L * R_prime);
           p.Rout/(p.C * R_prime),                        -1/(p.C * R_prime) ];
       
    B = [ p.Vin / p.L;  0 ];        % Control-to-Output Input
    C = [ (p.Rout * p.rC)/R_prime,  p.Rout/R_prime ];
    D = 0;

    plant = ss(A, B, C, D);
end



function [C_III, C_II, loops] = design_compensators(plant, d, p)
    % --- Type III Design (Voltage Mode) ---
    % Zeros at w0 (cancel LC), Pole 1 at ESR, Pole 2 at fsw/2
    % *Note: Preserved your specific tuning offset (+120k) for Pole 1
    w_z1 = d.w0;  
    w_z2 = d.w0;
    w_p1 = d.w_esr + 2*pi*19e3; % Adjusted to match your text file logic approx
    w_p2 = d.w_sw / 2;

    Shape_III = zpk([-w_z1, -w_z2], [0, -w_p1, -w_p2], 1);
    
    % Calculate K for 0dB at fc
    [mag_sys, ~] = bode(plant, d.w_c);
    [mag_comp, ~] = bode(Shape_III, d.w_c);
    K_III = 1 / (mag_sys * mag_comp);
    C_III = K_III * Shape_III;

    % --- Type II Design (For Comparison) ---
    w_z_II = d.w0 * 0.5;
    w_p_II = d.w_sw / 2;
    
    Shape_II = zpk(-w_z_II, [0, -w_p_II], 1);
    [mag_comp_II, ~] = bode(Shape_II, d.w_c);
    K_II = 1 / (mag_sys * mag_comp_II);
    C_II = K_II * Shape_II;

    % --- Create Loops ---
    loops.L_III = plant * C_III;
    loops.L_II  = plant * C_II;
    loops.CL_III = feedback(loops.L_III, 1);
    loops.CL_II  = feedback(loops.L_II, 1);
end



function comp = calculate_type3_components(d, p, C_III, plant)
    % Extract poles/zeros from the designed object to ensure match
    [z, pole, ~] = zpkdata(C_III, 'v');
    % Sort magnitude to identify them correctly
    z = sort(abs(z)); 
    pole = sort(abs(pole)); 
    
    % Map to design frequencies
    fz = z(1) / (2*pi);       % Double Zero location
    fp1 = pole(2) / (2*pi);   % First non-zero pole
    fp2 = pole(3) / (2*pi);   % High freq pole
    
    % 1. Arbitrary R1
    R1 = 10e3; 
    
    % 2. Calculate Gain Requirement
    [mag_plant, ~] = bode(plant, d.w_c);
    Desired_Gain = 1 / mag_plant;

    % 3. Calculate Components
    C1 = 1 / (2 * pi * fz * Desired_Gain * R1);
    R2 = 1 / (2 * pi * fz * C1);
    C2 = 1 / (2 * pi * fp1 * R2);
    
    % Input Feedforward
    C3 = (fp2 - fz) / (2 * pi * R1 * fp2 * fz);
    R3 = 1 / (2 * pi * C3 * fp2);
    
    % DC Setpoint
    R4 = R1 * p.Vref / (p.Vout - p.Vref);

    % Store in struct
    comp.R1 = R1; comp.R2 = R2; comp.R3 = R3; comp.R4 = R4;
    comp.C1 = C1; comp.C2 = C2; comp.C3 = C3;
    comp.fz = fz; comp.fp1 = fp1; comp.fp2 = fp2;
end


function analyze_results(plant, C_III, C_II, loops, d, p)
    opts = bodeoptions; 
    opts.FreqUnits = 'Hz'; 
    opts.Grid = 'on'; 
    opts.PhaseMatching = 'on';

    % --- Fig 1: Plant Frequency Response with Markers ---
    figure('Name', 'Buck Plant Gvd');
    margin(plant, opts);
    title(sprintf('Plant Gvd (rC=%.3f, rL=%.3f)', p.rC, p.rL));
    
    % Add Markers (using logic similar to your original code)
    hold on;
    ax = findall(gcf, 'type', 'axes');
    for k = 1:length(ax)
        axes(ax(k)); %#ok<LAXES>
        xline(d.f0, '--r', 'LC Res');
        xline(d.fz_esr, '--b', 'ESR Zero');
    end
    hold off;

    % --- Fig 2: Stability Comparison (Open Loop) ---
    figure('Name', 'Stability Analysis');
    margin(loops.L_II, opts);
    hold on;
    margin(loops.L_III, opts);
    legend('Type II (Current)', 'Type III (Voltage)');
    title('Loop Stability Comparison');
    hold off;

    figure('Name', 'Controller Analysis');
    margin(C_II, opts);
    hold on;
    margin(C_III, opts);
    legend('Type II (Current)', 'Type III (Voltage)');
    title('Loop Stability Comparison');
    hold off;

    % --- Fig 3: Transient Response (Closed Loop) ---
    figure('Name', 'Step Response');
    step(loops.CL_II);
    hold on;
    step(loops.CL_III);
    legend('Type II', 'Type III');
    grid on;
    title('Load Step Response');

    % --- Print Margins ---
    [Gm3, Pm3, ~, Wcp3] = margin(loops.L_III);
    [Gm2, Pm2, ~, Wcp2] = margin(loops.L_II);
    
    fprintf('\n=== Final Stability Results ===\n');
    fprintf('Type III: PM = %.2f deg @ %.2f kHz (GM = %.1f dB)\n', Pm3, Wcp3/2/pi/1000, 20*log10(Gm3));
    fprintf('Type II : PM = %.2f deg @ %.2f kHz (GM = %.1f dB)\n', Pm2, Wcp2/2/pi/1000, 20*log10(Gm2));
end