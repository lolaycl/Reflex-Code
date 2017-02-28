function [I01 TAU31 TAU41 ALPHA BETA I02 TAU32 TAU42 GAMMA DELTA IW1 TAU11 TAU21 NW1 IW2 TAU12 TAU22 NW2 I0B T0B] = current_wave(current_type)

if strcmp(current_type,'farhad_first_return_stroke')

    % first double-exponential function

    I01 = 0.0e3;

    TAU31 = 80.0e-6;

    TAU41 = 6e-6;

    ALPHA = 1/TAU31;

    BETA = 1/TAU41;

    % second double-exponential function

    I02 = 0.0e3;

    TAU32 = 100.0e-6;

    TAU42 = 6e-6;

    GAMMA = 1/TAU32;

    DELTA = 1/TAU42;

    % first Heidler function
    
    IW1 = 28.0e3;
    
    TAU11 = 1.8e-6;
    
    TAU21 = 95.0e-6;
    
    NW1 = 2;

    % second Heidler function

    IW2 = 0.0e0;

    TAU12 = 5.0e-6;

    TAU22 = 50.e-6;

    NW2 = 2;

    % Barbosa function
    
    I0B = 0;
    
    T0B = 0.12e-6;

elseif strcmp(current_type,'farhad_subsequent_return_stroke')

    % first double-exponential function

    I01 = 0.0e3;

    TAU31 = 80.0e-6;

    TAU41 = 6e-6;

    ALPHA = 1/TAU31;

    BETA = 1/TAU41;

    % second double-exponential function

    I02 = 0.0e3;

    TAU32 = 100.0e-6;

    TAU42 = 6e-6;

    GAMMA = 1/TAU32;

    DELTA = 1/TAU42;

    % first Heidler function

    IW1 = 10.7e3;

    TAU11 = 0.25e-6;

    TAU21 = 2.5e-6;

    NW1 = 2;

    % second Heidler function

    IW2 = 6.5e3;

    TAU12 = 2.1e-6;

    TAU22 = 230.e-6;

    NW2 = 2;

    % Barbosa function
    
    I0B = 0;
    
    T0B = 0.12e-6;

elseif strcmp(current_type,'farhad_subsequent_return_stroke_decay')

    % first double-exponential function

    I01 = 0.0e3;

    TAU31 = 80.0e-6;

    TAU41 = 6e-6;

    ALPHA = 1/TAU31;

    BETA = 1/TAU41;

    % second double-exponential function

    I02 = 0.0e0;

    TAU32 = 100.0e-6;

    TAU42 = 6e-6;

    GAMMA = 1/TAU32;

    DELTA = 1/TAU42;

    % first Heidler function

    IW1 = 10.7e3;

    TAU11 = 0.25e-6;

    TAU21 = 2.5e-6;

    NW1 = 2;

    % second Heidler function

    IW2 = 6.5e3;

    TAU12 = 2.1e-6;

    TAU22 = 180.e-6;

    NW2 = 2;

    % Barbosa function
    
    I0B = 0;
    
    T0B = 0.12e-6;

elseif strcmp(current_type,'du_first_return_stroke')

    % first double-exponential function

    I01 = 0.0e0;

    TAU31 = 100.0e-6;

    TAU41 = 6e-6;

    ALPHA = 1/TAU31;

    BETA = 1/TAU41;

    % second double-exponential function

    I02 = 0.0e0;

    TAU32 = 100.0e-6;

    TAU42 = 6e-6;

    GAMMA = 1/TAU32;

    DELTA = 1/TAU42;

    % first Heidler function

    IW1 = 28e3;

    TAU11 = 0.3e-6;

    TAU21 = 6.0e-6;

    NW1 = 2;

    % second Heidler function

    IW2 = 16.0e3;

    TAU12 = 10.0e-6;

    TAU22 = 50.e-6;

    NW2 = 2;

    % Barbosa function
    
    I0B = 0;
    
    T0B = 0.12e-6;

elseif strcmp(current_type,'du_subsequent_return_stroke')

    % first double-exponential function

    I01 = 0.0e0;

    TAU31 = 100.0e-6;

    TAU41 = 6e-6;

    ALPHA = 1/TAU31;

    BETA = 1/TAU41;

    % second double-exponential function

    I02 = 0.0e0;

    TAU32 = 100.0e-6;

    TAU42 = 6e-6;

    GAMMA = 1/TAU32;

    DELTA = 1/TAU42;

    % first Heidler function

    IW1 = 13e3;

    TAU11 = 0.15e-6;

    TAU21 = 3.0e-6;

    NW1 = 2;

    % second Heidler function

    IW2 = 7.0e3;

    TAU12 = 5.0e-6;

    TAU22 = 50.e-6;

    NW2 = 2;

    % Barbosa function
    
    I0B = 0;
    
    T0B = 0.12e-6;

elseif strcmp(current_type,'carlo_alberto_subsequent_return_stroke')

    % first double-exponential function

    I01 = 7.5e3;

    TAU31 = 100.0e-6;

    TAU41 = 6e-6;

    ALPHA = 1/TAU31;

    BETA = 1/TAU41;

    % second double-exponential function

    I02 = 0.0e0;

    TAU32 = 100.0e-6;

    TAU42 = 6e-6;

    GAMMA = 1/TAU32;

    DELTA = 1/TAU42;

    % first Heidler function

    IW1 = 9.9e3;

    TAU11 = 0.072e-6;

    TAU21 = 5.0e-6;

    NW1 = 2;

    % second Heidler function

    IW2 = 0.0e3;

    TAU12 = 14.0e-6;

    TAU22 = 95.e-6;

    NW2 = 2;

    % Barbosa function
    
    I0B = 0;
    
    T0B = 0.12e-6;

elseif strcmp(current_type,'fredholm_subsequent_return_stroke')

    % first double-exponential function

    I01 = 0.0e0;

    TAU31 = 100.0e-6;

    TAU41 = 6e-6;

    ALPHA = 1/TAU31;

    BETA = 1/TAU41;

    % second double-exponential function

    I02 = 0.0e0;

    TAU32 = 100.0e-6;

    TAU42 = 6e-6;

    GAMMA = 1/TAU32;

    DELTA = 1/TAU42;

    % first Heidler function

    IW1 = 8.5e3;

    TAU11 = 0.12e-6;

    TAU21 = 14.0e-6;

    NW1 = 2;

    % second Heidler function

    IW2 = 3.2e3;

    TAU12 = 14.0e-6;

    TAU22 = 95.e-6;

    NW2 = 2;

    % Barbosa function
    
    I0B = 0;
    
    T0B = 0.12e-6;

elseif strcmp(current_type,'marcos_horizontal_electric_field')

    % first double-exponential function

    I01 = 10.0e3;

    TAU31 = 1/3e4;

    TAU41 = 1e-7;

    ALPHA = 1/TAU31;

    BETA = 1/TAU41;

    % second double-exponential function

    I02 = 0.0e0;

    TAU32 = 100.0e-6;

    TAU42 = 6e-6;

    GAMMA = 1/TAU32;

    DELTA = 1/TAU42;

    % first Heidler function

    IW1 = 0.0e0;

    TAU11 = 0.12e-6;

    TAU21 = 14.0e-6;

    NW1 = 2;

    % second Heidler function

    IW2 = 0.0e0;

    TAU12 = 14.0e-6;

    TAU22 = 95.e-6;

    NW2 = 2;

    % Barbosa function
    
    I0B = 0;
    
    T0B = 0.12e-6;

elseif strcmp(current_type,'barbosa_linear_ramp_current_original')

    % first double-exponential function

    I01 = 0.0e3;

    TAU31 = 1/80;

    TAU41 = 1/0.6e7;

    ALPHA = 1/TAU31;

    BETA = 1/TAU41;

    % second double-exponential function

    I02 = 0.0e3;

    TAU32 = 100.0e-6;

    TAU42 = 6e-6;

    GAMMA = 1/TAU32;

    DELTA = 1/TAU42;

    % first Heidler function

    IW1 = 0.0e3;

    TAU11 = 0.5e-6;

    TAU21 = 100000.0e-6;

    NW1 = 2;

    % second Heidler function

    IW2 = 0.0e3;

    TAU12 = 1e-6;

    TAU22 = 60000e-6;

    NW2 = 2;

    % Barbosa function
    
    I0B = 12.7e3;
    
    T0B = 0.12e-6;

elseif strcmp(current_type,'abbas')
    
    % first double-exponential function
    
    I01 = 0.0;
    
    TAU31 = 0.0;
    
    TAU41 = 0.0;
    
    ALPHA = 1.0;
    
    BETA = 1.0;
    
    % second double-exponential function
    
    I02 = 0.0;
    
    TAU32 = 0.0;
    
    TAU42 = 0.0;
    
    GAMMA = 1.0;
    
    DELTA = 1.0;
    
    % first Heidler function
    
    IW1 = 9.5e3;
    
    TAU11 = 0.5e-6;
    
    TAU21 = 63.0e-6;
    
    NW1 = 2;
    
    % second Heidler function
    
    IW2 = 0e3;
    
    TAU12 = 5.5e-6;
    
    TAU22 = 180.0e-6;
    
    NW2 = 5;
    
    % Barbosa function
    
    I0B = 0.0;
    
    T0B = 0.0;
    
end