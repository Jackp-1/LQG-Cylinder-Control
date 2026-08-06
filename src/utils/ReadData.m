% Function Read Data
% Reads data from datpath and stores in workspace

Train_lift = readmatrix([datpath,'Train_lift.txt']);
Train_drag = readmatrix([datpath,'Train_drag.txt']);
Train_input = readmatrix([datpath,'Train_input.txt']);
Base_drag = readmatrix([datpath, 'Base_drag.txt']);
Base_lift = readmatrix([datpath, 'Base_lift.txt']);
Half_drag = readmatrix([datpath, 'Half_cyl_drag.txt']);
Mode_4 = readmatrix([datpath, '4Mode.txt']);
Mode_5 = readmatrix([datpath, '5Mode.txt']);
Mode_6 = readmatrix([datpath, '6Mode.txt']);
Mode_7 = readmatrix([datpath, '7Mode.txt']);
Mode_8 = readmatrix([datpath, '8Mode.txt']);
Mode_9 = readmatrix([datpath, '9Mode.txt']);
Mode_10 = readmatrix([datpath, '10Mode.txt']);
Mode_11 = readmatrix([datpath, '11Mode.txt']);
Mode_12 = readmatrix([datpath, '12Mode.txt']);
Mode_13 = readmatrix([datpath, '13Mode.txt']);
Mode_14 = readmatrix([datpath, '14Mode.txt']);
Mode_15 = readmatrix([datpath, '15Mode.txt']);
Mode_16 = readmatrix([datpath, '16Mode.txt']);
Mode_18 = readmatrix([datpath, '18Mode.txt']);
Half_Cyl_Drag = readmatrix([datpath, 'Half_cyl_drag.txt']);
Sin_Force = readmatrix([datpath, 'Sin_forcing.txt']);
Lift_constant_double = readmatrix([datpath, 'lift_constant_double.txt']);
Lift_constant_average = readmatrix([datpath, 'lift_constant_average.txt']);
Lift_constant_half = readmatrix([datpath, 'lift_constant_half.txt']);
Lift_constant_quarter = readmatrix([datpath, 'lift_constant_quarter.txt']);
Drag_constant_double = readmatrix([datpath, 'drag_constant_double.txt']);
Drag_constant_half = readmatrix([datpath, 'drag_constant_half.txt']);
Drag_constant_quarter = readmatrix([datpath, 'drag_constant_quarter.txt']);
Drag_constant_average = readmatrix([datpath, 'drag_constant_average.txt']);
R2e_5 = readmatrix([datpath, 'R2e-5.txt']);
R1e_5 = readmatrix([datpath, 'R1e-5.txt']);
R75e_7 = readmatrix([datpath, 'R75e-7.txt']);
R5e_6 = readmatrix([datpath, 'R5e-6.txt']);
R5e_5 = readmatrix([datpath, 'R5e-5.txt']);

R75e_7_noise = readmatrix([datpath, 'R75e-7_noise.txt']);
R75e_7_noise2 = readmatrix([datpath, 'R75e-7_noise2.txt']);
Noise = readmatrix([datpath, 'Noise.txt']);
Disturbance_Noise = readmatrix([datpath, 'Disturbance_noise_rejection.txt']);
Mode_11_Integral = readmatrix([datpath, '11Mode_integral.txt']);
Mode_11_Integral2 = readmatrix([datpath, '11Mode_Int2.txt']);
Ref_Tracking = readmatrix([datpath, 'Ref_tracking.txt']);

SRS_Base_Drag = readmatrix([datpath, 'SRS_Base_drag.txt']);
SRS_Base_Lift = readmatrix([datpath, 'SRS_Base_lift.txt']);
SRS_Control = readmatrix([datpath, 'SRS_control.txt']);