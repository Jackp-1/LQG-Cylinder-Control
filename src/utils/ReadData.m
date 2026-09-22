% Function Read Data
% Reads data from datpath and stores in workspace

Train_lift = readmatrix('Train_lift.txt');
Train_drag = readmatrix('Train_drag.txt');
Train_input = readmatrix( 'Train_input.txt' );
Base_drag = readmatrix( 'Base_drag.txt' );
Base_lift = readmatrix( 'Base_lift.txt' );
Half_drag = readmatrix( 'Half_cyl_drag.txt' );
Mode_4 = readmatrix( '4Mode.txt' );
Mode_5 = readmatrix( '5Mode.txt' );
Mode_6 = readmatrix( '6Mode.txt' );
Mode_7 = readmatrix( '7Mode.txt' );
Mode_8 = readmatrix( '8Mode.txt' );
Mode_9 = readmatrix( '9Mode.txt' );
Mode_10 = readmatrix( '10Mode.txt' );
Mode_11 = readmatrix( '11Mode.txt' );
Mode_12 = readmatrix( '12Mode.txt' );
Mode_13 = readmatrix( '13Mode.txt' );
Mode_14 = readmatrix( '14Mode.txt' );
Mode_15 = readmatrix( '15Mode.txt' );
Mode_16 = readmatrix( '16Mode.txt' );
Mode_18 = readmatrix( '18Mode.txt' );
Half_Cyl_Drag = readmatrix( 'Half_cyl_drag.txt' );
Sin_Force = readmatrix( 'Sin_forcing.txt' );
Lift_constant_double = readmatrix( 'lift_constant_double.txt' );
Lift_constant_average = readmatrix( 'lift_constant_average.txt' );
Lift_constant_half = readmatrix( 'lift_constant_half.txt' );
Lift_constant_quarter = readmatrix( 'lift_constant_quarter.txt' );
Drag_constant_double = readmatrix( 'drag_constant_double.txt' );
Drag_constant_half = readmatrix( 'drag_constant_half.txt' );
Drag_constant_quarter = readmatrix( 'drag_constant_quarter.txt' );
Drag_constant_average = readmatrix( 'drag_constant_average.txt' );
R2e_5 = readmatrix( 'R2e-5.txt' );
R1e_5 = readmatrix( 'R1e-5.txt' );
R75e_7 = readmatrix( 'R75e-7.txt' );
R5e_6 = readmatrix( 'R5e-6.txt' );
R5e_5 = readmatrix( 'R5e-5.txt' );

R75e_7_noise = readmatrix( 'R75e-7_noise.txt' );
R75e_7_noise2 = readmatrix( 'R75e-7_noise2.txt' );
Noise = readmatrix( 'Noise.txt' );
Disturbance_Noise = readmatrix( 'Disturbance_noise_rejection.txt' );
Mode_11_Integral = readmatrix( '11Mode_integral.txt' );
Mode_11_Integral2 = readmatrix( '11Mode_Int2.txt' );
Ref_Tracking = readmatrix( 'Ref_tracking.txt' );

SRS_Base_Drag = readmatrix('SRS_Base_drag.txt');
SRS_Base_Lift = readmatrix('SRS_Base_lift.txt');
SRS_Control = readmatrix('SRS_control_2.txt');