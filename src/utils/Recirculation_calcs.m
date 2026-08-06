base_recirc = readmatrix("recirculation_base")
base_recirc = sortrows(base_recirc,[1]);
control_recirc = readmatrix("recirculation_control")
control_recirc = sortrows(control_recirc,[1]);
half_recirc = readmatrix("recirculation_half")
half_recirc = sortrows(half_recirc,[1]);


base = 1.04-0.5; %recirculation length from cfd
control = 2.12-0.5
min = 4.65-0.5;