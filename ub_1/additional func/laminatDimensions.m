% inputs

%left
ODln_left=175e-3; IDln_left=135e-3; Lln_left=22e-3; % lock nut
ODsens_left=175e-3; IDsens_left=143e-3; Lsens_left = 55.5*1e-3;% left side sensor + spacer lamination
ODamb_left=175e-3; IDamb_left=145e-3; Lamb_left=84.0*1e-3; % left side AMB lamination

% right 1
Lamb_right1=84.0*1e-3; ODamb_right1=175e-3;  IDamb_right1=145e-3;  %radial AMB
Lsens_right1 = 96*1e-3;ODsens_right1=175e-3; IDsens_right1=143e-3; %sensor lamination_right side 1

Lspacer_right1 = 40.5*1e-3;ODspacer_right1=175e-3; IDspacer_right1=143e-3; %sensor lamination_right side 1
% right 2
Lamb_right2=84.0*1e-3; ODamb_right2=175e-3; IDamb_right2=141e-3;  %radial AMB
Lsens_right2 = 55.5*1e-3;ODsens_right2=175e-3; IDsens_right2=139e-3; %sensor lamination_right side 1
ODln_right=175e-3; IDln_right=135e-3; Lln_right=26.5e-3;

%% left side
[Lln_left 0 IDln_left 1 ODln_left 2]; addSection;       %locknut left
[Lsens_left 0 IDsens_left 1 ODsens_left 2]; addSection; %left side sensor + spacer lamination
[Lspacer_left 0 IDsens_left 1 ODsens_left 2]; addSection; %left side sensor + spacer lamination
[Lamb_left 0 IDamb_left 1 ODamb_left 2]; addSection;    % left AMB

%% right side

[Lamb_right1    0 IDamb_right1      1 ODamb_right1      2] ; addSection; % right AMB 1
[Lspacer_right1 0 IDspacer_right1   1 ODspacer_right1   2]; addSection; %left side sensor + spacer lamination
[Lsens_right1   0 IDsens_right1     1 ODsens_right1     2] ; addSection; % right side sensor 1 + spacer lamination
[Lspacer_right2 0 IDspacer_right2   1 ODspacer_right2   2]; addSection; %left side sensor + spacer lamination
[Lamb_right2    0 IDamb_right2      1 ODamb_right2  2] ; addSection; % right AMB 2
[Lspacer_right3 0 IDspacer_right3   1 ODspacer_right3   2]; addSection; %left side sensor + spacer lamination
[Lsens_right2   0 IDsens_right2     1 ODsens_right2 2] ; addSection; % right side sensor 2 + spacer lamination
[Lln_right      0 IDln_right        1 ODln_left 2]; addSection;       %locknut left

%% EXTRA BEARING STUFF FOR EXTENDED SHAFT
[Lln_left 0 IDln_left 1 ODln_left 2]; addSection;       %locknut left
[Lamb_right3  0 IDamb_right3  1 ODamb_right3  2] ; addSection; % right AMB 3
[Lsens_right3 0 IDsens_right3 1 ODsens_right3 2] ; addSection; % right side sensor 3 + spacer lamination
[Lln_left 0 IDln_left 1 ODln_left 2]; addSection;       %locknut left