%% ECE 6380 - Homework 7
% Caitlyn Caggia
% Problem 4

%% Table 5
% before and after slab space = 0.1, slab length = 0.2, epsilon = 4

% exact solutions
Rexact5 = [0.40338 60.245];
R5 = Rexact5(1)*exp(Rexact5(2)*j);
Texact5 = [0.91503 150.245];
T5 = Texact5(1)*exp(Texact5(2)*j);

% 8 cells
% 2 left, 4 in slab, 2 right
R8raw = [0.40356 60.2634]; % reflection coefficient
R8 = R8raw(1)*exp(R8raw(2)*j);
T8raw = [0.91495 150.2634]; % transmission coefficient
T8 = T8raw(1)*exp(T8raw(2)*j);
R8error = abs((R8 - R5))/abs(R5) * 100;
T8error = abs((T8 - T5))/abs(T5) * 100;

% 16 cells
% 5 left, 6 in slab, 5 right
R16raw = [0.40342 60.2486];
R16 = R16raw(1)*exp(R16raw(2)*j);
T16raw = [0.91502 150.2486];
T16 = T16raw(1)*exp(T16raw(2)*j);
R16error = abs((R16 - R5))/abs(R5) * 100;
T16error = abs((T16 - T5))/abs(T5) * 100;

% 32 cells
% 10 left, 12 in slab, 10 right
R32raw = [0.40338 60.2453];
R32 = R32raw(1)*exp(R32raw(2)*j);
T32raw = [0.91503 150.2453];
T32 = T32raw(1)*exp(T32raw(2)*j);
R32error = abs((R32 - R5))/abs(R5) * 100;
T32error = abs((T32 - T5))/abs(T5) * 100;

% 64 cells
% 21 left, 22 in slab, 21 right
R64raw = [0.40338 60.2451];
R64 = R64raw(1)*exp(R64raw(2)*j);
T64raw = [0.91503 150.2451];
T64 = T64raw(1)*exp(T64raw(2)*j);
R64error = abs((R64 - R5))/abs(R5) * 100;
T64error = abs((T64 - T5))/abs(T5) * 100;

fprintf('Table 5 Reflection')
ratio16to8R = R8error/R16error;
ratio32to16R = R16error/R32error;
ratio64to32R = R32error/R64error;
pR1 = log(ratio16to8R)/log(2)
pR2 = log(ratio32to16R)/log(2)
pR3 = log(ratio64to32R)/log(2)

fprintf('Table 5 Transmission')
ratio16to8T = T8error/T16error;
ratio32to16T = T16error/T32error;
ratio64to32T = T32error/T64error;
pT1 = log(ratio16to8T)/log(2)
pT2 = log(ratio32to16T)/log(2)
pT3 = log(ratio64to32T)/log(2)

%% Table 6
% before and after slab space = 0.05, slab length = 0.4, epsilon = 5
Rexact6 = [0.4824168 100.35461];
R6 = Rexact6(1)*exp(Rexact6(2)*j);
Texact6 = [0.8759418 10.35461];
T6 = Texact6(1)*exp(Texact6(2)*j);

% 10 cells
% 3 left, 4 slab, 3 right
R10raw = [0.49068 101.2204];
R10 = R10raw(1)*exp(R10raw(2)*j);
T10raw = [0.87134 11.2204];
T10 = T10raw(1)*exp(T10raw(2)*j);
R10error = abs((R10 - R6)/R6) * 100;
T10error = abs((T10 - T6)/T6) * 100;

% 20 cells
% 6 left, 8 slab, 6 right
R20raw = [0.48295 100.4126];
R20 = R20raw(1)*exp(R20raw(2)*j);
T20raw = [0.87565 10.4126];
T20 = T20raw(1)*exp(T20raw(2)*j);
R20error = abs((R20 - R6)/R6) * 100;
T20error = abs((T20 - T6)/T6) * 100;

% 40 cells
% 13 left, 14 slab, 13 right
R40raw = [0.48247 100.3609];
R40 = R40raw(1)*exp(R40raw(2)*j);
T40raw = [0.87591 10.3609];
T40 = T40raw(1)*exp(T40raw(2)*j);
R40error = abs((R40 - R6)/R6) * 100;
T40error = abs((T40 - T6)/T6) * 100;

% 80 cells
% 26 left, 28 slab, 26 right
R80raw = [0.48242 100.355];
R80 = R80raw(1)*exp(R80raw(2)*j);
T80raw = [0.87594 10.355];
T80 = T80raw(1)*exp(T80raw(2)*j);
R80error = abs((R80 - R6)/R6) * 100;
T80error = abs((T80 - T6)/T6) * 100;

fprintf('Table 6 Reflection')
ratio10to20R = R10error/R20error;
ratio20to40R = R20error/R40error;
ratio40to80R = R40error/R80error;
pR4 = log(ratio10to20R)/log(2)
pR5 = log(ratio20to40R)/log(2)
pR6 = log(ratio40to80R)/log(2)

fprintf('Table 6 Transmission')
ratio10to20T = T10error/T20error;
ratio20to40T = T20error/T40error;
ratio40to80T = T40error/T80error;
pT4 = log(ratio10to20T)/log(2)
pT5 =log(ratio20to40T)/log(2)
pT6 = log(ratio40to80T)/log(2)

