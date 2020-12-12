% ECE 6380 Homework 1.2
% Caitlyn Caggia

Rgamma = [  4 0.4343 68.27;
            8 0.4156 62.59;
            16 0.4068 60.86;
            32 0.4043 60.40;
            64 0.4036 60.28];
Rgexact = [0.40338 60.245];
Rgest = cell(size(Rgamma,1)-1, size(Rgamma,2)+1);

Rtau = [    4 0.9008 158.27;
            8 0.9096 152.59;
            16 0.9135 150.86;
            32 0.9147 150.40;
            64 0.9149 150.28];
Rtexact = [0.91503 150.245];
Rtest = cell(size(Rtau,1)-1, size(Rtau,2)+1);

for i = 1:size(Rgest,1)
    
    % Start new row
    Rgest{i,1} = [num2str(Rgamma(i,1)) '/' num2str(Rgamma(i+1,1))];
    Rtest{i,1} = [num2str(Rtau(i,1)) '/' num2str(Rtau(i+1,1))];
    
    % Extrapolate magnitude of gamma and tau
    Rgest{i,2} = (4*Rgamma(i+1,2) - Rgamma(i,2)) / 3;
    Rtest{i,2} = (4*Rtau(i+1,2) - Rtau(i,2)) / 3;
        
    % Extrapolate the angle of gamma and tau
    Rgest{i,3} = (4*Rgamma(i+1,3) - Rgamma(i,3)) / 3;
    Rtest{i,3} = (4*Rtau(i+1,3) - Rtau(i,3)) / 3;
    
    % Compute error for magnitude of gamma and tau
    Rgest{i,4} = abs(Rgexact(1)-Rgest{i,2})/Rgexact(1) * 100;
    Rtest{i,4} = abs(Rtexact(1)-Rtest{i,2})/Rtexact(1) * 100;
    
end

% Format table for extrapolated gamma values
Rgheaders = {'Number of Cells' '|gamma|' '<gamma' 'Error'};
Rg_table = [Rgheaders; Rgest];
Rg_table = [Rg_table; {'Exact' Rgexact(1) Rgexact(2) []}]

% Format table for extrapolated tau values
Rtheaders = {'Number of Cells' '|tau|' '<tau' 'Error'};
Rt_table = [Rtheaders; Rtest];
Rt_table = [Rt_table; {'Exact' Rtexact(1) Rtexact(2) []}]
