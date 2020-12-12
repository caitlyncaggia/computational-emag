% ECE 6380 Homework 1.1
% Caitlyn Caggia

cellNumberList =[10 20 40 80 160 320];
wavelength = 1.0;
E = cell(1,length(cellNumberList));

for i = 1:length(cellNumberList)
    % Create input files
    inputFileName = ['inputfil' num2str(cellNumberList(i)) '.txt'];
    fh = fopen(inputFileName, 'wt');
    
    % Number of nodes
    fprintf(fh, [num2str(cellNumberList(i)+1) '\n']);
    
    % Node positions
    nodes = linspace(0,wavelength,cellNumberList(i)+1);
    for j = 1:cellNumberList(i)+1
        fprintf(fh, [num2str(j) ' ' num2str(nodes(j)) '\n']);
    end
    
    % Relative permittivity
    for k = 1:cellNumberList(i)+1
        fprintf(fh, [num2str(k) ' 1.0\n']);
    end
    
    fclose(fh);
    
    % Create output files
    newE = HW1_1fem(inputFileName);
    E(1,i) = {newE};
    
end


% Plot magnitude of E
figure

subplot(3,2,1)
plot(abs(E{1}))
title('10 cells')
axis([0 cellNumberList(1) 0.98 1.02])

subplot(3,2,2)
plot(abs(E{2}))
title('20 cells')
axis([0 cellNumberList(2) 0.98 1.02])

subplot(3,2,3)
plot(abs(E{3}))
title('40 cells')
axis([0 cellNumberList(3) 0.99 1.01])

subplot(3,2,4)
plot(abs(E{4}))
title('50 cells')
axis([0 cellNumberList(4) 0.99 1.01])

subplot(3,2,5)
plot(abs(E{5}))
title('160 cells')
axis([0 cellNumberList(5) 0.995 1.005])

subplot(3,2,6)
plot(abs(E{6}))
title('320 cells')
axis([0 cellNumberList(6) 0.995 1.005])
