function HW2_1femtot(inputFile)

% use finite element formulas to solve the scalar Helmholtz
% equation for EM fields in a dielectric slab
%
% June 25, 2018   A. F. Peterson
% Modified September 5, 2018 by Caitlyn Caggia


% read mesh from file 'inputfil.txt'

n_nodes = dlmread(inputFile,'', [0,0,0,0]);

x=dlmread(inputFile,'', [1,1,n_nodes,1]);

nstart = n_nodes + 1;
nend = nstart + n_nodes - 2;
epsilon = dlmread(inputFile,'', [nstart,1,nend,1]);

% initialize variables

k0 = 2*pi;
n_unknowns = n_nodes;
Z=zeros(n_unknowns);
RHS=zeros(n_unknowns,1);

% fill global matrix

for irow=1:n_unknowns
    
    if (irow == 1)
        deltaR = x(2) - x(1);
        Z(irow,irow) = 1/deltaR - k0^2*epsilon(irow)*deltaR/3 + 1j*k0;
        Z(irow,irow+1) = -1/deltaR - k0^2*epsilon(irow)*deltaR/6;
        
    elseif(irow == n_unknowns)
        deltaL = x(irow) - x(irow-1);
        Z(irow,irow-1) = -1/deltaL - k0^2*epsilon(irow-1)*deltaL/6;
        Z(irow,irow) = 1/deltaL - k0^2*epsilon(irow-1)*deltaL/3 + 1j*k0;
        
    else
        deltaR = x(irow+1) - x(irow);
        deltaL = x(irow) - x(irow-1);
        Z(irow,irow-1) = -1/deltaL - k0^2*epsilon(irow-1)*deltaL/6;
        Z(irow,irow) = 1/deltaL + 1/deltaR...
            - k0^2*(epsilon(irow-1)*deltaL/3 + epsilon(irow)*deltaR/3);
        Z(irow,irow+1) = -1/deltaR - k0^2*epsilon(irow)*deltaR/6;
    end
end

% disp(Z);

% fill excitation vector (right hand side)

RHS(1) = 2*1j*k0*1;  % assumes that incident Ey(a)=1

% solve system of equations

E = Z\RHS;

% write fields to output file

fid = fopen('outputfil.txt', 'wt');

gamma = E(1)-1;
mag=abs(gamma);
phs=180*atan2(imag(gamma),real(gamma))/pi;
str = ['reflection coeff = ',num2str(mag),' ',num2str(phs)];
disp(str);
fprintf(fid,'%50s\n\n',str);
tau = E(n_unknowns);
mag=abs(tau);
phs=180*atan2(imag(tau),real(tau))/pi;
str = ['transmission coeff = ',num2str(mag),' ',num2str(phs)];
disp(str);
fprintf('\n\n');
fprintf(fid,'%44s\n\n',str);
for irow=1:n_unknowns
    mag=abs(E(irow));
    phs=180*atan2(imag(E(irow)),real(E(irow)))/pi;
    fprintf(fid,'%6d %15.14g %15.14g\n',irow, mag, phs);
    
end

fclose(fid);

end