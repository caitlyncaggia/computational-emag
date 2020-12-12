function [E] = HW1_1fem(inputFileName)

% use finite element formulas to solve the scalar Helmholtz
% equation for EM fields in a dielectric slab
%
% June 21, 2018   A. F. Peterson
% Modified August 28, 2018 by Caitlyn Caggia

 
% read mesh from file 'inputfil.txt'

 n_nodes = dlmread(inputFileName,'', [0,0,0,0]);

 x=dlmread(inputFileName,'', [1,1,n_nodes,1]);
 
 nstart = n_nodes + 1;
 nend = nstart + n_nodes - 1;
 epsilon = dlmread(inputFileName,'', [nstart,1,nend,1]);

% initialize variables

  k0 = 2*pi;
  n_unknowns = n_nodes;
  delta = x(2) - x(1);
  Z=zeros(n_unknowns);
  RHS=zeros(n_unknowns,1);

% specify Dirichlet boundary conditions

 Ka = 1;
 %Kb = exp(-1j*k0*x(n_nodes));
 
% fill global matrix

 Z(1,1)=1;
 Z(n_unknowns,n_unknowns)=-1./delta^2 + k0^2 * ...
     epsilon(n_unknowns-1)/3 + j*k0/delta;
 Z(n_unknowns, n_unknowns-1) = 1./delta^2 + k0^2 * epsilon(n_unknowns-1)/6;
 
 for irow=2:n_unknowns-1

   Z(irow,irow-1)=1./delta^2 + k0^2 * epsilon(irow-1)/6;
   Z(irow,irow)=-2./delta^2 + k0^2 * (epsilon(irow-1) + epsilon(irow))/3;
   Z(irow,irow+1)=1./delta^2 + k0^2 * epsilon(irow)/6;

 end

 disp(Z);
 
% fill excitation vector (right hand side)

 RHS(1)=Ka;
 RHS(n_unknowns)=0;

% solve system of equations

 E = Z\RHS;
 disp(E);

% write fields to output file
 outputFileName = ['outputfil' num2str(n_nodes-1) '.txt'];
 fid = fopen(outputFileName, 'wt');

 for irow=1:n_unknowns
         mag=abs(E(irow));
         phs=180*atan2(imag(E(irow)),real(E(irow)))/pi;
         fprintf(fid,'%6d %15.14g %15.14g\n',irow, mag, phs);
 end
 
 fclose(fid);
end