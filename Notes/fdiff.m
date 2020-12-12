function fdiff

% use finite difference formulas to solve the scalar Helmholtz
% equation for EM fields in a dielectric slab
%
% June 21, 2018   A. F. Peterson

 
% read mesh from file 'inputfil.txt'

 n_nodes = dlmread('inputfil.txt','', [0,0,0,0]);

 x=dlmread('inputfil.txt','', [1,1,n_nodes,1]);
 
 nstart = n_nodes + 1;
 nend = nstart + n_nodes - 1;
 epsilon = dlmread('inputfil.txt','', [nstart,1,nend,1]);

% initialize variables

  k0 = 2*pi;
  n_unknowns = n_nodes;
  delta = x(2) - x(1);
  Z=zeros(n_unknowns);
  RHS=zeros(n_unknowns,1);

% specify Dirichlet boundary conditions

 Ka = 1;
 Kb = exp(-1j*k0*x(n_nodes));
 
% fill global matrix

 Z(1,1)=1;
 Z(n_unknowns,n_unknowns)=1;
 
 for irow=2:n_unknowns-1;

   Z(irow,irow-1)=1./delta^2;
   Z(irow,irow)=-2./delta^2 + k0^2 * epsilon(irow);
   Z(irow,irow+1)=1./delta^2;

 end

 disp(Z);
 
% fill excitation vector (right hand side)

 RHS(1)=Ka;
 RHS(n_unknowns)=Kb;

% solve system of equations

 E = Z\RHS;

% write fields to output file
 
 fid = fopen('outputfil.txt', 'wt');

 for irow=1:n_unknowns;
         mag=abs(E(irow));
         phs=180*atan2(imag(E(irow)),real(E(irow)))/pi;
         fprintf(fid,'%6d %15.14g %15.14g\n',irow, mag, phs);

 end
end