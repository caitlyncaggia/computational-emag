function fdifftot1

% use finite difference formulas to solve the scalar Helmholtz
% equation for EM fields in a dielectric slab
%
% backwards derivative ABC
%
% June 28, 2018   A. F. Peterson

 
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
  phs=1j*k0*delta;
  Z=zeros(n_unknowns);
  RHS=zeros(n_unknowns,1);

  for irow=1:n_unknowns

      if (irow == 1)
   Z(irow,irow) = -1/delta^2  - 1j*k0/delta;
   Z(irow,irow+1) = 1/delta^2;
          
      elseif(irow == n_unknowns)
   Z(irow,irow-1) =  1/delta^2;
   Z(irow,irow) = -1/delta^2  - 1j*k0/delta;
          
      else
   Z(irow,irow-1)=1./delta^2;
   Z(irow,irow)=-2./delta^2 + k0^2 * epsilon(irow);
   Z(irow,irow+1)=1./delta^2;
      end
 end
 
% disp(Z);
 
% fill excitation vector (right hand side)

 RHS(1) = (exp(-phs)-1)/delta^2 -1j*k0/delta;  % assumes that incident Ey(a)=1
 RHS(n_unknowns) = ((exp(phs)-1)/delta^2 -1j*k0/delta) * ...
     exp(-1j*k0*x(n_nodes)); 

 disp(delta); disp(RHS(n_unknowns));
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
 fprintf(fid,'%44s\n\n',str);

 for irow=1:n_unknowns
         mag=abs(E(irow));
         phs=180*atan2(imag(E(irow)),real(E(irow)))/pi;
         fprintf(fid,'%6d %15.14g %15.14g\n',irow, mag, phs);

 end
end