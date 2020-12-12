function [TriFEMQ] = TriFEMQ

% compute potential and C/epsilon0 for 2D region modeled with curved
% triangular cells (quadratic)
%
% October 29, 2018   A. F. Peterson
 
% read fem mesh from file 'cylfil.txt'
% 
% mesh should be organized so that interior nodes appear first, followed
% by nodes on the outer boundary (the driven boundary), followed by nodes
% on the inner boundary

  global pcetond xy
  global icell
  
  nnodes = dlmread('cylfil.txt','', [0,0,0,0]);
  ncells = dlmread('cylfil.txt','', [0,1,0,1]);
  ninner = dlmread('cylfil.txt','', [0,2,0,2]);
  nouter = dlmread('cylfil.txt','', [0,3,0,3]);
 
  xy=dlmread('cylfil.txt','', [1,0,nnodes,1]);

  nstart=nnodes + 1;
  nend=nstart + ncells - 1;
  pcetond=dlmread('cylfil.txt','', [nstart,0,nend,5]);

  
  nstart=nend + 1;
  nend=nstart + ninner - 1;
%  pinner=dlmread('cylfil.txt','', [nstart,0,nend,0]);

  nstart=nend + 1;
  nend=nstart + nouter - 1;
%  pouter=dlmread('cylfil.txt','', [nstart,0,nend,0]);

  nstart=nend + 1;
  nend=nstart + ncells - 1;
  er=dlmread('cylfil.txt','', [nstart,0,nend,0]);
  
% initialize variables

  rerr = 1.0e-8;
  aerr = 1.0e-12;

  nunks = nnodes - ninner - nouter;
  
  Wtilda=zeros(nnodes);
  W=zeros(nunks);
  V=zeros(nunks,1);
  
% loop through the cells, filling global matrix one cell at a time

  for icell=1:ncells

%    compute 6 by 6 element matrix for cell 'icell'

     eleS = elemat(rerr,aerr);
                 
%    add contributions from cell 'icell' to global matrices

     for ii=1:6
        ig=pcetond(icell,ii); % 'ig' is the global node for 'ii'
        for jj=1:6
           jg=pcetond(icell,jj);  % 'jg' is the global node for 'jj'
 		   Wtilda(ig,jg) = Wtilda(ig,jg) + er(icell)*eleS(ii,jj);
           if(ig <= nunks) % test function at interior node
               if(jg <= nunks) % basis function at interior node
 			      W(ig,jg) = W(ig,jg) + er(icell)*eleS(ii,jj);
               elseif(jg <= nunks+nouter) % basis function on outer bnd
                  V(ig) = V(ig) - er(icell)*eleS(ii,jj);
               end
           end
        end
     end
  end

%  solve the system of equations to find the potential function
  
  Pot = W\V;

%  add potential functions on the boundaries to the list

  Pot(nunks+nouter+ninner)=0;
  nstart=nunks+1;
  nend=nunks+nouter;
  for ii=nstart:nend
      Pot(ii)=1;
  end
  nstart=nend+1;
  nend=nunks+nouter+ninner;
  for ii=nstart:nend
      Pot(ii)=0;
  end
  
%  compute the capacitance

  Cap = Pot.'*Wtilda*Pot;

%  write results to file 'potfil.txt'

   str = ['Capacitance/epsilon0 pul = ',num2str(Cap)];  disp(str);
  
   fid = fopen('potfil.txt', 'wt');

   fprintf(fid,'%s\n\n',str);
   str = 'node    Potential';
   fprintf(fid,'%s\n\n',str);
   for ii=1:nend
      fprintf(fid,'%6d %15.14g\n',ii, Pot(ii));
   end
     
   TriFEMQ = Cap;
end

% ----------------------------------------------------------------------

function [eleS,eleT] = elemat(rerr,aerr)

%  elemat: construct the 6 by 6 element matrices for the contributions of
%          basis and testing functions of the form
%                   grad T dot grad B    and    T times B
%          over a triangular cell
%
%          done by quadrature over integrands 'Sint' and 'Tint'

     global itest ibasis;

     eleS(6,6)=0;  eleT(6,6)=0;

%    compute 6 by 6 element matrices S, T

     for ii=1:6
        itest=ii;
        for jj=1:6
           ibasis=jj;
           eleS(ii,jj)=rtriad('Sint',rerr,aerr);
%          eleT(ii,jj)=rtriad('Tint',rerr,aerr);  (T not used in this code)

        end
     end
end
         
