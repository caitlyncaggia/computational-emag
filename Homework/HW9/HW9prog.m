function HW9prog

% compute resonant wavenumbers for 2D cavity modeled with tri cells
% using curl-conforming vector basis functions 
%
%  computes TM and TE wavenumbers
%
% December 3 2018   A. F. Peterson

 global pcetond xy;
 
% read mesh from file 'cylfil.txt'

 nnodes = dlmread('cylfil.txt','', [0,0,0,0]);
 ncells = dlmread('cylfil.txt','', [0,1,0,1]);
 nedges = dlmread('cylfil.txt','', [0,2,0,2]);
 ninter = dlmread('cylfil.txt','', [0,3,0,3]);   %   interior nodes

 xy=dlmread('cylfil.txt','', [1,1,nnodes,2]);

 nstart=nnodes + 1;
 nend=nstart + ncells - 1;
 pcetond=dlmread('cylfil.txt','', [nstart,1,nend,3]);

 nstart=nend + 1;
 nend=nstart + ncells - 1;
 pcetoed=dlmread('cylfil.txt','', [nstart,1,nend,3]);

% nstart=nend + 1;
% nend=nstart + ncells - 1;
% epsilon=dlmread('cylfil.txt','', [nstart,1,nend,1]);

% nstart=nend + 1;
% nend=nstart + nedges - 1;
% pedtoce=dlmread('cylfil.txt','', [nstart,1,nend,2]);

% nstart=nend + 1;
% nend=nstart + nedges - 1;
% pedtond=dlmread('cylfil.txt','', [nstart,1,nend,2]);

  fid = fopen('eigfil.txt', 'wt');

  str = 'vector Helmholtz eqn CT/LN FEM analysis of cavity with: ';
  fprintf(fid,'%s \n',str);
  disp(str);

  str = ['number of cells = ',num2str(ncells)];
  fprintf(fid,'%s \n',str);
  disp(str);
  str = ['number of nodes = ',num2str(nnodes)];
  fprintf(fid,'%s \n',str);
  disp(str);
  str = ['number of edges = ',num2str(nedges)];
  fprintf(fid,'%s \n\n',str);
  disp(str);

% initialize variables

  nunks = nedges;
  
  A=zeros(nunks);
  B=zeros(nunks);
  if(ninter > 0)
     C(ninter,ninter)=0;
     D(ninter,ninter)=0;
  end

%  echo input variables
%  str='nodes';
%  disp(str);
%  disp(xy);
  
% fill global matrix one cell at a time

   for icell=1:ncells

     [eleS,eleT]=elemat(icell);
     
%    correct edge-based entries for global orientation; basis points
%    from smaller node index to larger node index

      n1=pcetond(icell,1);
      n2=pcetond(icell,2);
      n3=pcetond(icell,3);

     for kk=1:3
         iedg=0;
         if ((kk == 1) && (n2 > n3) )
             iedg=1;
         elseif ((kk == 2) && (n3 > n1) )
             iedg=1;
         elseif ((kk == 3) && (n1 > n2) )
             iedg=1;
         end
         if (iedg == 1)

%           flip signs of row and column kk

            for jj=1:3
                eleS(kk,jj)=-eleS(kk,jj);
                eleS(jj,kk)=-eleS(jj,kk);
                eleT(kk,jj)=-eleT(kk,jj);
                eleT(jj,kk)=-eleT(jj,kk);
            end
         end
     end
                 
%    add contributions from cell 'i' to global matrix

     for jj=1:3
         jp=pcetoed(icell,jj);             
         for kk=1:3
             kp=pcetoed(icell,kk);
 			 B(jp,kp)=B(jp,kp)+eleT(jj,kk);
 			 A(jp,kp)=A(jp,kp)+eleS(jj,kk);
         end
     end
   end

%  str='A matrix';
%  disp(str);
%  disp(A);

%  str='B matrix';
%  disp(str);
%  disp(B);
 
   if(ninter > 0)
 
%     duplicate system, truncate, and compute eigenvalues for TE case
 
      for jj=1:ninter
         for kk=1:ninter
           C(jj,kk)=A(jj,kk);
           D(jj,kk)=B(jj,kk);
         end
      end
   
      E = sort(eig(C,D));
       
      str = 'TE resonant wavenumbers (real / imag): ';
      fprintf(fid,'%s \n\n',str);

      for ii=1:ninter
         reaE=real(sqrt(E(ii)));
         imaE=imag(sqrt(E(ii)));
         fprintf(fid,'%6d %15.14g %15.14g\n',ii, reaE,imaE);
      end
   end
   
% compute eigenvalues for TM case

  E = sort(eig(A,B));
 
  fprintf(fid,'\n\n');
  str = 'TM resonant wavenumbers (real / imag): ';
  fprintf(fid,'%s \n\n',str);

 for ii=1:nunks
     reaE=real(sqrt(E(ii)));
     imaE=imag(sqrt(E(ii)));
   fprintf(fid,'%6d %15.14g %15.14g\n',ii, reaE,imaE);
 end
     
 disp('output placed in EIGFIL.TXT');
 
end

% ----------------------------------------------------------------------

function [eleS,eleT] = elemat(icell)

%
%  elemat: construct the element matrix for the contributions of
%          vector basis and testing functions of the form
%
%          curl T dot curl B  (ele2)
%                 and 
%               T dot B       (ele1)
%
%          over a (non-curved) quadrilateral cell.
%
     global pcetond xy;

     eleS(3,3)=0;  eleT(3,3)=0;
     n1=pcetond(icell,1);
     n2=pcetond(icell,2);
     n3=pcetond(icell,3);

	 x(1)=xy(n1,1);  y(1)=xy(n1,2);
	 x(2)=xy(n2,1);  y(2)=xy(n2,2);
	 x(3)=xy(n3,1);  y(3)=xy(n3,2);

     b(1)=y(2)-y(3);
     b(2)=y(3)-y(1);
     b(3)=y(1)-y(2);

     c(1)=x(3)-x(2);
     c(2)=x(1)-x(3);
     c(3)=x(2)-x(1);
     
     w(1) = sqrt((x(3)-x(2))^2 + (y(3)-y(2))^2);
     w(2) = sqrt((x(1)-x(3))^2 + (y(1)-y(3))^2); 
     w(3) = sqrt((x(2)-x(1))^2 + (y(2)-y(1))^2);
     
     Area = abs(b(3)*c(1) - b(1)*c(3))*0.5;

     for m=1:3
        mp1=circ(m+1);
        mp2=circ(m+2);
        
        for n=1:3
           np1=circ(n+1);
           np2=circ(n+2);
            
           eleS(m,n)=(w(m)*w(n)/((Area^3)*4))* ... 
           (b(mp1)*c(mp2)-b(mp2)*c(mp1))*(b(np1)*c(np2)-b(np2)*c(np1));
            
           eleTij = 0; 
           for ii=1:2
              for jj=1:2

              if ii == jj 
                 alpha = 1; 
              else 
                 alpha = -1; 
              end 
              
              if circ(m+ii) == circ(n+jj) 
                 beta = 1/12; 
              else 
                 beta = 1/24; 
              end 

              mp3mi = circ(m+3-ii);
              np3mj = circ(n+3-jj);

              eleTij = eleTij + alpha*beta*(b(mp3mi)*b(np3mj) + ...
                                c(mp3mi)*c(np3mj));

              end

           end

           eleT(m,n) = eleTij * (w(m)*w(n)/(Area*2));

        end
     end
end

function x = circ(x)  %  correct to circular arithmetic
  if(x<1)
    x=x+3;
  elseif(x>3)
    x=x-3;
  end
end
