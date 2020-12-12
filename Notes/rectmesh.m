function rectmesh
%
% Generate triangular-cell model of rectangular domain 
%
% nodes on line and nodes on boundary are at the end of the list
%
% September 16, 2018
%
% A. F. Peterson

% initialization

 x(3000) = 0;
 y(3000) = 0;
 xb(2000) = 0;
 yb(2000) = 0;
 ppaton(3000,3)=0;
 pctnd(3000,3)=0;

%  input description of cavity

 bigx = input('give cavity x dimension  ');
 bigy = input('give cavity y dimension  ');

% er = input('give real part of epsilon-r for substrate  ');

 ncx = input('give number of cells along x  ');
 ncy = input('give number of cells along y  ');
 
  nbound1=0; % carry over from static capacitance code
 
%--------------------------------------------------------------------
% part one: generate nodes (x,y) and pointer 'nctnd'
%--------------------------------------------------------------------

	  xlow=-0.5*bigx;
	  ylow=0.0;
	  delx=bigx/ncx;
	  dely=bigy/ncy;
      
%     initialize the cell counter

      ncells=0;
      ninterior=0;
      pabove(ncx+1)=0;
      pbelow(ncx+1)=0;

%     initialize the nodes along lower edge of plate

      nbound2=ncx+1;
      for ix=1:nbound2
	     xb(ix)=xlow+(ix-1)*delx;
	     yb(ix)=ylow;
         pbelow(ix)=5000+ix;
      end

%   loop through domain by rows of cells

    for iy=1:ncy

%     add a new row of nodes

      if(iy == ncy)  % entire top of domain is a boundary
          
	     for ix=0:ncx
            nbound2=nbound2+1;
	        xb(nbound2)=xlow+ix*delx;
	        yb(nbound2)=ylow+iy*dely;
            pabove(ix+1)=5000+nbound2;
         end

      else %  treat normally
          
         nbound2=nbound2+1;
	     xb(nbound2)=xlow;
	     yb(nbound2)=ylow+iy*dely;
         pabove(1)=5000+nbound2;
         
	     for ix=1:ncx-1
            ninterior=ninterior+1;
	        x(ninterior)=xlow+ix*delx;
	        y(ninterior)=ylow+iy*dely;
            pabove(ix+1)=ninterior;
         end
 
         nbound2=nbound2+1;
	     xb(nbound2)=xlow+ncx*delx;
	     yb(nbound2)=ylow+iy*dely;
         pabove(ncx+1)=5000+nbound2;
      end
      
%     build connectivity -- cell to node
%
%     -- nodes below: pbelow(1:ncx)
%     -- nodes above: pabove(1:ncx)

         for ix=1:ncx
	        ncells=ncells+1;
	        pctnd(ncells,1)=pbelow(ix);
	        pctnd(ncells,2)=pbelow(ix+1);
	        pctnd(ncells,3)=pabove(ix);

	        ncells=ncells+1;
	        pctnd(ncells,1)=pbelow(ix+1);
	        pctnd(ncells,2)=pabove(ix+1);
	        pctnd(ncells,3)=pabove(ix);
         end
         
%        update pointer to nodes below  

         for ix=0:ncx
             pbelow(ix+1)=pabove(ix+1);
         end
    end
    
%--------------------------------------------------------------------
% part two:  renumber master node list and update 'ppaton'
%--------------------------------------------------------------------

%  have 'ninterior' nodes in list (x,y) and 'bound2' nodes in (xb,yb)
%
%  pceton has nbound2 nodes shifted by 5000 for convenience

   for ii=1:nbound2
       x(ninterior+nbound1+ii)=xb(ii);
       y(ninterior+nbound1+ii)=yb(ii);
   end
       
%  convert corner nodes into 'ppaton' 

   for ii=1:ncells
     for jj=1:3
       if(pctnd(ii,jj) > 5000)
           newindex=ninterior+nbound1+pctnd(ii,jj)-5000;
           ppaton(ii,jj)=newindex;
       else    
	       ppaton(ii,jj)=pctnd(ii,jj);
       end
     end
   end
      
   nodes=ninterior+nbound1+nbound2;
   
%--------------------------------------------------------------------
%  part three -- dump data to file
%--------------------------------------------------------------------

 fid = fopen('cylfil.txt', 'wt');

 fprintf(fid,'%6d %6d %6d %6d\n',nodes,ncells,nbound2,nbound1);

  for ii=1:nodes
      fprintf(fid,'%6d %15.14g %15.14g\n',ii,x(ii),y(ii));
  end
 
  for ii=1:ncells     
    fprintf(fid,'%6d %6d %6d %6d\n',ii,ppaton(ii,1),ppaton(ii,2),...
       ppaton(ii,3));
  end
 
  er2=1;
  for ii=1:ncells
   fprintf(fid,'%6d %15.14g\n',ii,er2);
  end

 disp('data placed in file: cylfil.txt');
 disp(' ');
 disp('number of nodes:   ');
 disp(nodes);
 disp('number of cells: ');
 disp(ncells);

end
