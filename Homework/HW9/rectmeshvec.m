function rectmeshvec
%
% Generate triangular-cell model of rectangular cavity with boundary
% edges at the end of the list
%
% January 19, 2013, adapted from quadmesh version April 4, 2012
% A. F. Peterson
% minor mods November 2018

% initialization

 global nedgs pedtop pedton;

 nsize = 1000;
 x = zeros(nsize);
 y = zeros(nsize);
 ppaton = zeros(nsize,3);
 ppatoe = zeros(nsize,3);
 pedtop = zeros(nsize,2);
 pedton = zeros(nsize,2);
 pctnd = zeros(nsize,6);
 point2 = zeros(nsize);

%  input description of cavity

 bigx = input('give cavity x dimension  ');
 bigy = input('give cavity y dimension  ');

% err = input('give real part of epsilon-r  ');
% eri = input('give imaginary part of epsilon-r  ');
% er(1)=err+1j*eri;

 ncx = input('give number of cells along x  ');
 ncy = input('give number of cells along y  ');

%--------------------------------------------------------------------
% part one: generate nodes and pointer 'npaton'
%--------------------------------------------------------------------

	  xlow=-0.5*bigx;
	  ylow=-0.5*bigy;
	  delx=bigx/ncx;
	  dely=bigy/ncy;
%     initialize the cell counter

      ncells=0;

%     initialize the nodes along left edge of plate

      ncy2=2*ncy;
      for iy=0:ncy2
	     x(iy+1)=xlow;
	     y(iy+1)=ylow+iy*dely*0.5;
      end
      nodes=ncy2+1;

	  nstart=1;
%	  nend=nodes;

%     loop through plate by columns of cells

      for ix=1:ncx

%     add a new column of cell nodes

	     for iy=0:ncy2
	        in=nodes+1+iy;
	        x(in)=xlow+delx*0.5;
	        y(in)=ylow+iy*dely*0.5;
         end
	     midst=nodes+1;
	     nodes=nodes+ncy2+1;
	     midend=nodes;
%
	     for iy=0:ncy2
	        in=nodes+1+iy;
	        x(in)=xlow+delx;
	        y(in)=ylow+iy*dely*0.5;
         end
	     nodes=nodes+ncy2+1;
	     nright=midend+1;

%     build connectivity -- cell to node
%
%     -- nodes on left run from 'nstart' to 'nend'
%     -- nodes in center run from 'midst' to 'midend'
%     -- nodes on right run from 'nright' to 'nodes'

         for iy=1:ncy
	        ncells=ncells+1;
	        pctnd(ncells,1)=nstart;
	        pctnd(ncells,2)=nright;
	        pctnd(ncells,3)=nright+2;
	        pctnd(ncells,4)=nright+1;
	        pctnd(ncells,5)=midst+1;
	        pctnd(ncells,6)=midst;

	        ncells=ncells+1;
	        pctnd(ncells,1)=nstart;
	        pctnd(ncells,2)=nright+2;
	        pctnd(ncells,3)=nstart+2;
	        pctnd(ncells,4)=midst+2;
	        pctnd(ncells,5)=nstart+1;
	        pctnd(ncells,6)=midst+1;

	        nstart=nstart+2;
	        midst=midst+2;
	        nright=nright+2;
         end

         nstart=nodes-ncy2;
		 xlow=xlow+delx;
      end

%      convert corner nodes into 'ppaton' to use legacy code
%      that follows!

      for i=1:ncells
	     ppaton(i,1)=pctnd(i,1);
	     ppaton(i,2)=pctnd(i,2);
	     ppaton(i,3)=pctnd(i,3);
      end
      
%--------------------------------------------------------------------
% part two:  generate pointers associated with the edges
%--------------------------------------------------------------------

      nedgs=0;

%   initialize pointers to -1 to flag edge-of-model

      for ii=1:nsize
         pedtop(ii,2) = -1;
      end

%   each patch has three associated edges -- add edges to database
%   if necessary and construct pointer 'ppatoe'

      for ii=1:ncells
         edgenum = adedg(ii,ppaton(ii,1),ppaton(ii,2));
         ppatoe(ii,1)=edgenum;
         edgenum = adedg(ii,ppaton(ii,2),ppaton(ii,3));
         ppatoe(ii,2)=edgenum;
         edgenum = adedg(ii,ppaton(ii,3),ppaton(ii,1));
         ppatoe(ii,3)=edgenum;
      end

      if(nedgs > nsize) 
	     disp('WARNING:  nedgs exceeds limit of arrays ')
      end

%--------------------------------------------------------------------
%  part three -- sort edges to put boundary edges at end of list
%--------------------------------------------------------------------

 nbound=0;
 ninter=0;
 for ii=1:nedgs
     if(pedtop(ii,2) == -1)
         nbound=nbound+1;
         point2(ii)=1000+nbound;
     else
         ninter=ninter+1;
         point2(ii)=ninter;
     end
 end
 for ii=1:nedgs
     if(point2(ii) > 1000)
     point2(ii)=point2(ii)-1000+ninter;
     end
 end
 
%   point2 points from original edge # to new edge # (interior edges first)
%
%   rearrange other pointer arrays

 for ii=1:ncells
     for jj=1:3
     index = ppatoe(ii,jj);
     ppatoe(ii,jj)=point2(index);
     end
 end
 
 pedton2=pedton;
 for ii=1:nedgs
     for jj=1:2
         index=point2(ii);
         pedton(index,jj)=pedton2(ii,jj);
     end
 end

 pedtop2=pedtop;
 for ii=1:ncells
    for jj=1:2
         index=point2(ii);
    pedtop(ii,jj)=pedtop2(index,jj);
    end
 end

%--------------------------------------------------------------------
%  part six -- adjust cell to edge pointer so it aligns with ppaton
%--------------------------------------------------------------------
 
% clean up 'ppatoe' so edges are synchronized with nodes in 'ppaton'

  for ii=1:ncells
      n1=ppaton(ii,1);
      n2=ppaton(ii,2);
      n3=ppaton(ii,3);

      pstart(1:3)=ppatoe(ii,1:3);
      pnew=zeros(3,1);
      
%  find edge 1 between n2 and n3

    for jj=1:3
      nd1=pedton(pstart(jj),1);
      nd2=pedton(pstart(jj),2);
      if((nd1 == n2) || (nd1 == n3))
          if((nd2 == n3) || (nd2 == n2))
              pnew(1)=ppatoe(ii,jj);
          end
      end
    end
%  find edge 2 between n3 and n1

    for jj=1:3
      nd1=pedton(pstart(jj),1);
      nd2=pedton(pstart(jj),2);
      if((nd1 == n3) || (nd1 == n1))
          if((nd2 == n1) || (nd2 == n3))
              pnew(2)=ppatoe(ii,jj);
          end
      end
    end
      
%  find edge 3 between n1 and n2

    for jj=1:3
      nd1=pedton(pstart(jj),1);
      nd2=pedton(pstart(jj),2);
      if((nd1 == n1) || (nd1 == n2))
          if((nd2 == n2) || (nd2 == n1))
              pnew(3)=ppatoe(ii,jj);
          end
      end
    end
    
    for jj=1:3
        if(pnew(jj) == 0)
            disp('error in edge sorting for patch ');
            disp(ii);
        else
            ppatoe(ii,jj)=pnew(jj);   
        end
    end
  end

%--------------------------------------------------------------------
%  part seven -- compute edge lengths
%--------------------------------------------------------------------

  min=99999;
  max=0;
  sum=0;
  length(nedgs)=0;
  for ii=1:nedgs
      n1=pedton(ii,1);
      n2=pedton(ii,2);
      length(ii)=sqrt((x(n2)-x(n1))^2+(y(n2)-y(n1))^2);
      if(length(ii) > max)
          max=length(ii);
      end
      if(length(ii) < min)
          min=length(ii);
      end
      sum=sum+length(ii);
  end
  average=sum/nedgs;
  
  disp('min, max, and average edge lengths:');
  disp(min);  disp(max);  disp(average);
  disp('ratio of max length to min length:');
  sum=max/min;
  disp(sum);

%--------------------------------------------------------------------
%  part eight -- compute interior angles
%--------------------------------------------------------------------

  nptchs=ncells;
  ang(nptchs,3)=0;
  for ii=1:ncells
      n1=ppaton(ii,1);
      n2=ppaton(ii,2);
      n3=ppaton(ii,3);
      
      v12x=x(n2)-x(n1);
      v12y=y(n2)-y(n1);
      mag12=sqrt(v12x^2+v12y^2);
      
      v13x=x(n3)-x(n1);
      v13y=y(n3)-y(n1);
      mag13=sqrt(v13x^2+v13y^2);
      
      v32x=x(n2)-x(n3);
      v32y=y(n2)-y(n3);
      mag32=sqrt(v32x^2+v32y^2);
      
%     ang1 is between edge 12 and edge 13

      dot=v12x*v13x + v12y*v13y;
      ang(ii,1)=acos(dot/(mag12*mag13))*180/pi;
      
%     ang2 is between edge 21 and edge 23

      dot=v12x*v32x + v12y*v32y;
      ang(ii,2)=acos(dot/(mag12*mag32))*180/pi;
      
%     ang3 is between edge 32 and edge 31

      dot=-(v32x*v13x + v32y*v13y);
      ang(ii,3)=acos(dot/(mag32*mag13))*180/pi;
      
%      str=[num2str(ang(ii,1)),'  ',num2str(ang(ii,2)),'  ',...
%           num2str(ang(ii,3))];
%      disp(str);
  end

  min=99999;
  max=0;
  for ii=1:ncells
      for jj=1:3
          angle=ang(ii,jj);
          if(angle < min)
              min=angle;
          end
          if(angle > max)
              max=angle;
          end
      end
  end
  
  disp('min and max interior angles:');
  disp(min);  disp(max);
  disp('ratio of max angle to min angle:');
  sum=max/min;
  disp(sum);

%--------------------------------------------------------------------
%  part four -- dump data to file
%--------------------------------------------------------------------

 fid = fopen('cylfil.txt', 'wt');

 disp('CYLFIL.TXT contains nnods,nptchs,nedgs,ninn');
 disp('                    x,y coordinates');
 disp('                    ppaton(1:nptchs,3)');
 disp('                    ppatoe(1:nptchs,3)');
 disp('                    er(1:nptchs)');
 disp('                    pedtop(1:nedgs,2)');
 disp('                    pedton(1:nedgs,2)');
 disp('  ');
 disp('Run PLOTTRIMESH.M to generate a plot of the mesh');

 fprintf(fid,'%6d %6d %6d %6d\n',nodes,ncells,nedgs,ninter);

 for ii=1:nodes
   fprintf(fid,'%6d %15.14g %15.14g\n',ii, x(ii),y(ii));
 end

 for ii=1:ncells
   fprintf(fid,'%6d %6d %6d %6d\n',ii,ppaton(ii,1),ppaton(ii,2),...
       ppaton(ii,3));
 end

 for ii=1:ncells
   fprintf(fid,'%6d %6d %6d %6d\n',ii,ppatoe(ii,1),ppatoe(ii,2),...
       ppatoe(ii,3));
 end
 
 err = 1;
 for ii=1:ncells
  fprintf(fid,'%6d %15.14g\n',ii,err);
 end

 for ii=1:nedgs
  fprintf(fid,'%6d %6d %6d\n',ii,pedtop(ii,1),pedtop(ii,2));
 end
 
 for ii=1:nedgs
   fprintf(fid,'%6d %6d %6d\n',ii,pedton(ii,1),pedton(ii,2));
 end

% make copy in file 'trifile.txt' for plotter
 
 gid = fopen('trifile.txt', 'wt');
 
  fprintf(gid,'%6d %6d\n',nodes,ncells);

 for ii=1:nodes
   fprintf(gid,'%6d %15.14g %15.14g\n',ii,x(ii),y(ii));
 end
  
 for ii=1:ncells
   fprintf(gid,'%6d %6d %6d %6d\n',ii,ppaton(ii,1),ppaton(ii,2),...
       ppaton(ii,3));
 end

 disp('data placed in file: cylfil.txt');
 disp(' ');
 disp('number of nodes:   ');
 disp(nodes);
 disp('number of edges:   ');
 disp(nedgs);
 disp('number of patches: ');
 disp(ncells);

end

%
%--------------------------------------------------------------------
%     subroutine to create edge-related database
%--------------------------------------------------------------------
%

function[edgenum] = adedg(npat,n1,n2)

%  add the edge between nodes 'n1' and 'n2' to the database
%
%     npat: one of the patches associated with this edge
%     n1,n2: nodes associated with this edge
%     nedgs: number of edges already in the pointers

%     for new edge, add new row to pointers 'pedtop' and 'pedton'

%     if this edge is already in the database, add new info
%     about patch 'np' to the existing pointers

      global nedgs pedtop pedton;

%    find out if the edge is already in the database

      edgenum = 0;
      for k=1:nedgs
        if pedton(k,1) == n1 || pedton(k,2) == n1
          if pedton(k,1) == n2 || pedton(k,2) == n2
          edgenum=k;
          end
        end
      end

      if(edgenum == 0)           %  new edge

          nedgs=nedgs+1;
          pedton(nedgs,1)=n1;
          pedton(nedgs,2)=n2;
          pedtop(nedgs,1)=npat;
          edgenum=nedgs;
          
         else               %  edge already in database

          pedtop(edgenum,2)=npat;

      end
end
      

