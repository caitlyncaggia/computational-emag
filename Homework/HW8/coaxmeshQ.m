function coaxmeshQ

% Generate quadratic-curved triangular-cell model of coaxial 
% cross section
%
% Code generates some edge-related pointers; not all are dumped to file
% in present version
%
% A. F. Peterson  October 25, 2018
%
% initialization

      global nnods  x  y;
      global ppaton;
      global nlast plast;
      global nnew pnew;
      global pedtop pedton;
      global nedgs;
      global nptchs;
      global pouted;

 nlmat = 4000;
 nsmat = 20;
 plast(nsmat)=0;

 nptch(nsmat)=0;
 pouted(nsmat)=0;
 x(nlmat)=0;
 y(nlmat)=0;
 pnew(nsmat)=0;
 ppaton = zeros(nlmat,3);
 pedtop = zeros(nlmat,2);
 pedton = zeros(nlmat,2);
 innerbnd(nsmat)=0;
 outerbnd(nsmat)=0;
 er(nsmat)=0;

%--------------------------------------------------------------------
% part one:  nodes on inner boundary for second-order mesh
%--------------------------------------------------------------------

 r = input('give inner cylinder radius in wavelengths  ');
 % err = input('give real part of epsilon-r  ');
 % eri = input('give imaginary part of epsilon-r  ');
 % er(1)=err+1j*eri;

 ncellinner = input('give number of cells along inner radius  ');
 nnods = 2*ncellinner;

%     define the starting edge nodes and set up pointer 'pstart'

      theta=2*pi/nnods;
      for ii=1:nnods
        x(ii)=r*cos(theta*(ii-1));
        y(ii)=r*sin(theta*(ii-1));
        plast(ii)=ii;
        innerbnd(ii)=ii;
      end
      ninnerbnd=nnods;
      
      nptchs=0;
      nlast=nnods;
            
%--------------------------------------------------------------------
% part two:     generate layers
%--------------------------------------------------------------------

      nlays = input('number of layers desired  ');

      if(nlays > 0)
      for ii=1:nlays
         str = 'for coating # ';
         disp(str);
         disp(ii);
         ri = input('give desired radius of circle:  ');

         er(ii) = input('give relative permittivity ');
%        err = input('give real relative permittivity ');
%        eri = input('give imag relative permittivity ');       
%        er(ii) = err + 1j*eri;

%     at this point in the computation, 'nlast' represents the number
%     of node points along the outer perimeter of the preceding layer,
%     'plast' points to these node numbers in sequential order from
%     the positive x-axis

%     'nflag' is 1 if a preceding node happens to sit on the y-axis,
%     which is the case if (nlast/4) is EVEN

        ndummy = fix((nlast*0.25)+.1);
        if(rem(ndummy,2) == 0)
            disp('preceding node on y axis');
           nflag=1;
           else
            disp('preceding node NOT on y axis');
           nflag=0;
        end

%       'nnew' is the number of points to be generated along radius 'r'

        nnew=nlast+4+4*nflag;
        str = ['generating new points on annulus: ',num2str(nnew)];
        disp(str);

%       assign new nodes along radius 'ri'

        for jj=1:nnew
           pnew(jj)=nnods+jj;
           theta=(jj-1)*2*pi/nnew;
           x(pnew(jj))=ri*cos(theta);
           y(pnew(jj))=ri*sin(theta);
           outerbnd(jj)=pnew(jj);
        end
        nouterbnd=nnew;

%       wrap 'plast' and 'pnew' around the circle for simplicity

        plast(nlast+1)=plast(1);
        plast(nlast+2)=plast(2);
        pnew(nnew+1)=pnew(1);
        pnew(nnew+2)=pnew(2);
        nnods=nnods+nnew;

%       create pointers from new patches to nodes, add mid-side nodes
%       in center of annulus as required

        ndummy = annpat(nflag);
        nptch(ii)=nptchs;

%       we have finished annulus 'i'; need to replace 'plast' with 'pnew'

        nlast=nnew;
        for jj=1:nlast
           plast(jj)=pnew(jj);
        end

      end
      end

%      disp(x(1:nnods));
%      disp('total number of nodes = ');
%      disp(nnods);
%      disp(ppaton(1:nptchs,1:3));
      
%--------------------------------------------------------------------
% part three:  renumber nodes so boundary nodes are at the end of the list
%--------------------------------------------------------------------

% At this point, the pointer 'innerbnd(1:ninnerbnd)' has the list of
% nodes on the inner boundary, while 'outerbnd(1:nounterbnd)' has the
% list of nodes on the outer boundary.  These should be numbered
% consecutively.

% the inner boundary nodes are the first 'ninnerbnd' in the node list,
% but the outer boundary nodes are not at the end

    disp('nodes on inner and outer bondaries: ');
    disp(ninnerbnd);
    disp(nouterbnd);
    disp('starting and ending nodes for inner bondaries: ');
    disp(innerbnd(1));
    disp(innerbnd(ninnerbnd));
    disp('starting and ending nodes for outer bondaries: ');
    disp(outerbnd(1));
    disp(outerbnd(nouterbnd));
    
% place outer boundary nodes at the end of the list

      nstart=outerbnd(1);
      for ii=1:nouterbnd
          x(nnods+ii)=x(nstart+ii-1);
          y(nnods+ii)=y(nstart+ii-1);
      end
      
% shift nodes back in the list to fill up the gap

      for ii=nstart:nnods
          x(ii)=x(ii+nouterbnd);
          y(ii)=y(ii+nouterbnd);
      end  
      
% fix cell to node pointer
      
      nend = outerbnd(nouterbnd);   %   nend = nstart+nouterbnd
      nextra=nnods-nend;
      str = ['adding index to outer boundary nodes: ',num2str(nextra)];
      disp(str);
      str = ['subtracting from other nodes: ',num2str(nouterbnd)];
      disp(str);
      
      for ii=1:nptchs
          for jj=1:6
             ndummy=ppaton(ii,jj);
             if((ndummy >= nstart) && (ndummy <= nend))
                ppaton(ii,jj)=ndummy+nextra;
             elseif(ndummy > nend)
                ppaton(ii,jj)=ndummy-nouterbnd;
             end
          end
      end
      
      for ii=1:nouterbnd
          outerbnd(ii)=outerbnd(ii)+nextra;
      end

% place inner boundary nodes at the very end of the node list

      for ii=1:ninnerbnd
          x(nnods+ii)=x(ii);
          y(nnods+ii)=y(ii);
      end
      
% shift nodes back in the list to fill up the gap

      for ii=1:nnods
          x(ii)=x(ii+ninnerbnd);
          y(ii)=y(ii+ninnerbnd);
      end  
         
      % fix cell to node pointer
      
      nextra=nnods-ninnerbnd;
      str = ['adding index to inner boundary nodes: ',num2str(nextra)];
      disp(str);
      str = ['subtracting from other nodes: ',num2str(ninnerbnd)];
      disp(str);

      for ii=1:nptchs
          for jj=1:6
             ndummy=ppaton(ii,jj);
             if(ndummy <= ninnerbnd)
                ppaton(ii,jj)=ndummy+nextra;
             else
                ppaton(ii,jj)=ndummy-ninnerbnd;
             end
          end
      end
 
      for ii=1:ninnerbnd
          innerbnd(ii)=innerbnd(ii)+nextra;
      end

      for ii=1:nouterbnd
          outerbnd(ii)=outerbnd(ii)-ninnerbnd;
      end

%--------------------------------------------------------------------
% part four:  generate pointers associated with the edges
%--------------------------------------------------------------------

      nedgs=0;

%   initialize pointers to -1 to flag edge-of-model

      for ii=1:nlmat
      pedtop(ii,2) = -1;
      end

%     each patch has three associated edges -- add data to database

      for ii=1:nptchs
      ndummy = adedge(ii,ppaton(ii,1),ppaton(ii,2));
      ndummy = adedge(ii,ppaton(ii,2),ppaton(ii,3));
      ndummy = adedge(ii,ppaton(ii,3),ppaton(ii,1));
      end

%--------------------------------------------------------------------
%  part five -- compute edge lengths as though straight
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
%  part six -- compute interior angles as though straight
%--------------------------------------------------------------------

  ang(nptchs,3)=0;
  for ii=1:nptchs
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
  for ii=1:nptchs
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
%  part seven -- dump relevant data to file
%--------------------------------------------------------------------

 fid = fopen('cylfil.txt', 'wt');

 fprintf(fid,'%6d %6d %6d %6d\n',nnods,nptchs, ...
                                 ninnerbnd,nouterbnd);

 for ii=1:nnods
   fprintf(fid,'%15.14g %15.14g\n',x(ii),y(ii));
 end
  
 for ii=1:nptchs
   fprintf(fid,'%6d %6d %6d %6d %6d %6d\n',ppaton(ii,1), ppaton(ii,2),...
       ppaton(ii,3), ppaton(ii,4), ppaton(ii,5), ppaton(ii,6));
 end

 for ii=1:ninnerbnd
     fprintf(fid,'%6d \n',innerbnd(ii));
 end
 
  for ii=1:nouterbnd
     fprintf(fid,'%6d \n',outerbnd(ii));
  end
  
 nflag=1;
 for ii=1:nptchs
   if(ii > nptch(nflag))
     nflag=nflag+1;
   end
%   err=real(er(nflag));
%   eri=imag(er(nflag));
%   fprintf(fid,'%15.14g %15.14g\n',err,eri);
   fprintf(fid,'%15.14g\n',er(nflag));
 end
  
% for ii=1:nedgs;
%   fprintf(fid,'%6d %6d\n',pedtop(ii,1),pedtop(ii,2));
% end

 %for ii=1:nedgs;
 %  fprintf(fid,'%6d %6d\n',pedton(ii,1),pedton(ii,2));
 %end

 disp('data placed in file: cylfil.txt');
 disp(' ');
 disp('number of nodes:   ');
 disp(nnods);
 disp('number of edges:   ');
 disp(nedgs);
 disp('number of patches: ');
 disp(nptchs);
 disp('  ');
 disp('Data in file:  nnods,nedgs,nptchs,ninnerbnd,nouterbnd');
 disp('               x, y list');
 disp('               ppaton');
 disp('               list of nodes on inner boundary');
 disp('               list of nodes on outer boundary');
 disp('               permittivity of each patch');
 
end

%--------------------------------------------------------------------
%--------------------------------------------------------------------

function[np1] = annpat(nflag)

%  annpat: add patch pointers to the database for an annulus
%          between nodes in 'plast' and 'pnew'
%     --> this version for 6-node quadratic triangles
%     --> this version adds mid-side nodes along center of annulus

%  nflag=1 if there is a outer corner node directly on the y-axis

      global ppaton;
      global plast;
      global nnew pnew;
      global nptchs;
      global nnods  x  y;

      disp('call to annpat with nptchs = ');
      disp(nptchs);
      disp('call to annpat with nnods = ');
      disp(nnods);
      
      np1=1;
      np2=1;

%     compute the number of patch-pairs per quarter-plane

      nquarp=fix((nnew+.1)*0.125)-1;
      str = ['# of patch pairs per quarter = ',num2str(nquarp)];
      disp(str);

%        need new node newA between plast(np1) and pnew(np2) to start
         nnods=nnods+1;
         newA=nnods;
         nsaved=newA;
         x(newA)=0.5*(x(plast(np1))+x(pnew(np2)));
         y(newA)=0.5*(y(plast(np1))+y(pnew(np2)));

      for ii=1:4           %     work in quadrant ii
  
         for k=1:nquarp

%           add pointers to new patches to database

%           need new node newB between plast(np1) and pnew(np2+2)
            nnods=nnods+1;
            newB=nnods;
            x(newB)=0.5*(x(plast(np1))+x(pnew(np2+2)));
            y(newB)=0.5*(y(plast(np1))+y(pnew(np2+2)));
            ppaton(nptchs+1,1)=plast(np1);
            ppaton(nptchs+1,2)=pnew(np2);
            ppaton(nptchs+1,3)=pnew(np2+2);
            ppaton(nptchs+1,4)=pnew(np2+1);
            ppaton(nptchs+1,5)=newB;
            ppaton(nptchs+1,6)=newA;
  
%           need new node newA between pnew(np2+2) and plast(np1+2)
            nnods=nnods+1;
            newA=nnods;
            x(newA)=0.5*(x(plast(np1+2))+x(pnew(np2+2)));
            y(newA)=0.5*(y(plast(np1+2))+y(pnew(np2+2)));
            ppaton(nptchs+2,1)=plast(np1);
            ppaton(nptchs+2,2)=pnew(np2+2);
            ppaton(nptchs+2,3)=plast(np1+2);
            ppaton(nptchs+2,4)=newA;
            ppaton(nptchs+2,5)=plast(np1+1);
            ppaton(nptchs+2,6)=newB;
            
            nptchs=nptchs+2;
            np1=np1+2;
            np2=np2+2;
         end

%        add additional outer cell at end of quadrant

         if(ii < 4)
            nnods=nnods+1;
            newB=nnods;
            x(newB)=0.5*(x(plast(np1))+x(pnew(np2+2)));
            y(newB)=0.5*(y(plast(np1))+y(pnew(np2+2)));
         else
            newB=nsaved; % last cell of annulus wraps around
         end
         nptchs=nptchs+1;
         ppaton(nptchs,1)=plast(np1);
         ppaton(nptchs,2)=pnew(np2);
         ppaton(nptchs,3)=pnew(np2+2);
         ppaton(nptchs,4)=pnew(np2+1);
         ppaton(nptchs,5)=newB;
         ppaton(nptchs,6)=newA;
         np2=np2+2;
         newA=newB;
         
         if rem(ii,2) == 1
            if nflag == 0

% special case: add an inner patch that straddles the y-axis
% (should only occur on first annulus and not always)

            nnods=nnods+1;
            newB=nnods;
            x(newB)=0.5*(x(plast(np1+2))+x(pnew(np2)));
            y(newB)=0.5*(y(plast(np1+2))+y(pnew(np2)));
            nptchs=nptchs+1;
            ppaton(nptchs,1)=plast(np1);
            ppaton(nptchs,2)=pnew(np2);
            ppaton(nptchs,3)=plast(np1+2);
            ppaton(nptchs,4)=newB;
            ppaton(nptchs,5)=plast(np1+1);
            ppaton(nptchs,6)=newA;
            np1=np1+2;
            end
         end
      end
 
      disp('return from annpat with nptchs = ');
      disp(nptchs);
      disp('return from annpat with nnods = ');
      disp(nnods);

      return

end

%--------------------------------------------------------------------
%--------------------------------------------------------------------

function[nnew] = adedge(np,n1,n2)

%  add the edge between nodes 'n1' and 'n2' to the database
%
%     np: one of the patches associated with this edge
%     n1,n2: nodes associated with this edge
%     nedgs: number of edges already in the pointers

%     for new edge, add new row to pointers 'pedtop' and 'pedton'

%     if this edge is already in the database, add new info
%     about patch 'np' to the existing pointers

      global pedtop pedton;
      global nedgs;

%    find out if the edge is already in the database

      ne = 0;
      for k=1:nedgs
        if pedton(k,1) == n1 || pedton(k,2) == n1
          if pedton(k,1) == n2 || pedton(k,2) == n2
          ne=k;
          end
        end
      end

      if(ne == 0)           %  new edge

          nedgs=nedgs+1;
          pedton(nedgs,1)=n1;
          pedton(nedgs,2)=n2;
          pedtop(nedgs,1)=np;
          nnew=1;

      else                       %  edge already in database

          pedtop(ne,2)=np;
          nnew=0;

      end
      return
end




