function coaxmesh

% Generate triangular-cell model of coaxial cross section
%
% Adapted from August 1 2011 full circular domain   
%
% Code generates some edge-related pointers; not all are dumped to file
% in present version
%
% A. F. Peterson  October 9, 2017; mesh quality added July 30, 2018
%
% initialization

      global x y;
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
% part one:  nodes on inner boundary
%--------------------------------------------------------------------

 r = input('give inner cylinder radius in wavelengths  ');
 % err = input('give real part of epsilon-r  ');
 % eri = input('give imaginary part of epsilon-r  ');
 % er(1)=err+1j*eri;

 nnods = input('give number of nodes on inner radius  ');

%     define the starting edge nodes and set up pointer 'pstart'

      theta=2*pi/nnods;
      for ii=1:nnods
        x(ii)=r*cos(theta*(ii-1));
        y(ii)=r*sin(theta*(ii-1));
        plast(ii)=ii;
        innerbnd(ii)=ii;
      end
      ninnerbnd=nnods;

%     compute total area of the inner cylinder cross section and scale
%     nodes so that the areas agree with the desired circle

      theta=pi/nnods;
      area=nnods*r*r*cos(theta)*sin(theta);
      anew=pi*r*r;

%      str = 'total area of inner cylinder cross section is: ';
%      disp(str);
%      disp(area);

%      str = 'desired area based on radius is: ';
%      disp(str);
%      disp(aread);

%      jj = input('want different area for inner cylinder (1=yes) ');
%      if(jj == 1)
%         anew = input('give new area  ');

%        scale nodes to match new area

         scale=sqrt(anew/area);
         for ii=1:nnods
            x(ii)=x(ii)*scale;
            y(ii)=y(ii)*scale;
         end
%      end
      
      nptchs=0;
      nlast=nnods;
            
%--------------------------------------------------------------------
% part two:     generate model of layers
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
%     which is the case if (nlast/2) is EVEN

        ndummy = fix((nlast*0.5)+.1);
        if(rem(ndummy,2) == 0)
%            disp('preceding node on y axis');
           nflag=1;
           else
%            disp('preceding node NOT on y axis');
           nflag=0;
        end

%       'nnew' is the number of points to be generated along radius 'r'

        nnew=nlast+2+2*nflag;
%        disp('generating new points on annulus: ');
%        disp(nnew);

%       polygon radius is scaled to get desired area correct

        theta=pi/nnew;
        r=ri*sqrt(pi/(nnew*cos(theta)*sin(theta)));

%       assign new nodes along radius 'r'

        for jj=1:nnew
           pnew(jj)=plast(nlast)+jj;
           theta=(jj-1)*2*pi/nnew;
           x(pnew(jj))=r*cos(theta);
           y(pnew(jj))=r*sin(theta);
        end

%       wrap 'plast' and 'pnew' around the circle for simplicity

        plast(nlast+1)=plast(1);
        pnew(nnew+1)=pnew(1);

%       create pointers from new patches to nodes

        ndummy = annpat(nflag);

%       we have finished annulus 'i'; need to update 'nnods' and replace
%       'plast' with 'pnew'

        nptch(ii)=nptchs;
        nnods=nnods+nnew;
        nlast=nnew;
        for jj=1:nlast
           plast(jj)=pnew(jj);
        end

      end
      end

      for ii=1:nlast
          outerbnd(ii)=plast(ii);
      end
      nouterbnd=nlast;
      
%--------------------------------------------------------------------
% part three:  renumber nodes so boundary nodes are at the end of the list
%--------------------------------------------------------------------

    disp('nodes on inner and outer bondaries: ');
    disp(ninnerbnd);
    disp(nouterbnd);
    
% place inner boundary nodes at the end

      for ii=1:ninnerbnd
          x(nnods+ii)=x(ii);
          y(nnods+ii)=y(ii);
      end
      for ii=1:nnods
          x(ii)=x(ii+ninnerbnd);
          y(ii)=y(ii+ninnerbnd);
      end          
      
      for ii=1:nptchs
          for jj=1:3
             ndummy=ppaton(ii,jj);
             if(ndummy > ninnerbnd)
                ppaton(ii,jj)=ndummy-ninnerbnd;
             else
                ppaton(ii,jj)=ndummy+nnods-ninnerbnd;
             end
          end
      end

      for ii=1:ninnerbnd
          innerbnd(ii)=innerbnd(ii)+nnods-ninnerbnd;
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
% part five: clean up 'ppaton' so nodes appear in CCW orientation
%--------------------------------------------------------------------

  for ii=1:nptchs
      n1=ppaton(ii,1);
      n2=ppaton(ii,2);
      n3=ppaton(ii,3);
      xmid=(x(n1)+x(n2)+x(n3))/3;
      ymid=(y(n1)+y(n2)+y(n3))/3;
      
      phi1=atan2(y(n1)-ymid,x(n1)-xmid);
      if(phi1 < 0) 
          phi1=phi1+2*pi;
      end
      
      phi2=atan2(y(n2)-ymid,x(n2)-xmid);
      if(phi2 < 0) 
          phi2=phi2+2*pi;
      end
      
      phi3=atan2(y(n3)-ymid,x(n3)-xmid);
      if(phi3 < 0) 
          phi3=phi3+2*pi;
      end
      
      if(phi1 < phi2)
        if(phi1 < phi3)
          if(phi2 > phi3)
%              disp('changing order in PPATON');
              ntemp=n2;
              n2=n3;
              n3=ntemp;
          end
        end
      else
          if(phi2 > phi3)
%              disp('changing order in PPATON');
              ntemp=n1;
              n1=n3;
              n3=ntemp;
          elseif(phi3 > phi1)
%              disp('changing order in PPATON');
              ntemp=n1;
              n1=n2;
              n2=ntemp;
          end
      end
      ppaton(ii,1)=n1;
      ppaton(ii,2)=n2;
      ppaton(ii,3)=n3;      
  end

%--------------------------------------------------------------------
%  part six -- compute edge lengths
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
%  part seven -- compute interior angles
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
%  part eight -- dump relevant data to file
%--------------------------------------------------------------------

 fid = fopen('cylfil.txt', 'wt');

 disp('CYLFIL.TXT contains nnods,nptchs,ninnerbnd,nouterbnd');
 disp('                    x,y coordinates');
 disp('                    ppaton(1:nptchs,3)');
 disp('                    er(1:nptchs)');
 disp('                    innerbnd(1:ninnerbnd),  outerbnd(1:nouterbnd)');
 
 fprintf(fid,'%6d %6d %6d %6d\n',nnods,nptchs, ...
                                     ninnerbnd,nouterbnd);

 for ii=1:nnods
   fprintf(fid,'%6d %15.14g %15.14g\n',ii,x(ii),y(ii));
 end
  
 for ii=1:nptchs
   fprintf(fid,'%6d %6d %6d %6d\n',ii,ppaton(ii,1),ppaton(ii,2),...
       ppaton(ii,3));
 end

 % make copy in file 'trifile.txt' for plotter
 
 gid = fopen('trifile.txt', 'wt');
 
  fprintf(gid,'%6d %6d\n',nnods,nptchs);

 for ii=1:nnods
   fprintf(gid,'%6d %15.14g %15.14g\n',ii,x(ii),y(ii));
 end
  
 for ii=1:nptchs
   fprintf(gid,'%6d %6d %6d %6d\n',ii,ppaton(ii,1),ppaton(ii,2),...
       ppaton(ii,3));
 end

 nflag=1;
 for ii=1:nptchs
   if(ii > nptch(nflag))
     nflag=nflag+1;
   end
%   err=real(er(nflag));
%   eri=imag(er(nflag));
%   fprintf(fid,'%15.14g %15.14g\n',err,eri);
   fprintf(fid,'%6d %15.14g\n',ii,er(nflag));
 end
 
  for ii=1:ninnerbnd
     fprintf(fid,'%6d \n',innerbnd(ii));
 end
 
  for ii=1:nouterbnd
     fprintf(fid,'%6d \n',outerbnd(ii));
  end
  
% for ii=1:nedgs;
%   fprintf(fid,'%6d %6d\n',pedtop(ii,1),pedtop(ii,2));
% end

 %for ii=1:nedgs;
 %  fprintf(fid,'%6d %6d\n',pedton(ii,1),pedton(ii,2));
 %end

 disp(' ');
 disp('number of nodes:   ');
 disp(nnods);
 disp('number of edges:   ');
 disp(nedgs);
 disp('number of patches: ');
 disp(nptchs);

end

%--------------------------------------------------------------------
%--------------------------------------------------------------------

function[np1] = annpat(nflag)

%  annpat: add patch pointers to the database for an annulus
%          between nodes in 'plast' and 'pnew'

%  nflag=1 if there is a node directly on the y-axis

      global ppaton;
      global plast;
      global nnew pnew;
      global nptchs;

%      disp('call to annpat with nptchs = ');
%      disp(nptchs);
      
      np1=1;
      np2=1;

%     compute the number of patch-pairs per quarter-plane

      nquarp=fix((nnew+.1)*0.25)-1;

      for ii=1:4           %     work in quadrant ii

         for k=1:nquarp

%           add pointers to new patches to database

            ppaton(nptchs+1,1)=plast(np1);
            ppaton(nptchs+1,2)=pnew(np2);
            ppaton(nptchs+1,3)=pnew(np2+1);
            ppaton(nptchs+2,1)=plast(np1);
            ppaton(nptchs+2,2)=plast(np1+1);
            ppaton(nptchs+2,3)=pnew(np2+1);
            nptchs=nptchs+2;
            np1=np1+1;
            np2=np2+1;
         end

%        add additional cell at end of quadrant

         nptchs=nptchs+1;
         ppaton(nptchs,1)=plast(np1);
         ppaton(nptchs,2)=pnew(np2);
         ppaton(nptchs,3)=pnew(np2+1);
         np2=np2+1;

         if rem(ii,2) == 1
            if nflag == 0

%           add another patch that straddles the y-axis

            nptchs=nptchs+1;
            ppaton(nptchs,1)=plast(np1);
            ppaton(nptchs,2)=plast(np1+1);
            ppaton(nptchs,3)=pnew(np2);
            np1=np1+1;
            end
         end
      end
 
%      disp('return from annpat with nptchs = ');
%      disp(nptchs);

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




