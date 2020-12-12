function FDTDsin

% 2D propagation in a free-space parallel plate waveguide terminated with a
% short circuit load, sinusoidal steady state n=1 mode excitation
%
% August 16, 2018   A. F. Peterson

  mu = pi*4.0e-7;
  epsilon = 8.854e-12;
  freq = 1.0e06;
  c = 2.998e8;
  sugdel = (c/freq) * 0.075;
  str=['suggested cell dimension = ',num2str(sugdel)];  disp(str);
  
  L = 2500; %input('Give guide length ');
  A = 2500; %input('Give guide width  ');
  Nz = 100; %input('Number of cells along length ');
  Nx = 100; %input('Number of cells along width ');
  delz = L/Nz;
  delx = A/Nx;
  Nsamp = round((Nx+1)/2);
  Nzp1=Nz+1;

% initialize variables
  
  Hy=zeros(Nx,Nz);
  Ex=zeros(Nx,Nzp1);
  Ez=zeros(Nx,Nz);

% determine time step

  deltmax = 1/c/sqrt((1/delz)^2+(1/delx)^2);
  str = ['CFL time step = ',num2str(deltmax)];  disp(str);
  delt = 3e-8; %input('Give time step  ');
  
  dtomudz = delt/(mu*delz);
  dtomudx = delt/(mu*delx);
  dtoepdz = delt/(epsilon*delz);
  dtoepdx = delt/(epsilon*delx);

% march in time
%
% Ex(jj,ii) is used to denote Ex(jj,ii-1/2)
% Ez(jj,ii) is used to denote Ez(jj-1/2,ii)
  
  fid = fopen('snapshots.txt', 'wt');
  
   for kk = 1:250

%    update source   ------------------------------------------------------

     t = kk*delt;
     for jj = 1:Nx
%         x=(jj-0.5)*delx;
        Ex(jj,1)=sin(2*pi*freq*t);
%         Ex(jj,1)=cos(pi*x/A)*sin(2*pi*freq*t);
     end
%     disp(Ex(Nsamp,1));
     
%    update Hy at time kk+1/2  (Nx by Nz samples)  ------------------------

     for ii = 1:Nz     
       for jj = 1:Nx
         
         if(jj == 1) % zero Ez field at wall of guide jj
                
            Hy(jj,ii)=Hy(jj,ii)-dtomudz*(Ex(jj,ii+1)-Ex(jj,ii)) ...
                               +dtomudx*(Ez(jj+1,ii)          );
                           
         elseif(jj == Nx)  % zero Ez field at wall of guide jj+1
                
            Hy(jj,ii)=Hy(jj,ii)-dtomudz*(Ex(jj,ii+1)-Ex(jj,ii)) ...
                               +dtomudx*(           -Ez(jj,ii));
                           
         else % out in the middle of the mesh somewhere
                
            Hy(jj,ii)=Hy(jj,ii)-dtomudz*(Ex(jj,ii+1)-Ex(jj,ii)) ...
                               +dtomudx*(Ez(jj+1,ii)-Ez(jj,ii));
                           
         end
       end
     end

%    update Ex at time kk  (Nx by Nz-1 samples)  --------------------------

     for ii = 2:Nz     
       for jj = 1:Nx
                         
         Ex(jj,ii)=Ex(jj,ii)-dtoepdz*(Hy(jj,ii)-Hy(jj,ii-1));
                           
       end 
     end
       
%    update Ez at time kk (Nx-1 by Nz samples)  ---------------------------

     for ii = 1:Nz     
       for jj = 2:Nx
                 
         Ez(jj,ii)=Ez(jj,ii)+dtoepdx*(Hy(jj,ii)-Hy(jj-1,ii));  
            
       end
         
     end

% store snapshot down centerline in file 'fid'  ---------------------------

      str = ['result at time = ',num2str(t)];
      fprintf(fid,'%s \n',str);
      for ii=1:Nz
         z = delz * (ii-1);
         fprintf(fid,'%15.14g %15.14g\n',z,Ex(Nsamp,ii));
      end
      fprintf(fid,'%15.14g %15.14g\n',L,0.0);
      fprintf(fid,'\n\n');
      
      if t > 5.999e-6 && t < 6.001e-6
          figure
          z = delz .* (0:Nz);
          plot(z, Ex(Nsamp,:))
      end

   end
  
end

