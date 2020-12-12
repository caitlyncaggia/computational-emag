function FDTDpml

% 2D propagation in a parallel plate waveguide w/ a dielectric slab
% 
% short circuit walls at z=0, x=0, x=W;   PML termination at z=L end
%
% one period sin line source excitation at user specified location
%
% dielectric slab at user specified location between z1 and z2
%
% PML region should be placed from some z3 to z=L
%
% September 23, 2018   A. F. Peterson
 
  L = input('Give guide length ');
  A = input('Give guide width  ');
  Nz = input('Number of cells along length ');
  Nx = input('Number of cells along width ');
  delz = L/Nz;
  delx = A/Nx;
  Nsamp = round((Nx+1)/2);
  Nzm1=Nz-1;

  Sz = input('Cell index of source along length? ');
  Sx = input('Cell index of source along width?  ');

% initialize variables
  
  BigT = 1.0e-6; % transient excitation of 1 microsecond duration
  mu = pi*4.0e-7;
  epsilon = 8.854e-12;
  c = 2.998e8;
  eta=sqrt(mu/epsilon);
  Hy=zeros(Nx,Nz);
  Ex=zeros(Nx,Nz);
  Ez=zeros(Nx,Nz);
  Dz=zeros(Nx,Nz);
  Dzmem=zeros(Nx,Nz);
  Jx=zeros(Nx,Nz);

% Ex(jj,ii) is used to denote Ex(jj,ii-1/2)
% Ez(jj,ii) is used to denote Ez(jj-1/2,ii)

  fid = fopen('snapshots.txt', 'wt');

% erx(ii) samples of permittivity at Ex-field locations
% erz(ii) samples of permittivity at Ez-field locations

  epsx=ones(Nz,1);
  epsz=ones(Nz,1);
  zstart = input('starting cell index in z for slab? ');
    zend = input('  ending cell index in z for slab? ');
   epslb = input('          epsilon-r for slab? ');

  str=['slab between ',num2str(zstart),' and ',num2str(zend)];
  fprintf(fid,'%s \n',str);
  str=['with epsilon-r = ',num2str(epslb)];
  fprintf(fid,'%s \n',str);
  fprintf(fid,'\n');
  
  for ii=zstart:zend
    if(ii == zstart) % values of epsx get average, epsz in slab
        epsx(ii)=(epslb+1)/2;
        epsz(ii)=epslb;
    elseif(ii == zend) % values of epsx get average, epsz out of slab
        epsx(ii)=(epslb+1)/2;
    else   % values of epsx and epsz are in slab
        epsx(ii)=epslb;
        epsz(ii)=epslb;
    end
  end
  
  %  PML region  sigmah(ii) for hy and ez
  %              sigmae(ii) for ex
  
  sigmax=zeros(Nz,1);
  sigmaz=zeros(Nz,1);
  zstart = input('starting cell index in z for PML? ');
    zend = input('  ending cell index in z for PML? ');
    if(zend > Nzm1)
        disp('error -- past the end of the mesh');
        zend = input('  ending cell index in z for PML? ');
    end
  maxsig = input('         max sigmaz for PML? ');

  str=['PML between ',num2str(zstart),' and ',num2str(zend)];
  fprintf(fid,'%s \n',str);

   num=(zend-zstart);
   for ii=zstart:zend
        sigmax(ii)=maxsig*((ii-zstart)/num)^2;
        sigmaz(ii)=maxsig*((ii-zstart+0.5)/num)^2;
        fprintf(fid,'%6d %15.14g\n',ii,sigmax(ii));
   end 

% determine time step based on background permittivity

  deltmax = 1/c/sqrt((1/delz)^2+(1/delx)^2);
  str = ['CFL time step = ',num2str(deltmax)];  disp(str);
  delt = input('Give time step  ');
  
  dtomudz = delt/(mu*delz);
  dtomudx = delt/(mu*delx);
  dtoepdz = delt/(epsilon*delz);

% march in time
    
   for kk = 1:500

%    update source   ------------------------------------------------------

     t = kk*delt;
     
         if(( t >= 0) && (t <= BigT))
         Jx(Sx,Sz)=125*sin(2*pi*t/BigT)/eta/delx/delz;
         else
         Jx(Sx,Sz)=0;
         end
     
%    store Dz from previous time step  ------------------------------------

     for ii=1:Nz
       for jj = 1:Nx
         Dzmem(jj,ii)=Dz(jj,ii); % store previous values for PML
       end
     end
     
%    update Hy at time kk+1/2  (Nx by Nz samples)  ------------------------

     for ii = 1:Nzm1     
      pmlcon=sigmaz(ii)*delt/(2*epsilon*epsz(ii)); %  zero outside PML
      for jj = 1:Nx
         
         if(jj == 1) % zero Ez field at wall of guide jj
                
            Hy(jj,ii)=Hy(jj,ii)*(1-pmlcon)/(1+pmlcon) ...
                     -dtomudz*(Ex(jj,ii+1)-Ex(jj,ii))/(1+pmlcon) ...
                     +dtomudx*(Ez(jj+1,ii)          )/(1+pmlcon);
                           
         elseif(jj == Nx)  % zero Ez field at wall of guide jj+1
                
            Hy(jj,ii)=Hy(jj,ii)*(1-pmlcon)/(1+pmlcon) ...
                     -dtomudz*(Ex(jj,ii+1)-Ex(jj,ii))/(1+pmlcon) ...
                     +dtomudx*(           -Ez(jj,ii))/(1+pmlcon);
                           
         else % out in the middle of the mesh somewhere
                
            Hy(jj,ii)=Hy(jj,ii)*(1-pmlcon)/(1+pmlcon) ...
                     -dtomudz*(Ex(jj,ii+1)-Ex(jj,ii))/(1+pmlcon) ...
                     +dtomudx*(Ez(jj+1,ii)-Ez(jj,ii))/(1+pmlcon);
                           
         end
       end
     end
     
%    Hy and Ez at ii=Nz are always zero to terminate mesh (PMC?)

%    update Ex at time kk  (Nx by Nz-1 samples)  --------------------------

     for ii = 2:Nz  % don't update Ex(x,1), keep at zero for PEC
      pmlcon=sigmax(ii)*delt/(2*epsilon*epsx(ii)); %  zero outside PML
      for jj = 1:Nx
                         
         Ex(jj,ii)=Ex(jj,ii)*(1-pmlcon)/(1+pmlcon) ...
                 -dtoepdz*(Hy(jj,ii)-Hy(jj,ii-1))/epsx(ii)/(1+pmlcon)...
                 -(delt/epsilon)*Jx(jj,ii)/epsx(ii)/(1+pmlcon);                           
       end 
     end
       
%    update Dz at time kk (Nx-1 by Nz samples)  ---------------------------

     for ii = 1:Nzm1  % value of Hy at Nz is zero   
       for jj = 2:Nx
                 
         Dz(jj,ii)=Dz(jj,ii)+delt*(Hy(jj,ii)-Hy(jj-1,ii))/delx;  
            
       end
     end
     
%    update Ez at time kk (Nx-1 by Nz samples)  ---------------------------

     for ii = 1:Nzm1    % value of Hy, Ez at Nz is zero       
      pmlcon=sigmaz(ii)*delt/(epsilon*epsz(ii))^2; % zero outside PML
      for jj = 2:Nx
                 
        Ez(jj,ii)=Ez(jj,ii)+(Dz(jj,ii)-Dzmem(jj,ii))/(epsilon*epsz(ii)) ...
                    +pmlcon*(Dz(jj,ii)+Dzmem(jj,ii))/2;  
            
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

   end
  
end

