function [TMCTLNcavity] = TMCTLNcavity

%

global pceton xy
global icell
global pcetoe
global er
global pedtoc
global pedton

nnodes = dlmread('cylfil.txt','', [0,0,0,0]);
ncells = dlmread('cylfil.txt','', [0,1,0,1]);
nedges = dlmread('cylfil.txt','', [0,2,0,2]);
% ninnernodes = dlmread('cylfil.txt','', [0,3,0,3]);

xy=dlmread('cylfil.txt','', [1,1,nnodes,2]);

nstart=nnodes + 1;
nend=nstart + ncells - 1;
pceton=dlmread('cylfil.txt','', [nstart,1,nend,3]);

nstart=nend + 1;
nend=nstart + ncells - 1;
pcetoe=dlmread('cylfil.txt','', [nstart,1,nend,3]);

nstart=nend + 1;
nend=nstart + ncells - 1;
er=dlmread('cylfil.txt','', [nstart,1,nend,1]);

nstart=nend + 1;
nend=nstart + nedges - 1;
pedtoc=dlmread('cylfil.txt','', [nstart,1,nend,2]);

nstart=nend + 1;
nend=nstart + nedges - 1;
pedton=dlmread('cylfil.txt','', [nstart,1,nend,2]);

% initialize variables
nunks = nedges;

W=zeros(nunks);
Y=zeros(nunks);

% loop through the cells, filling global matrix one cell at a time

for icell=1:ncells
    
    [eleS,eleT] = elemat(icell);
    node1 = pceton(icell,1);
    node2 = pceton(icell,2);
    node3 = pceton(icell,3);
    if node1>node2
        eleS(3,:) = -eleS(3,:);
        eleS(:,3) = -eleS(:,3);
        eleT(3,:) = -eleT(3,:);
        eleT(:,3) = -eleT(:,3);
    elseif node2>node3
        eleS(3,:) = -eleS(3,:);
        eleS(:,3) = -eleS(:,3);
        eleT(3,:) = -eleT(3,:);
        eleT(:,3) = -eleT(:,3);
    elseif node3>node1
        eleS(3,:) = -eleS(3,:);
        eleS(:,3) = -eleS(:,3);
        eleT(3,:) = -eleT(3,:);
        eleT(:,3) = -eleT(:,3);
    end
    
    
    %    add contributions from cell 'icell' to global matrices
    
    for ii=1:3
        for jj=1:3
            ig=pcetoe(icell,ii); % 'ig' is the global node for 'ii'
            jg=pcetoe(icell,jj);  % 'jg' is the global node for 'jj'
            W(ig,jg) = W(ig,jg) + er(icell)*eleS(ii,jj);
            Y(ig,jg) = Y(ig,jg) + er(icell)*eleT(ii,jj);
        end
    end
    
end


%  write results to file 'eigfil.txt'

fid = fopen('eigfil.txt', 'wt');

E = sort(eig(W,Y)); % use [V,E] = eig(W,Y) to get eigenvectors as well

str = 'TM resonant wavenumbers: ';
fprintf(fid,'%s \n',str);

for ii=1:nunks
    reaE=real(sqrt(E(ii)));
    imaE=imag(sqrt(E(ii)));
    fprintf(fid,'%6d %15.14g %15.14g\n',ii, reaE,imaE);
end
end

% ----------------------------------------------------------------------

function [eleS,eleT] = elemat(icell)

%  elemat: construct the 6 by 6 element matrices for the contributions of
%          basis and testing functions of the form
%                   grad T dot grad B    and    T times B
%          over a triangular cell
%
%          done by quadrature over integrands 'Sint' and 'Tint'

global pceton;
global xy;

eleS(3,3)=0;  eleT(3,3)=0;

%cell-to-node pointers
node1 = pceton(icell,1);
node2 = pceton(icell,2);
node3 = pceton(icell,3);

%x,y coordinates
x = [xy(node1,1); xy(node2,1); xy(node3,1)];
y = [xy(node1,2); xy(node2,2); xy(node3,2)];

%edge lengths
w = [sqrt((x(3)-x(2))^2 + (y(3)-y(2))^2); ...
    sqrt((x(3)-x(1))^2 + (y(3)-y(1))^2); ...
    sqrt((x(2)-x(1))^2 + (y(2)-y(1))^2)];

%new coordinate system
b = [y(2)-y(3); y(3)-y(1); y(1)-y(2)];
c = [x(3)-x(2); x(1)-x(3); x(2)-x(1)];

%area
A = abs(b(3)*c(1) - b(1)*c(3))*0.5;
    for m = 1:3
        for n = 1:3
            m1 = modind(m+1);
            m2 = modind(m+2);
            n1 = modind(n+1);
            n2 = modind(n+2);
            
            eleS(m,n) = ((w(m)*w(n))/(4*A^3))* ...
                (b(m1)*c(m2) - b(m2)*c(m1)) *(b(n1)*c(n2) - b(n2)*c(n1));
            eleT(m,n) = 0;
            for i = 1:2
               for j = 1:2
                   if(modind(m+i) == modind(n+j))
                       beta = 1/12;
                   else
                       beta = 1/24;
                   end
                   
                   if (i == j)
                       alpha = 1;
                   else
                       alpha = -1;
                   end
                   m3 = modind(m+3-i);
                   n3 = modind(n+3-j);
                   eleT(m,n) = eleT(m,n) + ...
                       alpha*beta*(b(m3)*b(n3)+c(m3)*c(n3));
               end
            end
            eleT(m,n) = eleT(m,n)*w(m)*w(n)/(2*A);

        end
    end
end

function [ind] = modind(oldInd)
    ind = mod((oldInd-1),3)+1;
end
    
