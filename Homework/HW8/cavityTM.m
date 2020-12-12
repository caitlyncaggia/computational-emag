function cavityTM
% compute propagation constants for waveguide cross section modeled with
% triangular cells
%
% September 23, 2018 A. F. Peterson

global pcetond xy;
global icell;

% read mesh from file 'cylfil.txt'
nnodes = dlmread('cylfil.txt','', [0,0,0,0]);
ncells = dlmread('cylfil.txt','', [0,1,0,1]);
ninner = dlmread('cylfil.txt','', [0,2,0,2]);
nouter = dlmread('cylfil.txt','', [0,3,0,3]);
xy=dlmread('cylfil.txt','', [1,0,nnodes,1]);
nstart=nnodes + 1;
nend=nstart + ncells - 1;
pcetond=dlmread('cylfil.txt','', [nstart,0,nend,5]);

nstart=nend + 1;
nend=nstart + ncells - 1;
er=dlmread('cylfil.txt','', [nstart,0,nend,0]);

% initialize variables
rerr = 1.0e-8;
aerr = 1.0e-12;
nunks = nnodes - ninner - nouter;
W=zeros(nunks);
Y=zeros(nunks);

% fill global matrix one cell at a time
for icell=1:ncells
    %    compute 6 by 6 element matrix for cell 'icell'

    [eleS, eleT] = elemat(rerr,aerr);
     
    % add contributions from cell 'icell' to global matrix
    for ii=1:6
        ig=pcetond(icell,ii);
        for jj=1:6
            jg=pcetond(icell,jj);
            if(ig <= nunks)
                if(jg <= nunks)
                    W(ig,jg)=W(ig,jg) + er(icell)*eleS(ii,jj);
                    Y(ig,jg)=Y(ig,jg) + er(icell)*eleT(ii,jj);
                end
            end
        end
    end
end

fid = fopen('eigfil.txt', 'wt');
E = eig(W,Y); % use [V,E] = eig(W,Y) to get eigenvectors as well

str = 'TM resonant wavenumbers: ';
fprintf(fid,'%s \n',str);
for ii=1:nunks
    reaE=real(sqrt(E(ii)));
    imaE=imag(sqrt(E(ii)));
    fprintf(fid,'%6d %15.14g %15.14g\n',ii, reaE,imaE);
end

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
           eleT(ii,jj)=rtriad('Tint',rerr,aerr);
        end
     end
end
