function mesh2

%  generate a 1D FEM mesh for a dielectric slab with air on both sides
%  -- second order connectivity matrix

%  October 20, 2018  A. F. Peterson

fid = fopen('inputfil.txt', 'wt');

str = 'mesh for a 3-region problem,'; disp(str);
str = 'all dimensions are in free-space wavelengths'; disp(str);

sizeL =  input(' give size of air region to left of slab: ');
sizeS =  input('     give size of dielectric slab region: ');
sizeR =  input('give size of air region to right of slab: ');

cellsL = input(' give number of cells spanning left region: ');
cellsS = input('        give number of cells spanning slab: ');
cellsR = input('give number of cells spanning right region: ');

deltaL=0.5*sizeL/cellsL;
deltaS=0.5*sizeS/cellsS;
deltaR=0.5*sizeR/cellsR;

epsilonS =  input('give relative permittivity of slab: ');

n_cells = cellsL + cellsS + cellsR;
n_nodes = 2*n_cells + 1;

fprintf(fid,'%6d %6d\n',n_nodes,n_cells);

% generate list of nodes

 nodesL=2*cellsL+1;
 for ii=1:nodesL
     x = (ii-1)*deltaL;
     fprintf(fid,'%6d %15.14g\n',ii, x);
 end
 nodesS=2*cellsS+1;
 for ii=2:nodesS
     xs = x + (ii-1)*deltaS;
     jj=2*cellsL + ii;
     fprintf(fid,'%6d %15.14g\n',jj, xs);
 end
 nodesR=2*cellsR+1;
 for ii=2:nodesR
     x = xs + (ii-1)*deltaR;
     jj=2*(cellsL+cellsS) + ii;
     fprintf(fid,'%6d %15.14g\n',jj, x);
 end
 
% generate cell-to-node connectivity matrix
 
pcell_to_node(n_cells,3)=0;
for ii=1:cellsL
  pcell_to_node(ii,1)=2*ii-1;
  pcell_to_node(ii,2)=2*ii;
  pcell_to_node(ii,3)=2*ii+1;
end
for ii=1:cellsS
  pcell_to_node(cellsL+ii,1)=2*cellsL + 2*ii - 1;
  pcell_to_node(cellsL+ii,2)=2*cellsL + 2*ii;
  pcell_to_node(cellsL+ii,3)=2*cellsL + 2*ii + 1;
end
for ii=1:cellsR
  pcell_to_node(cellsL+cellsS+ii,1)=2*(cellsL + cellsS) + 2*ii-1;
  pcell_to_node(cellsL+cellsS+ii,2)=2*(cellsL + cellsS) + 2*ii;
  pcell_to_node(cellsL+cellsS+ii,3)=2*(cellsL + cellsS) + 2*ii+1;
end
for ii=1:n_cells
   fprintf(fid,'%6d %6d %6d %6d\n',ii, pcell_to_node(ii,1), ...
       pcell_to_node(ii,2), pcell_to_node(ii,3));
end

 % generate a list of relative permittivities by cell
 
 x = 1.0;
 for ii=1:cellsL
     fprintf(fid,'%6d %15.14g\n',ii, x);
 end
 for ii=1:cellsS
     jj=cellsL+ii;
     fprintf(fid,'%6d %15.14g\n',jj, epsilonS);
 end
 for ii=1:cellsR
     jj=cellsL+cellsS+ii;
     fprintf(fid,'%6d %15.14g\n',jj, x);
 end

 str = 'mesh has been placed in file INPUTFIL.TXT'; disp(str);

end

