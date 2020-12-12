function mesh1D

%  generate a 1D finite difference/element mesh

%  June 21, 2018  A. F. Peterson

fid = fopen('inputfil.txt', 'wt');

n_nodes = input('give number of nodes ');

size =  input('give interval size ');
delta = size/(n_nodes-1);

fprintf(fid,'%6d\n',n_nodes);

 for irow=1:n_nodes
     x = (irow-1)*delta;
     fprintf(fid,'%6d %15.14g\n',irow, x);
 end
 x = 1.0;
 for irow=1:n_nodes
     fprintf(fid,'%6d %15.14g\n',irow, x);
 end

end

