function data = readm(fname,m)
 %
 % read a matrix of data with m values on each line
 %
fid = fopen(fname);
data = fscanf(fid,'%g',[m,inf]);
data = data';
status = fclose(fid);

