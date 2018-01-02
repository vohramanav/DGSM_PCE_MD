close all
clear all

set1 = load('energy_set1.txt');
set2 = load('energy_set2.txt');
com = zeros(80,3);
np = 5; % N for each set where N(d+1) computations were performed

f = 1; % counter for the combined energy file
k = 1; % counter for individual energy files

while f <= 80
  com(f:f+(np-1),1) = f:f+(np-1);
  com(f:f+(np-1),2:3) = set1(k:k+(np-1),2:3);
  com(f+np:(f+np)+(np-1),1) = f+np:(f+np)+(np-1);
  com(f+np:(f+np)+(np-1),2:3) = set2(k:k+(np-1),2:3);
  f = f+2*np;
  k = k+np;
end

fid = fopen('energy_combined.txt','w');
fmt = '%02d %10.3f %10.3f\n';

for i = 1:80
  fprintf(fid,fmt,com(i,1),com(i,2),com(i,3));
end

fclose(fid);
