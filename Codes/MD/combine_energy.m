close all
clear all

set1 = load('energy_set1.txt');
set2 = load('energy_set2.txt');
set3 = load('energy_set3.txt');
set4 = load('energy_set4.txt');
com = zeros(160,3);
np = 5; % N for each set where N(d+1) computations were performed
tot = 160;

f = 1; % counter for the combined energy file
k = 1; % counter for individual energy files

while f <= tot

% set1
  in_f = f; in_k = k;
  fin_f = in_f+(np-1); fin_k = in_k+(np-1);   

  com(in_f:fin_f,1) = in_f:fin_f;
  com(in_f:fin_f,2:3) = set1(in_k:fin_k,2:3);

% set2
  in_f = fin_f+1;
  fin_f = in_f+(np-1);

  com(in_f:fin_f,1) = in_f:fin_f;
  com(in_f:fin_f,2:3) = set2(in_k:fin_k,2:3);
  
% set3
  in_f = fin_f+1;
  fin_f = in_f+(np-1);

  com(in_f:fin_f,1) = in_f:fin_f;
  com(in_f:fin_f,2:3) = set3(in_k:fin_k,2:3);
  
% set4
  in_f = fin_f+1;
  fin_f = in_f+(np-1);

  com(in_f:fin_f,1) = in_f:fin_f;
  com(in_f:fin_f,2:3) = set4(in_k:fin_k,2:3);
  
  f = f+4*np;
  k = k+np;
end

fid = fopen('energy_combined.txt','w');
fmt = '%02d %10.3f %10.3f\n';

for i = 1:tot
  fprintf(fid,fmt,com(i,1),com(i,2),com(i,3));
end

fclose(fid);
