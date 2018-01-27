close all
clear all

set1 = load('energy_setwise/energy_set1.txt');
set2 = load('energy_setwise/energy_set2.txt');
set3 = load('energy_setwise/energy_set3.txt');
set4 = load('energy_setwise/energy_set4.txt');
set5 = load('energy_setwise/energy_set5.txt');
set6 = load('energy_setwise/energy_set6.txt');
tot = 240;
com = zeros(tot,3);
np = 5; % N for each set where N(d+1) computations were performed

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
%  
%% set3
  in_f = fin_f+1;
  fin_f = in_f+(np-1);

  com(in_f:fin_f,1) = in_f:fin_f;
  com(in_f:fin_f,2:3) = set3(in_k:fin_k,2:3);
%  
%% set4
  in_f = fin_f+1;
  fin_f = in_f+(np-1);

  com(in_f:fin_f,1) = in_f:fin_f;
  com(in_f:fin_f,2:3) = set4(in_k:fin_k,2:3);
  
%% set5
  in_f = fin_f+1;
  fin_f = in_f+(np-1);

  com(in_f:fin_f,1) = in_f:fin_f;
  com(in_f:fin_f,2:3) = set5(in_k:fin_k,2:3);
  
%% set6
  in_f = fin_f+1;
  fin_f = in_f+(np-1);

  com(in_f:fin_f,1) = in_f:fin_f;
  com(in_f:fin_f,2:3) = set6(in_k:fin_k,2:3);
  
  f = f+6*np;
  k = k+np;
end

fid = fopen('combined_energy.txt','w');
fmt = '%02d %10.3f %10.3f\n';

for i = 1:tot
  fprintf(fid,fmt,com(i,1),com(i,2),com(i,3));
end

fclose(fid);
