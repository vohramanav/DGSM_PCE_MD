close all
clear all

set1 = load('Lr_setwise/Lr_set1.txt');
set2 = load('Lr_setwise/Lr_set2.txt');
set3 = load('Lr_setwise/Lr_set3.txt');
set4 = load('Lr_setwise/Lr_set4.txt');
tot = 160;
com = zeros(tot,3);
np = 5; % N for each set where N(d+1) computations were performed

f = 1; % counter for the combined energy file
k = 1; % counter for individual energy files

while f <= tot

% set1
  in_f = f; in_k = k;
  fin_f = in_f+(np-1); fin_k = in_k+(np-1);   

  com(in_f:fin_f,1) = set1(in_k:fin_k,1);

% set2
  in_f = fin_f+1;
  fin_f = in_f+(np-1);

  com(in_f:fin_f,1) = set2(in_k:fin_k,1);
%  
%% set3
  in_f = fin_f+1;
  fin_f = in_f+(np-1);

  com(in_f:fin_f,1) = set3(in_k:fin_k,1);
%  
%% set4
  in_f = fin_f+1;
  fin_f = in_f+(np-1);

  com(in_f:fin_f,1) = set4(in_k:fin_k,1);
  
  f = f+4*np;
  k = k+np;
end

fid = fopen('combined_Lr.txt','w');
fmt = '%5.3f\n';

for i = 1:tot
  fprintf(fid,fmt,com(i,1));
end

fclose(fid);
