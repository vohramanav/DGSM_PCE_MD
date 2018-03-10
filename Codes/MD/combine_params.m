close all
clear all

set1 = load('params_setwise/params_set1.txt');
set2 = load('params_setwise/params_set2.txt');
set3 = load('params_setwise/params_set3.txt');
set4 = load('params_setwise/params_set4.txt');
set5 = load('params_setwise/params_set5.txt');
set6 = load('params_setwise/params_set6.txt');
set7 = load('params_setwise/params_set7.txt');
tot = 280;
com = zeros(tot,7);
np = 5; % N for each set where N(d+1) computations were performed
nsets = 7; % number of sets

f = 1; % counter for the combined energy file
k = 1; % counter for individual energy files

while f <= tot

% set1
  in_f = f; in_k = k;
  fin_f = in_f+(np-1); fin_k = in_k+(np-1);   

  com(in_f:fin_f,:) = set1(in_k:fin_k,:);

% set2
  in_f = fin_f+1;
  fin_f = in_f+(np-1);

  com(in_f:fin_f,:) = set2(in_k:fin_k,:);
%  
%% set3
  in_f = fin_f+1;
  fin_f = in_f+(np-1);

  com(in_f:fin_f,:) = set3(in_k:fin_k,:);
%  
%% set4
  in_f = fin_f+1;
  fin_f = in_f+(np-1);

  com(in_f:fin_f,:) = set4(in_k:fin_k,:);
  
%% set5
  in_f = fin_f+1;
  fin_f = in_f+(np-1);

  com(in_f:fin_f,:) = set5(in_k:fin_k,:);
  
%% set6
  in_f = fin_f+1;
  fin_f = in_f+(np-1);

  com(in_f:fin_f,:) = set6(in_k:fin_k,:);
  
%% set7
  in_f = fin_f+1;
  fin_f = in_f+(np-1);

  com(in_f:fin_f,:) = set7(in_k:fin_k,:);
  
  f = f+nsets*np;
  k = k+np;
end

fid = fopen('combined_params.txt','w');
fmt = '%14.8f %14.8f %14.8f %14.8f %14.8f %14.8f %14.8f\n';

for i = 1:tot
  fprintf(fid,fmt,com(i,1),com(i,2),com(i,3),com(i,4),com(i,5),com(i,6),com(i,7));
end

fclose(fid);
