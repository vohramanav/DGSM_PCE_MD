
fid = fopen('params_5_5.txt','w');
fmt = '%14.8f %14.8f %14.8f %14.8f %14.8f %14.8f %14.8f\n';

for i = 1:40
  fprintf(fid,fmt,params_5_5(i,1),params_5_5(i,2),params_5_5(i,3),params_5_5(i,4),...
                  params_5_5(i,5),params_5_5(i,6),params_5_5(i,7));  
end

fclose(fid);
