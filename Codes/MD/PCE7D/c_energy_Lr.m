close all
clear all

e = load('raw_energy/energy_1_80.txt');
Lr = load('raw_Lr/Lr_1_80.txt');
p = load('raw_params/params_1_80.txt');

s = size(e,1);

ref_e = zeros(s,3);
ref_Lr = zeros(s,1);
ref_p = zeros(s,7);
cou = 0;

for i = 1:s
  if (e(i,2) < 0 && e(i,3) > 0)  
      av = (abs(e(i,2))+abs(e(i,3)))./2;
      ch = (abs(abs(e(i,3)) - abs(e(i,2))).*100)./av;
      if (ch < 10.0)
          cou = cou + 1;
          ref_e(cou,:) = e(i,:);
          ref_Lr(cou,:) = Lr(i,:);
	  ref_p(cou,:) = p(e(i,1),:);
      end
  end
end

ref_e = ref_e(1:cou,:);
ref_Lr = ref_Lr(1:cou,:);
ref_p = ref_p(1:cou,:);

fid1 = fopen('ref_e_1_80.txt','w');
fid2 = fopen('ref_Lr_1_80.txt','w');
fid3 = fopen('ref_p_1_80.txt','w');

for i = 1:cou
  fprintf(fid1,'%d %10.3f %10.3f\n',ref_e(i,1),ref_e(i,2),ref_e(i,3));
  fprintf(fid2,'%5.3f\n',ref_Lr(i,1));
  fprintf(fid3,'%15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n',ref_p(i,1),ref_p(i,2),ref_p(i,3),ref_p(i,4),ref_p(i,5),...
                ref_p(i,6),ref_p(i,7));
end
