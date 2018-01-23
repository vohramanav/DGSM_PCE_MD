close all;
clear all;

q1_1 = -15881.349;
q1_2 = 15849.305;
lr = 49.264;
%q2_1 = -19072.368;
%q2_2 = 19058.147;

q1_avg = (abs(q1_1)+abs(q1_2))./2.0;
%q2_avg = (abs(q2_1)+abs(q2_2))./2.0;
dt = 0.0005;
ps2s = 1e-12;
lc = 5.43;
ang2m = 1e-10;
w = 22.*lc;
l1 = 50.*lc;
l2 = 50.*lc;
ld = 500000;
ev2j = 1.602e-19;

q1 = (q1_avg*ev2j)./(2.0.*dt.*ld.*ps2s.*w.*w.*ang2m.*ang2m);
%dTdx1 = 50./((l1./2).*ang2m);
cte1 = 1.5938;
dTdx1 = (cte1.*50.0)./(lr.*lc.*ang2m);
%cte2 = 1.5938;
%dTdx2 = cte2./(lc.*ang2m);
%q2 = (q2_avg*ev2j)./(2.0.*dt.*ld.*ps2s.*w.*w.*ang2m.*ang2m);
%dTdx2 = 50./((l2./2).*ang2m);

k1 = q1./dTdx1;
%k2 = q2./dTdx2;
