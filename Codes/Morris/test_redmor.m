% Building a PC Surrogate for the reduced Morris test function and
% computing Sobol sensitivity indices.

close all;
clearvars;
uqlab;

% Input Model
MOpts.mFile = 'red_mor';
myModel = uq_createModel(MOpts);

IOpts.Marginals(1).Name = 'x1';
IOpts.Marginals(1).Type = 'Uniform';
IOpts.Marginals(1).Parameters = [0, 1]; 

IOpts.Marginals(2).Name = 'x2';
IOpts.Marginals(2).Type = 'Uniform';
IOpts.Marginals(2).Parameters = [0, 1]; 

IOpts.Marginals(3).Name = 'x3';
IOpts.Marginals(3).Type = 'Uniform';
IOpts.Marginals(3).Parameters = [0, 1]; 

IOpts.Marginals(4).Name = 'x4';
IOpts.Marginals(4).Type = 'Uniform';
IOpts.Marginals(4).Parameters = [0, 1]; 

my_Input = uq_createInput(IOpts);
%fid = fopen('itr_indices.txt','w');
%fmt = '%s %s %s %s %s %s %s %s %s'
%s1 ='n_itr';s2 ='STi_x1';s3='STi_x2';s4='STi_x3';s5='STi_x4';
%s6='UB1_x1';s7='UB1_x2';s8='UB1_x3';s9='UB1_x4';
%fprintf(fid,fmt,s1,s2,s3,s4,s5,s6,s7,s8,s9);
%fprintf(fid,'\n');
%fclose(fid);
%
%fmt = '%d %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f\n';

%for ns = 5:5:50
% Create a PCE of the model
PCEOpts.Type = 'Metamodel';
PCEOpts.MetaType = 'PCE';
PCEOpts.FullModel = myModel;
PCEOpts.Degree = 1:10;
PCEOpts.Method = 'LARS';
PCEOpts.ExpDesign.NSamples = 35;
myPCE = uq_createModel(PCEOpts);

% Sensitivity Analysis
%SobolOpts.Type = 'Sensitivity';
%SobolOpts.Method = 'Sobol';
% Specify the maximum order of the Sobol' indices to be calculated
%SobolOpts.Sobol.Order = 1;
% Specify the Sample size of each variable. Note that the total cost will
% equal (M+2)*SampleSize.
%SobolOpts.Sobol.SampleSize = 1e5;
%SobolOpts.Sobol.SampleSize = ns;
%PCESobolAnalysis = uq_createAnalysis(SobolOpts);
%PCESobolResults = PCESobolAnalysis.Results;
%STi = PCESobolResults.Total;

% Compute DGSM indices
X = myPCE.ExpDesign.X;
V = myPCE.PCE.Moments.Var;
dX = 1e-5;
G = red_mor(X);
dim = size(X,2);
nsam = size(X,1);
Gdx = zeros(nsam,dim);
nu = zeros(dim,1);
ub1 = zeros(dim,1);

for i = 1:dim
    X_new = X;
    X_new(:,i) = X_new(:,i) + dX;
    Gdx(:,i) = red_mor(X_new);
    dGdx = (Gdx(:,i) - G)./dX;
    dGdx2 = dGdx.^2;
    S = sum(dGdx2);
    nu(i) = S./nsam;
    ub1(i) = nu(i)./((pi.^2).*V);
end
%
%fid = fopen('itr_indices.txt','a');
%fprintf(fid,fmt,ns,STi(1),STi(2),STi(3),STi(4),ub1(1),ub1(2),ub1(3),ub1(4));
%fprintf(fid,'\n');
%fclose(fid);
%end

% ERROR Analysis

%data = load('val_pts.mat');
%vp = data.val_pts;
%Y_PCE = uq_evalModel(vp);
%Y_Model = red_mor(vp);
%NR = sum(((Y_Model-Y_PCE).^2)).^(0.5);
%DR = sum((Y_Model.^2)).^(0.5);
%E = float(NR)./float(DR);















