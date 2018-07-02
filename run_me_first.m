addpath(pwd);

cd proposed/;
addpath(genpath(pwd));
cd ..;

cd AltMinAlg/;
addpath(genpath(pwd));
cd ..;

cd manopt/;
addpath(genpath(pwd));
cd ..;

% if you need CVX, please uncomment below.
%cd cvx/;
%addpath(genpath(pwd));
%cd ..;


[version, release_date] = hpopt_version();
fprintf('##########################################################\n');
fprintf('###                                                    ###\n');
fprintf('###           Welcome to HybridPrecodingOpt            ###\n');
fprintf('###      (version:%s, released:%s)        ###\n', version, release_date);
fprintf('###                                                    ###\n');
fprintf('##########################################################\n');