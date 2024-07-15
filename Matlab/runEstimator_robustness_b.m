function runEstimator_robustness_b(id)

%Argument id matters for the random starting points
eSubsVals = [1.25;(1.25+2.3)/2;2.3];
gammaVals = (eSubsVals-1)./eSubsVals;

rng(1234+id);


%%% Estimating with regional variation
%Logging
diary OFF
diaryFileName = ['estimationOutput/estimation_rob_b' num2str(id) '.txt'];
delete(diaryFileName);
diary(diaryFileName);
%Loading regional moments and creating another estimator
load('../CleanData/moments_regions.mat');
ME = ModelEstimator(moments,V,pmeDataSummary);
%No further random shifts in starting points
ME.centerX = [1.519483984809266   -0.05746128964018259     0.6345595333994322     -4.840192719706102     0.2748069312196685       1.27676233415676   -0.09986608619501373     0.4674784486931917    -0.2558724996252222     0.5442189398921765      1.863367287023733  0.0006913340565242012      2.018431383085041     -6.397982725196995   -0.05923795416751185       1.47392314513044   -0.09280245819679893    -0.6381696897930241     -1.851903340765297     0.5352884630818725      1.507147147572896    -0.3414143620665406   -0.09618692899202715      -18.9265233194618    0.08373299925924343      1.877114232210329   -0.03766199307480714       2.59402481899324     -7.677164151569853    -0.6339114914213518     0.1428299019156748     0.9277190461179347     -1.238019830127258     -1.267044143094966      6.307463421799779     -7.264891546289503      1.270469152861867     -11.52364482980506      18.08448190532703]';
    if id == 1
        ME.rangeSampleX = 0;
        ME.initialSampleDraws = 1;
    else        
        ME.rangeSampleX = 0.4;
        ME.initialSampleDraws = 1;
    end

gamma = gammaVals(2);
for iR = 1:ME.numRegions
    ME.regionModels{iR}.gamma = gamma;
    ME.regionModels{iR}.bD = - ME.regionModels{iR}.bD .* ...
        (ME.regionModels{iR}.r./ME.regionModels{iR}.lambda_for+1);
end
disp('Gamma:');
disp(gamma);

%Use derivative-based method
ME.algOption = ModelEstimator.alg_fminunc;
%Estimate
results = ME.run();
diary off
end

