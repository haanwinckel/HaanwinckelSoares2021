function runEstimator(id)

%Argument id matters for the random starting points
eSubsVals = [1.25;(1.25+2.3)/2;2.3];
gammaVals = (eSubsVals-1)./eSubsVals;

rng(1234+id);


%%% Estimating with regional variation
%Logging
diary OFF
diaryFileName = ['estimationOutput/estimation_' num2str(id) '.txt'];
delete(diaryFileName);
diary(diaryFileName);
%Loading regional moments and creating another estimator
load('../CleanData/moments_regions.mat');
ME = ModelEstimator(moments,V,pmeDataSummary);
%No further random shifts in starting points
ME.centerX = [1.596549050429667   -0.01758208624871862     0.6163109647415214     -4.223549819227054     0.2164484380698613      1.532213748660098    0.01802601408642855      1.767962450320705     0.8970893561680364     0.4511371630044891      1.848544965818439   -0.01906823881017811      2.057504242499291     -6.859206348021765    -0.2467667049062195        1.4714341365529   -0.08735163337521784    -0.4798526165314614     -1.883004097128369      0.658742716961244      1.446874574672487    -0.3749078460251892    -0.2352003075854106     -18.90788015815535     0.3193971606430518      1.878026990910316   -0.01251142184318767      1.950702286374133     -7.133906486210638    -0.4594563112240468     0.1625102422860964       0.92314676054629     -1.124613536447893     -1.197604972299227      6.174156839527336     -7.234845020141777      1.127332259506962     -11.51851564292062      16.59013063660213]';
    if id <= 150
        gammaCat = 2;
    elseif id <= 175
        gammaCat = 1;
    else
        gammaCat = 3;
    end
    if id == 101 || id == 151 || id == 176
        ME.rangeSampleX = 0;
        ME.initialSampleDraws = 1;
    else        
        ME.rangeSampleX = 0.1;
        ME.initialSampleDraws = 5;
    end

gamma = gammaVals(gammaCat);
for iR = 1:ME.numRegions
    ME.regionModels{iR}.gamma = gamma;
end
disp('Gamma:');
disp(gamma);

%Use derivative-based method
ME.algOption = ModelEstimator.alg_fminunc;
%Estimate
results = ME.run();
diary OFF
end

