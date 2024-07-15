function runBootstrapEstimation(id)

%Argument id matters for the random starting points
rng(1234+id);

%%% FIRST PART:
%%% Re-estimate models

x_point_est = [1.566696665648196   -0.01757280635684886     0.6315309255540186     -4.232890802380131     0.3016655930448394      1.514682216068157     0.0580088489270596      1.754421541560044     0.8703427612678325      0.499368991610191      1.893519053315942   0.003937642559902204       1.97149396367282     -6.811023828514404    -0.1563851149559012      1.509080347621573   -0.08309889142654189    -0.5063474784942394     -1.916832267159649     0.5407799496439279      1.528703189392817    -0.3288360549480647    -0.2744963977687636     -18.86538500194821     0.1484355388197094      1.867398227696617    -0.0556488760024214       1.97029428610251     -7.192107425440216    -0.4508354820785889     0.1342419283620955     0.9010179200136184     -1.202924325648357     -1.174556392291442      6.199160943301402     -7.192761797132046      1.112493651694865     -11.56852700441887       16.6610526581012]';

%Logging
diary OFF
diaryFileName = ['estimationOutput/bsEst_' num2str(id) '.txt'];
delete(diaryFileName);
diary(diaryFileName);
%Loading regional moments and creating another estimator
load('../CleanData/moments_regions.mat');
if max(abs(imag(bsMoments(:,id)))) > 1e-6
    error('Imaginary moment!');
else %Sanitize input
    bsMoments(:,id) = real(bsMoments(:,id));
    zeroTargets = bsMoments(:,id)==0;
    bsMoments(zeroTargets,id) = 1e-4;
    for r = 1:6
        fn = fieldnames(bsDataSummaries{id}{r});
        for fid = 1:length(fn)
            f = fn{fid};
            bsDataSummaries{id}{r}.(f) = real(...
                bsDataSummaries{id}{r}.(f));
            if strcmp(f(1:2),'fs') || strcmp(f,'informal') ...
                    ||  strcmp(f(1:2),'tr') ...
                    ||  strcmp(f(1:2),'th') ...
                    ||  strcmp(f(1:2),'la')
                bsDataSummaries{id}{r}.(f) = max(...
                    bsDataSummaries{id}{r}.(f), 1e-4);
            end
        end
    end
end
ME = ModelEstimator(bsMoments(:,id),V,bsDataSummaries{id});
%Set starting point
ME.centerX = x_point_est;
ME.rangeSampleX = 0;
ME.initialSampleDraws = 1;
%Use derivative-based method
ME.algOption = ModelEstimator.alg_fminunc;
%Estimate and save results
results = ME.run();
if isempty(results) %An error with initial point;
    %try again with other starting points
    ME.rangeSampleX = 0.1;
    ME.initialSampleDraws = 20;
    results = ME.run();
end

end

