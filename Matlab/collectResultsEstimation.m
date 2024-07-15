
N_max=200;
N_display = 5;
gamma = nan(N_max,1);
ranges = {101:150;151:175;176:200};
header = {'Main gamma';'Low subs';'High subs'};
for r = 1:3
results = inf(N_max,1);
countevals = zeros(N_max,1);
for i = ranges{r}
    filename = ['./estimationOutput/estimation_' num2str(i) '.txt'];
    if length(dir(filename)) == 0
        continue;
    end
    fhandle = fopen(filename,'rt');
    thisline = fgetl(fhandle);
    gamma(i) = str2double(thisline((end-5):end));
    while true
        thisline = fgetl(fhandle);
        if ~ischar(thisline); break; end  %end of file
        
        lineStartCheck = 'Loss = '; 
        beginning = thisline(1:min(...
            length(thisline),length(lineStartCheck)));
        if strcmp(beginning,lineStartCheck)
            countevals(i) = countevals(i) + 1;
            lastChar = min(length(lineStartCheck)+6,length(thisline));
            val = str2double(thisline(...
                (length(lineStartCheck)+1):lastChar));
            if val < results(i)
                results(i) = val;
            end
        end
    end
    fclose(fhandle);
end
[sorted_results,order] = sort(results);
res = [sorted_results order countevals(order)];
disp([header{r} '; ' num2str(N_display) ' best results']);
disp('Loss fun / Id / Fun evals');
disp(res(1:N_display,:));
end