function [x,fv] = binarySearch(fun,low,f_low,high,f_high,tolX,tolFun,verbose,level)
    minLevelSimpleBisection = 20;
    if nargin < 8
        verbose = false;
    end
    if nargin < 9
        level = 0;
    end
    if verbose
        disp(['Level ' num2str(level) ': entering [' num2str([low high],'%.18f') ...
            '], f: [' num2str([f_low f_high]) ']']);
    end
    if f_low*f_high > 0
        error('binarySearch: invalid initial range.');
    elseif abs(f_low) < tolFun
        x = low;
        fv = f_low;
        if verbose
            disp('Exit left');
        end
        return;
    elseif (abs(f_high) < tolFun) || (high-low < tolX)
        x = high;
        fv = f_high;
        if verbose
            disp('Exit right');
        end
        return;
    end
    
    if level >= minLevelSimpleBisection
        mid = (low+high)/2;
    else
        mid = (low*abs(f_high)+high*abs(f_low))/(abs(f_high)+abs(f_low));
    end
    f_mid = fun(mid);
    if f_low*f_mid < 0
        [x,fv] = binarySearch(fun,low,f_low,mid,f_mid,tolX,tolFun,verbose,level+1);
    else
        [x,fv] = binarySearch(fun,mid,f_mid,high,f_high,tolX,tolFun,verbose,level+1);
    end
end