function [results, clusters] = Gergo_psth_sorter(g, psth_spx)
% Sorts PSTH based on response onset and response direction (inh or exc)

for ii = 1:size(psth_spx,1)
    datasegment = psth_spx(ii,g.roi);
    if all(datasegment == datasegment(1))
        results.onsetIdx(ii,1) = -1;
        results.offsetIdx(ii,1) = length(datasegment); 
        if datasegment(1)>=g.exctestvalue
            results.respDir(ii,1) = -1;
        elseif datasegment(1)<=-g.inhtestvalue
             results.respDir(ii,1) = 1;
        else
             results.respDir(ii,1) = 0;
        end
    elseif max(datasegment)>=g.exctestvalue
        results.respDir(ii,1) = -1;
        results.onsetIdx(ii,1) = find(datasegment>=g.exctestvalue,1,'first'); % begining of the response
        results.offsetIdx(ii,1) = find(datasegment>=g.exctestvalue,1,'last'); % end of response
        responseMean = mean(datasegment(results.onsetIdx(ii,1):results.offsetIdx(ii,1)));
        k = 0;
        while k == 0
            if responseMean < g.exctestvalue/2
                datasegment(results.offsetIdx(ii,1):end) = NaN;
                results.offsetIdx(ii,1) = find(datasegment>=g.exctestvalue,1,'last'); % end of response
                responseMean = mean(datasegment(results.onsetIdx(ii,1):results.offsetIdx(ii,1)));
            else 
                k = 1;
            end
        end
    elseif max(datasegment)<g.exctestvalue
        datasegment_inverse = -smoothdata(datasegment);
        if max(datasegment_inverse)>=g.inhtestvalue
            results.respDir(ii,1) = 1;
            results.onsetIdx(ii,1) = find(datasegment_inverse>=g.inhtestvalue,1,'first'); % begining of the response
            results.offsetIdx(ii,1) = find(datasegment_inverse>=g.inhtestvalue,1,'last'); % end of response
            responseMean = mean(datasegment_inverse(results.onsetIdx(ii,1):results.offsetIdx(ii,1)));
            k = 0;
            while k == 0
                if responseMean < g.inhtestvalue/2
                    datasegment_inverse(results.offsetIdx(ii,1):end) = NaN;
                    results.offsetIdx(ii,1) = find(datasegment_inverse>=g.inhtestvalue,1,'last'); % end of response
                    responseMean = mean(datasegment_inverse(results.onsetIdx(ii,1):results.offsetIdx(ii,1)));
                else 
                    k = 1;
                end
            end
        else
            results.respDir(ii,1) = 0;
            results.onsetIdx(ii,1) = -1;
            results.offsetIdx(ii,1) = length(datasegment);
        end
    end
end


edges = 0:50:g.test_time*1000;
onsetTime_d = discretize(results.onsetIdx*g.bin_time*1000, edges); % round to 10 ms values.
offsetTime_d = discretize(results.offsetIdx*g.bin_time*1000, edges); % round to 10 ms values.
durationTime_d = discretize((results.offsetIdx-results.onsetIdx)*g.bin_time*1000, edges); % round to 10 ms values.

clusters(results.respDir == -1 & onsetTime_d == 1 & durationTime_d == 1,1) = 1;
clusters(results.respDir == -1 & onsetTime_d == 1 & durationTime_d > 1,1) = 2;
clusters(results.respDir == -1 & onsetTime_d > 1 ,1) = 2;
clusters(results.respDir == 1 & onsetTime_d == 1 & durationTime_d == 1,1) = 3;
clusters(results.respDir == 1 & onsetTime_d == 1 & durationTime_d > 1,1) = 3;
clusters(results.respDir == 1 & onsetTime_d > 1 ,1) = 3;
clusters(results.respDir == 1 & isnan(onsetTime_d)) = 5;   
clusters(results.respDir == 0) = 5;