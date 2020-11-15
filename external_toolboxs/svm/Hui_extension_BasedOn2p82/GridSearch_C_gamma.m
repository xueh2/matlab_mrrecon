
function [goodC, goodGamma] = GridSearch_C_gamma(label, data, foldnum,...
    C_low, C_high, C_step, ...
    gamma_low, gamma_high, gamma_step, ...
    finesearch_flag, fine_range, fine_step)

% perform a grid search of the C and gamma for svm classification
% exponentially growing sequences

allC = C_low:C_step:C_high;
numC = length(allC);
allGamma = gamma_low:gamma_step:gamma_high;
numGamma = length(allGamma);

accuracyNum = 0;
goodC = -1;
goodGamma = -1;

disp('loose grid search ...');

for i = 1:numC
    for j = 1:numGamma
        
        disp(['current test:' num2str(j + (i-1)*numGamma) '/' num2str(numC*numGamma)]);
        
        currentC = 2^allC(i);
        currentGamma = 2^allGamma(j);
        libsvm_options = ['-c ' num2str(currentC) ' -g ' num2str(currentGamma)];
        
        anum = CrossValidation(label, data, foldnum, libsvm_options);
        
        if ( anum > accuracyNum )
            accuracyNum = anum;
            goodC = allC(i);
            goodGamma = allGamma(j);            
        end
        disp('% ============================================== %');
        disp('');
    end
end

disp(['Loose search final: ' num2str(accuracyNum/size(label,1)) ',' num2str(goodC) ',' num2str(goodGamma) ]);

if ( finesearch_flag )
    disp('find grid search ...');

    allC = goodC-fine_range:fine_step:goodC+fine_range;
    numC = length(allC);
    allGamma = goodGamma-fine_range:fine_step:goodGamma+fine_range;;
    numGamma = length(allGamma);

    accuracyNum = 0;
    goodC = -1;
    goodGamma = -1;

    disp('loose grid search ...');

    for i = 1:numC
        for j = 1:numGamma

            disp(['current test:' num2str(j + (i-1)*numGamma) '/' num2str(numC*numGamma)]);

            currentC = 2^allC(i);
            currentGamma = 2^allGamma(j);
            libsvm_options = ['-c ' num2str(currentC) ' -g ' num2str(currentGamma)];

            anum = CrossValidation(label, data, foldnum, libsvm_options);

            if ( anum > accuracyNum )
                accuracyNum = anum;
                goodC = allC(i);
                goodGamma = allGamma(j);            
            end
            disp('% ============================================== %');
            disp('');
        end
    end
    disp(['fine search final: ' num2str(accuracyNum/size(label,1)) ',' num2str(goodC) ',' num2str(goodGamma) ]);
end
return