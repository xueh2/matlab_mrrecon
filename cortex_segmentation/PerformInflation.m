
function PerformInflation(Input, Output, NumberOfIterations, RelaxationFactor, ThresRatio, target_stats, MaxofIterations)
% perform inflation until the second order statistics are less or close to the input values
% four second order statistics can be used, [ICI, ECI, GLN, MLN];
% -1 in the input values means the corresponding second order statistics are not in use

% Input = vtkfile_Medium;
% Output = vtkfile_Medium_inflated;
% NumberOfIterations = 5;
% RelaxationFactor = 0.5;
% ThresRatio = 0.01;
% target_stats = [ICI_Simple ECI_Simple MLN_Simple];

num = length(target_stats);
TargetStats = [];
for i=1:num
    if ( target_stats(i)==-1 )
        continue;
    end
    TargetStats = [TargetStats target_stats(i)];
end

Inflation_TriangleMesh3(Input, Output, NumberOfIterations, RelaxationFactor);
[ICI_inflated, ECI_inflated, GLN_inflated, MLN_inflated, meanC, gaussC] = Compute_SecondOrderStatistics_OutlierRemoval(Output);
% [ICI_inflated, ECI_inflated, GLN_inflated, MLN_inflated, meanC, gaussC] = SecondOrder_Statistics_vtk(Output);

source_stats = [ICI_inflated ECI_inflated GLN_inflated MLN_inflated];
SourceStats = [];
for i=1:num
    if ( target_stats(i)==-1 )
        continue;
    end
    SourceStats = [SourceStats source_stats(i)];
end

index = 1;
while ( ~StopInflation(ThresRatio, TargetStats, SourceStats) )
    disp('----------------------------------------------------------------');
    if ( index*NumberOfIterations >= MaxofIterations )
        break;
    end
    
    Inflation_TriangleMesh3(Output, Output, NumberOfIterations, RelaxationFactor);
    
    [ptIDs, meanC, gaussC] = Curvature_PolyData(Output);

    [ICI_inflated, ECI_inflated, GLN_inflated, MLN_inflated, meanC, gaussC] = Compute_SecondOrderStatistics_OutlierRemoval(Output);

%     [ICI_inflated, ECI_inflated, GLN_inflated, MLN_inflated, meanC, gaussC] = SecondOrder_Statistics_vtk(Output);
    
    source_stats = [ICI_inflated ECI_inflated GLN_inflated MLN_inflated];
    SourceStats = [];
    for i=1:num
        if ( target_stats(i)==-1 )
            continue;
        end
        SourceStats = [SourceStats source_stats(i)];
    end

    index = index + 1;
    errors = source_stats - target_stats
    disp([num2str(index*NumberOfIterations) ' iteration of the inflation has been performed ... ']);
end
return