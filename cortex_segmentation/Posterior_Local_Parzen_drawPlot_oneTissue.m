
function [post_local, a, plot_a] = Posterior_Local_Parzen_drawPlot_oneTissue(post, mix, x, minimalP, ...
    typeofselection, maxN, Hmin, Hmax, tissueK, H, post_local, a, plot_a)
% use the parzen local probability density estimation to refine the posterior estimation
% plot_a is the class probability from the parzon method

ndata = size(x, 1);
numofclass = size(post, 2);

a = zeros(size(mix.priors), 'single');

tissue_c = zeros(numofclass, 1);

minX = min(x);
maxX = max(x);
testX = [minX:maxX]';
% plot_a = zeros([maxX-minX+1, 1], 'single');

for i = 1:ndata
    label = find(post(i,:) == max(post(i,:)));
    tissue_c(label) = tissue_c(label)+1;
end

for k = tissueK:tissueK
    
    disp(['estimate the ' num2str(k) ' tissue ... '])
    
    tissue_count = tissue_c(k);
    tissue_I = zeros([tissue_count, 1], 'single');
    tissue_label = zeros([tissue_count, 1], 'single');

    index_tissue = 1;
    for i = 1:ndata
        label = find(post(i,:) == max(post(i,:)));
        if ( label == k ) 
            tissue_I(index_tissue) = x(i);
            tissue_label(index_tissue) = i;
            index_tissue = index_tissue + 1;
        end
    end

%     if ( (k==1) | (k==5) | (k==4)) % favor the gm and low intensity wm
%         typeofselection = ['probability'];
%     else
%         typeofselection = ['uniform'];
%     end
    
    [selected_tissue, selected_tissuelabel] = SelectSamples(tissue_I, tissue_label, post(:,k), minimalP, typeofselection, maxN);

    disp('estimate parzen optimal BW ... ');
    minI = min(selected_tissue)
    maxI = max(selected_tissue)
    selected_tissue2 = (selected_tissue - minI) ./ (maxI - minI);


%     H_min=0.01;
%     H_max=0.1;
    if ( H ~= -1 ) 
        K=20;
        cross=10;
        h_opt_parz = parzen(selected_tissue2,[Hmin,Hmax],K,cross);
    else
        h_opt_parz = H;
    end
    
    parzen_P = parzenprob_mex(single(selected_tissue2), single((tissue_I-minI) ./ (maxI-minI)), h_opt_parz);
    for i = 1:tissue_count
        a(tissue_label(i), k) = parzen_P(i);
    end
    
    % for plot purpose
    plot_a(:,k) = parzenprob_mex(single(selected_tissue2), single((testX-minI) ./ (maxI-minI)), h_opt_parz);;
end

post_local = mix.priors.*a;
s = sum(post_local, 2); % p(x)
if any(s==0)
   warning('Some zero posterior probabilities')
   % Set any zeros to one before dividing
   zero_rows = find(s==0);
   s = s + (s==0);
   post_local(zero_rows, :) = 1/mix.ncentres;
end
post_local = post_local./(s*ones(1, mix.ncentres));

return