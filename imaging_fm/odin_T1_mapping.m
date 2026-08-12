function [T1_map, psir_im_all] = odin_T1_mapping(data, TIs)
% t1_maps = odin_T1_mapping(data, TIs)

RO = size(data, 1);
E1 = size(data, 2);
SLC = size(data, 3);

N = size(data, 4);

T1_map = zeros(RO, E1, SLC);

[a, ind] = sort(TIs);
TIs = a;

data = data(:,:,:,ind);

psir_im_all = data;
psir_im_all(:,:,1,:) = 0;
psir_im_all(:,:,end,:) = 0;

for slc=2:SLC-1
    disp(['processing slice ' num2str(slc) ' out of ' num2str(SLC)])
    im = data(:,:, slc, :);
    im = squeeze(im);
    
    base_im = im(:,:, end);
   
    phs_im = base_im ./ abs(base_im);
    
    phs_removed_im = im .* repmat(conj(phs_im), [1 1 size(im, 3)]);
    psir_im = real(phs_removed_im);
    psir_im_all(:,:,slc,:) = psir_im;
    
    % get initial parameters
    
    A = psir_im(:,:,end);
    B = 2*A;
      
    TIS = zeros(RO, E1);
    
    for e1=1:E1
        for ro=1:RO
            T1S(ro, e1) = interp1(squeeze(psir_im(ro, e1, :)), TIs, 0, 'linear');
            if isnan(T1S(ro, e1)) | T1S(ro, e1)<=TIs(1) | T1S(ro, e1)>1000
                a = squeeze(psir_im(ro, e1, :));
                [v, ind] = min(abs(a));
                T1S(ro, e1) = TIs(ind);
            end
        end
    end
    T1S = T1S / log(2);
       
    % [S, J] = phase_sensitive_inversion_recovery_3param(params, TI)
    
    T1 = T1S;
    res = zeros(RO, E1, 3);
    
    options = optimset('MaxFunEvals', 4000, 'MaxIter', 4000, 'TolFun',1e-4,'Display','off');
    % options = optimset('MaxFunEvals', 2000, 'MaxIter', 200, 'TolFun',1e-4,'PlotFcns',@optimplotx);
    
    for e1=1:E1
        %disp(['processing - ' num2str(e1)])
        for ro=1:RO
            if (T1S(ro, e1) > 2)
                y_exp = squeeze(psir_im(ro, e1, :));
    
                model = @(params, TIs) params(1) - params(2) * exp(-TIs(:)./params(3));
                obj_fun = @(params) sum((y_exp - model(params, TIs)).^2);
    
                params0 = [A(ro, e1) B(ro, e1), T1S(ro, e1)];
                [opt_params, fval] = fminsearch(obj_fun, params0, options);
    
                res(ro, e1, :) = opt_params;
                T1(ro, e1) = opt_params(3) * (opt_params(2)/opt_params(1) - 1);
                if (T1(ro, e1)<1) | (T1(ro, e1)>2000)
                    T1(ro, e1) = 0;
                end
    
                if 0
                    y_hat = model(opt_params, TIs);
        
                    figure;
                    hold on
                    plot(y_exp);
                    plot(y_hat, 'r')
                    hold off
                    box on
                    grid on
                end
            end        
        end
    end

    T1_map(:,:,slc) = T1;
end
