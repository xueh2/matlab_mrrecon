function [T1, PD] = t1fit(TI, data, method);
% function [T1, PD] = t1fit(TI, data, method);
%
% function to estimate T1 from IR or SR data given inversion or saturation recovery times (TI)
% using levenberg-marquardt curve fit.
%
% method is either:
%     PSIR (default)
%     magnitude
%     MOLLI
%     SR

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  NIH NHLBI                          *
%     ***************************************

params0=[200, 500]; % initial value

if nargin < 3
    method = 'PSIR'; % default
end

switch method
    case {'PSIR','psir'}
        options = optimset('Algorithm','levenberg-marquardt','Jacobian','on','TolFun',1e-8,'Display','off');
        fit = lsqcurvefit('phase_sensitive_inversion_recovery', params0, TI, data,[],[],options);
    case {'MAGIR','magnitude', 'magn'}
        options = optimset('Algorithm','levenberg-marquardt','TolFun',1e-8,'Display','off');
        fit = lsqcurvefit('magnitude_inversion_recovery', params0, TI, abs(data),[],[],options);
    case {'MOLLI','molli'}
        data = abs(data);
        [TI,index]=sort(TI);
        data = data(index);        
        options = optimset('Algorithm','levenberg-marquardt','Jacobian','on','TolFun',1e-8,'Display','off');
        for i = 1:length(TI)-1
            tmpdata = data;
            tmpdata(1:i)= -tmpdata(1:i);
            [fit(:,i), resnorm(i)] = lsqcurvefit('phase_sensitive_inversion_recovery', params0, TI, tmpdata,[],[],options);
        end
        index = min(find(resnorm==min(resnorm)));
        fit = fit(:,index);
    case {'SR', 'sr'}
        options = optimset('Algorithm','levenberg-marquardt','TolFun',1e-8,'Display','off');
        fit = lsqcurvefit('saturation_recovery', params0, TI, abs(data),[],[],options);
end

T1 =  fit(2); % estimate of T1

if nargout ==2
    PD = fit(1);
end





