% function res = FullWavelet(filterType, filterSize, wavScale)
% % res = FullWavelet(Filtertype, filterSize, wavScale)
% %
% % implements a wavelet operator
% %
% % (c) Michael Lustig 2007
% 
% res.adjoint = 0;
% res.wavScale = wavScale;
% res = class(res,'FullWavelet');
classdef FullWavelet
   properties
      adjoint = 0;
      wavScale = 0;
   end
   
    methods
      function td = FullWavelet(adjoint_in,wavScale_in)
         if nargin > 0
            td.adjoint = adjoint_in;
            td.wavScale = wavScale_in;
         end
      end
      
      function res = mtimes(a,b)
          res = a;
          res.adjoint = a.adjoint*b;
          res.wavScale = a.wavScale*b;
      end
      
      function disp(td)
         disp('sasdf');
      end % disp
   end
   
end