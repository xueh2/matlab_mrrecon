%Class to simplify the use of linear operators
%If a linear operator and its transpose is defined by a function, give them
%to this class so that you can use them on other linear operators or on the
%signals simply using multiplications and transposes.

%BEWARE: this is only efficient to handle the operators, it may slow
%down computation

%BEWARE (2): it uses the transpose operator, not the inverse operator! it
%may be the same in some cases, but it makes a big difference most of the
%time.

%Example (if waveletTransform and inverseWaveletTransform were defined
%functions):
%Phi=SimpleOperator(@(x) waveletTransform(x),@(x) inverseWaveletTransform(x))
%x=cos(1:128)';
%y=Phi*x
%x=Phi'*y;
%Identity=Phi'*Phi
%x=Identity*x

classdef SimpleOperator
    %SimpleOperator easy to handle operators
    %   operator = SimpleOperator(directOperator,transposeOperator)
    %   containing the coefficients of descending powers of x.

    % The following properties can be set only by class methods
    properties (SetAccess = private)
        dirFct
        transFct
    end

    methods

        function op=SimpleOperator(directOperator,transposeOperator)
            %SimpleOperator easy to handle operators
            %   operator = SimpleOperator(directOperator,transposeOperator)
            %   containing the coefficients of descending powers of x.
            op.dirFct=directOperator;
            op.transFct=transposeOperator;
        end

        function transOp = ctranspose(op)
            %transpose of the operator
            transOp = SimpleOperator(op.transFct,op.dirFct);
        end

        function y = mtimes(a,b)
            %operator on a signal or concatenation with another operator
            if isa(b,'SimpleOperator')
                if isa(a,'double')
                    y=SimpleOperator(@(z) a*b.dirFct(z),@(z) a*b.transFct(z));
                else
                    if isa(b,'SimpleOperator')
                        y=SimpleOperator(@(z) a.dirFct(b.dirFct(z)),@(z) b.transFct(a.transFct(z)));
                    else
                        fprintf('No implementation for this type of product with SimpleOperator')
                    end
                end
            else
                y=a.dirFct(b);
            end

        end

        function y = plus(a,b)
            %sum of operators
            if isa(b,'SimpleOperator')
                y=SimpleOperator(@(z) a.dirFct(z)+b.dirFct(z),@(z) b.transFct(z)+a.transFct(z));
            else
                fprintf('No implementation for this type of addition with SimpleOperator')
            end
        end
        
        function y = minus(a,b)
            %difference of operators
            if isa(b,'SimpleOperator')
                y=SimpleOperator(@(z) a.dirFct(z)-b.dirFct(z),@(z) a.transFct(z)-b.transFct(z));
            else
                fprintf('No implementation for this type of addition with SimpleOperator')
            end
        end
        
        function y = uminus(a)
            %difference of operators
            if isa(b,'SimpleOperator')
                y=SimpleOperator(@(z) -a.dirFct(z),@(z) -a.transFct(z));
            else
                fprintf('No implementation for this type of addition with SimpleOperator')
            end
        end
        
        
    end
end
