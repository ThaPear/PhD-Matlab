classdef ABCDMatrix
    properties
        A
        B
        C
        D
    end
    methods
        function obj = ABCDMatrix(A, B, C, D)
            switch(nargin)
                case 1 % The first argument is an ABCD matrix.
                    obj.A = A.A;
                    obj.B = A.B;
                    obj.C = A.C;
                    obj.D = A.D;
                case 4 % Simply 4 specified paremeters.
                    obj.A = A;
                    obj.B = B;
                    obj.C = C;
                    obj.D = D;
                otherwise
                    error('%s::ABCDMatrix:\n\tInvalid number of constructor arguments.', mfilename);
            end
        end
        
        function [result] = mul(obj, ABCD)
            % Places the given matrix in series with this object.
            
            %        |------|    |------|
            %     -- |      | -- |      | -- 
            % -->    | obj  |    | ABCD |    -->
            %     -- |      | -- |      | -- 
            %        |------|    |------|
            
            % |--------------|
            % | AA+BC  AB+BD |
            % | CA+DC  CB+DD |
            % |--------------|
            
            Anew = obj.A .* ABCD.A + obj.B .* ABCD.C;
            Bnew = obj.A .* ABCD.B + obj.B .* ABCD.D;
            Cnew = obj.C .* ABCD.A + obj.D .* ABCD.C;
            Dnew = obj.C .* ABCD.B + obj.D .* ABCD.D;
            result = ABCDMatrix(Anew, Bnew, Cnew, Dnew);
        end
        
%         function [result] = parallel(obj, ABCD)
%             % Places the given matrix in parallel with this object.
%             
%             %              |------|
%             %          /-- |      | --
%             %         /    | obj  |    -->
%             %        / /-- |      | -- 
%             %     --/ /    |------|
%             % -->    X  
%             %     --\ \    |------|
%             %        \ \-- |      | --
%             %         \    | ABCD |    -->
%             %          \-- |      | -- 
%             %              |------|
%             
%             
%             % |--------------|
%             % | AA+BC  AB+BD |
%             % | CA+DC  CB+DD |
%             % |--------------|
%             
% %             Anew = obj.A .* ABCD.A + obj.B .* ABCD.C;
% %             Bnew = obj.A .* ABCD.B + obj.B .* ABCD.D;
% %             Cnew = obj.C .* ABCD.A + obj.D .* ABCD.C;
% %             Dnew = obj.C .* ABCD.B + obj.D .* ABCD.D;
% %             result = ABCDMatrix(Anew, Bnew, Cnew, Dnew);
%         end
    end
end