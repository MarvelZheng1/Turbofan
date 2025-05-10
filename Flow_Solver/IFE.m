classdef IFE
    properties
    end
    methods(Static)
        function T_ratio = T_ratio_from_P_ratio(P_ratio, gamma)
            T_ratio = (P_ratio)^((gamma-1)/gamma)
        end
    end
end
