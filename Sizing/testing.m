clear;clc;clf























































% Turbofan = struct( ...
%     "Info", [],...
%     "Compressor", struct( ...
%         "HP", repmat(struct( ...
%             "T01", [], ...
%             "P01", []), ...
%             num_HP_compressors, 1), ...
%         "LP", repmat(struct( ...
%             "T01", [], ...
%             "P01", []), ...
%             num_LP_compressors, 1) ...
%         ), ...
%     "Combustor", [], ...
%     "Turbine", struct( ...
%         "HP", repmat(struct( ...
%             "T01", [], ...
%             "P01", []), ...
%             num_HP_turbines, 1), ...
%         "LP", repmat(struct( ...
%             "T01", [], ...
%             "P01", []), ...
%             num_LP_turbines, 1) ...
%         ) ...
%     );