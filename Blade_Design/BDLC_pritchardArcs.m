% Generates the LE, TE, UGT circular arcs
function blade = BDLC_pritchardArcs(x, y, betas, R_LE, R_TE, res)
    [blade.LEx, blade.LEy] = arc(x(3) + R_LE * cos(pi/2 - betas(3)), y(3) - R_LE * sin(pi/2 - betas(3)), R_LE, x(3), x(4), y(3), y(4), 1, res);
    [blade.TEx, blade.TEy] = arc(x(1) + R_TE * cos(pi/2 - betas(1)), y(1) - R_TE * sin(pi/2 - betas(1)), R_TE, x(5), x(1), y(5), y(1), 1, res);
    % [blade.UGx, blade.UGy] = arc(x(6), y(6), R_new, x(1), x(2), y(1), y(2), 10000, res);
    blade.x = x;
    blade.y = y;
    blade.betas = betas;
end

% Helper function for generating arc coordinates
function [x_arc, y_arc] = arc(x,y,r, x1, x2, y1, y2, multiplier, res)
    start = atan2(y1-y,x1-x);
    stop = atan2(y2-y,x2-x);
    if start > stop
        stop = stop + 2 * pi;
    end
    theta = start:pi/(100*multiplier*res):stop;
    x_arc = r * cos(theta) + x;
    y_arc = r * sin(theta) + y;
end