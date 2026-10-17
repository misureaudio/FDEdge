classdef test_interp2d < matlab.unittest.TestCase
    % M5: bilinear interpolation (move of the PoC's idx_w + bilinear
    % combination).
    %
    % WEIGHT CONVENTION (verbatim PoC behavior): the combination
    % w*G(i1) + (1-w)*G(i2) with w = (val - v(i1))/(v(i2) - v(i1))
    % evaluates the bilinear surface at the query point REFLECTED through
    % the cell center: exact at cell centers and at the range edges; at
    % an exact interior node (i,j) it returns G(i+1, j+1). The M8 golden
    % replay requires this verbatim behavior (see the interp2d docstring).
    methods (Test)
        function exactAtMidpoints(t)
            % For a grid bilinear in the coordinates, exact at every
            % cell center (the reflection is a no-op there).
            axisRows = (0.20:0.10:0.60)';
            axisCols = (0.55:0.05:0.85)';
            G = axisRows*axisCols' + 3*axisRows + 7*axisCols' + 2;   % 5x7
            for i = 1:4
                for j = 1:6
                    qRow = (axisRows(i) + axisRows(i+1))/2;
                    qCol = (axisCols(j) + axisCols(j+1))/2;
                    d = interp2d(G, axisRows, axisCols, qRow, qCol);
                    dExact = qRow*qCol + 3*qRow + 7*qCol + 2;
                    t.verifyLessThan(abs(d - dExact), 1e-14);
                end
            end
        end

        function nodeConvention(t)
            % At an exact interior node (i,j) the PoC combination returns
            % the upper-right neighbor G(i+1, j+1); at the last row/col
            % (i1 == i2) it is exact. Bit-exact (w = 0 or i1 == i2).
            axisRows = [0.2 0.4 0.6]';
            axisCols = [0.55 0.7 0.85]';
            G = 10*axisRows + axisCols';   % implicit expansion: 3x3
            t.verifyEqual(interp2d(G, axisRows, axisCols, 0.2, 0.55), G(2,2));
            t.verifyEqual(interp2d(G, axisRows, axisCols, 0.4, 0.7), G(3,3));
            t.verifyEqual(interp2d(G, axisRows, axisCols, 0.2, 0.85), G(2,3));
            t.verifyEqual(interp2d(G, axisRows, axisCols, 0.6, 0.85), G(3,3));
        end

        function outOfRangeClipped(t)
            % Out-of-range queries are clipped to the range edges (no
            % extrapolation): the upper corner is exact; the lower corner
            % maps to the (2,2) cell by the PoC weight convention.
            axisRows = (0.20:0.2:0.60)';
            axisCols = (0.5:0.1:0.9)';
            G = 10*axisRows + axisCols';   % implicit expansion: 3x5
            t.verifyEqual(interp2d(G, axisRows, axisCols, 1.0, 0.9), G(3,5));
            t.verifyEqual(interp2d(G, axisRows, axisCols, 0.0, 0.4), G(2,2));
            t.verifyEqual(interp2d(G, axisRows, axisCols, 0.4, -0.2), G(3,2));
            % never extrapolates beyond the grid's value range
            d = interp2d(G, axisRows, axisCols, -5, 5);
            t.verifyGreaterThanOrEqual(d, min(G(:)));
            t.verifyLessThanOrEqual(d, max(G(:)));
        end

        function bothOrientations(t)
            % Works in the dgrid orientation (rows = alphas, cols =
            % beta_ref) and the dbeta orientation (rows = xes, cols =
            % betas) - the axes are explicit, so the same call covers
            % both. Cell-center queries: exact for a grid linear in the
            % coordinates.
            alphas = (0.55:0.05:0.85)';
            beta_ref = [0.55 0.65 0.75 0.85]';
            G1 = 10*alphas + beta_ref';    % implicit expansion: 7x4
            qRow = (alphas(2) + alphas(3))/2;       % 0.625, cell center
            qCol = (beta_ref(2) + beta_ref(3))/2;   % 0.70, cell center
            t.verifyLessThan(abs(interp2d(G1, alphas, beta_ref, qRow, qCol) - ...
                (10*qRow + qCol)), 1e-14);
            xes = (0.20:0.10:0.60)';
            betas = alphas;
            G2 = 100*xes + betas';         % implicit expansion: 5x7
            qRow2 = (xes(2) + xes(3))/2;        % 0.35
            qCol2 = (betas(3) + betas(4))/2;    % 0.70
            t.verifyLessThan(abs(interp2d(G2, xes, betas, qRow2, qCol2) - ...
                (100*qRow2 + qCol2)), 1e-14);
        end
    end
end
