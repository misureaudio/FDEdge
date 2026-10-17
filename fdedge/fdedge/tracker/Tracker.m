classdef Tracker
    % Per-frame FD edge tracking pipeline (M7 orchestrator).
    %
    % Per-profile signal contract (update(signal, xguess)):
    %   Profile 'pocB'  : 1025 points on x = (0:1024)'/1024   (one-sided edge)
    %   Profile 'beta'  : 1025 points on x = (0:1024)'/1024   (one-sided edge)
    %   Profile 'riesz' : 1024 points on xr = (0:1023)'/1024  (periodic, symmetric edge)
    %   Profile 'hybrid': 1025 points on x = (0:1024)'/1024   (Riesz bank computed
    %                      internally on signal(1:1024); GL on all 1025 points)
    %   Part A is NOT a tracking profile: use fdedge.analyzeMatchedOrder on the
    %   2049-point non-normalized grid t = (0:2048)'.
    %
    % update(signal, xguess) runs the per-profile per-frame block line-for-line
    % from the PoC (FD_EdgeDetector_PoC/edge_fd_*.m) with local-function calls
    % replaced by fdedge.* calls (same operation order - standing rule D14):
    %   poc    = edge_fd_poc.m Part B loop body (halfCross at level
    %            0.5*quantile(Rn,0.90) + scalar delta from calib.scalarHalfCross)
    %   beta   = edge_fd_beta_track.m per-frame block (GL bank -> c1 window ->
    %            exponentBank -> slopePeak at the estimated order -> debias)
    %   riesz  = edge_fd_riesz_track.m per-frame block (Riesz bank -> 11 probes
    %            rieszGfit per alpha -> exponentBank -> cuspMin -> debias)
    %   hybrid = edge_fd_hybrid_track.m two-pass block (estExponent at xguess ->
    %            cuspMin at ia1 -> estExponent at the fresh position -> debias)
    %
    % DEBIAS NOTES (byte-identity with the golden replays):
    %   - beta keeps the PoC's interp2(beta_ref, alphas, dgrid, b, b) call
    %     VERBATIM - interp2 is a built-in, not a local function, so the
    %     standing rule does not replace it. Convention verified on R2024a
    %     (probe): the FIRST grid vector pairs with the SECOND dimension of
    %     Z (meshgrid convention), so the PoC's call is dimensionally and
    %     semantically correct (row axis = alphas, column axis = beta_ref).
    %     Do NOT "fix" the axis order - the M9 byte gate depends on the
    %     PoC's behavior. (interp2d - the manual idx_w move - differs from
    %     interp2 at exact interior nodes (returns the upper-right neighbor);
    %     continuous bhat_raw queries never land on nodes, so formatted
    %     output is unaffected either way.)
    %   - riesz/hybrid use interp2d (M5) with the verbatim PoC weight
    %     convention (reflected through the cell center - see interp2d).
    %
    % KF: predict/update from frame 2 (the PoC pass-2 structure). poc: 2-state;
    % others: 4-state with the PoC's explicit P0 built from the raw variances:
    %   2-state: P0 = [R 0; 0 R]
    %   4-state: P0 = [Rx 0 0 0; 0 Rx*1e2 0 0; 0 0 Rb 0; 0 0 0 Rb*1e4]
    % The state is initialized from the first two raw estimates exactly as the
    % PoC (X = [raw1; raw2-raw1; rawb1; rawb2-rawb1]), so frame-2 output equals
    % the PoC's kf(2) and every later frame matches the PoC loop.
    % The replay driver sets kfR to the PoC's raw variances (var(raw - truth),
    % which the driver knows - computed in a raw pre-pass) BEFORE the loop for
    % a PoC-identical KF; the default kfR is a placeholder.
    %
    % Returns struct('xhat',...,'bhat',...,'vhat',...,'ohat',...,'xraw',...,'braw',...)
    % (vhat/ohat NaN at frame 1; bhat/braw NaN for the position-only poc
    % profile; xraw/braw are the raw (debiased) estimates - the driver manages
    % xguess from them with the PoC's const-velocity formula).
    %
    % CALL FORM (MATLAB classdef value semantics):  [st, trk] = trk.update(signal, xguess)
    % - update returns the UPDATED object (the same pattern as KalmanCV);
    %   without it the frame counter / KF state would not persist.
    properties
        cfg; calib; grid; kf; n; kfR; xraw1; braw1; xrawPrev; brawPrev
    end
    methods
        function obj = Tracker(cfg, calib)
            obj.cfg = cfg;
            obj.calib = calib;
            obj.n = 0;
            obj.kf = [];
            obj.xraw1 = NaN; obj.braw1 = NaN;
            obj.xrawPrev = NaN; obj.brawPrev = NaN;
            M = cfg.M;
            switch lower(cfg.Profile)
                case 'riesz'
                    obj.grid = (0:M-1)'/M;      % periodic (Riesz PoC convention)
                    obj.kfR = [1e-6 0; 0 1e-2]; % placeholder (see class header)
                case 'poc'
                    obj.grid = (0:M)'/M;
                    obj.kfR = [1e-6];           % placeholder
                otherwise
                    obj.grid = (0:M)'/M;
                    obj.kfR = [1e-6 0; 0 1e-2]; % placeholder
            end
        end

        function [st, obj] = update(obj, signal, xguess)
            % D4 grid-contract validation.
            if numel(signal) ~= numel(obj.grid)
                error('fdedge:gridLen', ...
                    'fdedge.tracker.Tracker: profile %s expects a signal of length %d; got %d', ...
                    obj.cfg.Profile, numel(obj.grid), numel(signal));
            end
            signal = signal(:);
            obj.n = obj.n + 1;
            n = obj.n;
            [xraw, braw] = rawFrame(obj, signal, xguess, n);
            % Kalman filter (from frame 2; PoC pass-2 structure).
            if n == 1
                xhat = xraw; bhat = braw; vhat = NaN; ohat = NaN;
                obj.xraw1 = xraw; obj.braw1 = braw;
            else
                if n == 2
                    if is2state(obj.cfg)
                        obj.kf = KalmanCV(2, obj.cfg.qV, obj.cfg.qW, obj.kfR);
                        obj.kf = obj.kf.init([obj.xraw1; xraw - obj.xraw1], ...
                            [obj.kfR(1) 0; 0 obj.kfR(1)]);
                    else
                        obj.kf = KalmanCV(4, obj.cfg.qV, obj.cfg.qW, obj.kfR);
                        Rx = obj.kfR(1,1); Rb = obj.kfR(2,2);
                        obj.kf = obj.kf.init([obj.xraw1; xraw-obj.xraw1; obj.braw1; braw-obj.braw1], ...
                            [Rx 0 0 0; 0 Rx*1e2 0 0; 0 0 Rb 0; 0 0 0 Rb*1e4]);
                    end
                end
                obj.kf = obj.kf.predict();
                if is2state(obj.cfg)
                    obj.kf = obj.kf.update([xraw]);
                    xhat = obj.kf.state(1); bhat = NaN;
                    vhat = obj.kf.state(2); ohat = NaN;
                else
                    obj.kf = obj.kf.update([xraw; braw]);
                    xhat = obj.kf.state(1); bhat = obj.kf.state(3);
                    vhat = obj.kf.state(2); ohat = obj.kf.state(4);
                end
            end
            obj.xrawPrev = xraw; obj.brawPrev = braw;
            st = struct('xhat', xhat, 'bhat', bhat, 'vhat', vhat, 'ohat', ohat, ...
                        'xraw', xraw, 'braw', braw);
        end

        function [xraw, braw] = rawFrame(obj, signal, xguess, n)
            % Per-profile per-frame block (line-for-line PoC move; local
            % function calls replaced by fdedge.* - same operation order).
            cfg = obj.cfg; M = cfg.M; x = obj.grid;
            alphas = cfg.alphas;
            switch lower(cfg.Profile)
                case 'poc'
                    % edge_fd_poc.m Part B loop body (per-frame part).
                    fs = smoothGauss(signal, cfg.sGrid);
                    Rn = caputoL1(x, fs, cfg.beta);
                    Pn = quantile(Rn, 0.90);
                    xraw = halfCross(Rn, x, 0.5*Pn);
                    xraw = xraw - obj.calib;   % scalar delta (shift invariance)
                    braw = NaN;
                case 'beta'
                    % edge_fd_beta_track.m per-frame block.
                    fs = smoothGauss(signal, cfg.sGrid);
                    Rbank = zeros(M+1, numel(alphas));
                    for ia = 1:numel(alphas)
                        Rbank(:,ia) = caputoGL(x, fs, alphas(ia));
                    end
                    ka = max(2, round((xguess + cfg.w0)*M) + 1);
                    kb = min(M+1, round((xguess + cfg.w1)*M) + 1);
                    if ~isfinite(xguess) || ka >= kb
                        % degenerate guess: carry over the previous frame's
                        % estimate (PoC; the driver updates xguess from the
                        % returned raw)
                        if n > 1
                            braw = obj.brawPrev; xraw = obj.xrawPrev;
                        else
                            braw = alphas(1); xraw = xguess;
                        end
                        return
                    end
                    dvec = x(ka:kb) - xguess;
                    c1 = zeros(size(alphas));
                    for ia = 1:numel(alphas)
                        y = Rbank(ka:kb, ia);
                        pf = polyfit(log(dvec), log(max(y, 1e-6)), 1);
                        c1(ia) = pf(1);   % slope only (polyfit returns [slope, intercept])
                    end
                    [braw, ~] = exponentBank(c1, alphas);
                    % refined position: slope peak at the ESTIMATED order,
                    % search restricted to +/-0.02 around the guess
                    Rstar = caputoGL(x, fs, braw);
                    xraw = slopePeak(Rstar, x, cfg.ell, xguess - 0.02, xguess + 0.02);
                    beta_ref = [0.55 0.65 0.75 0.85]';  % gridSlopePeak's axis (PoC)
                    dcal = interp2(beta_ref, alphas, obj.calib, braw, braw); % verbatim PoC call
                    xraw = xraw - dcal;
                case 'riesz'
                    % edge_fd_riesz_track.m per-frame block.
                    fs = smoothGauss(signal, cfg.sGrid);
                    Rbank = rieszLaplacian(fs, alphas);
                    if ~isfinite(xguess)
                        % degenerate guess: carry over the previous frame's
                        % estimate (PoC)
                        if n > 1
                            braw = obj.brawPrev; xraw = obj.xrawPrev;
                        else
                            braw = alphas(1); xraw = xguess;
                        end
                        return
                    end
                    xguess = min(max(xguess, cfg.xLo), cfg.xHi);   % keep wide probes in-domain
                    i0  = round(xguess*M);
                    ip  = round(cfg.dpr *M);
                    ipw = round(cfg.dprw*M);
                    % right-side probes, modular indices (periodic response
                    % depends only on the circular distance)
                    iR  = mod(i0 + ip  - 1, M) + 1;
                    iRw = mod(i0 + ipw - 1, M) + 1;
                    ghat = zeros(size(alphas));
                    for ia = 1:numel(alphas)
                        [gh, ok] = rieszGfit(Rbank(iR, ia), cfg.dpr);
                        if ~ok
                            [gh, ok] = rieszGfit(Rbank(iRw, ia), cfg.dprw);
                        end
                        if ~ok
                            gh = 0;
                        end
                        ghat(ia) = gh;
                    end
                    [braw, iaStar] = exponentBank(ghat, alphas);
                    xraw = cuspMin(Rbank(:,iaStar), x, xguess - 0.02, xguess + 0.02);
                    % debias (bilinear in xref and beta; alpha is the discrete
                    % bank index) - calibration grid axes as in the PoC
                    xes = (0.20:0.10:0.60)';
                    betas = alphas;
                    xcl = min(max(xguess, xes(1)), xes(end));
                    bcl = min(max(braw, betas(1)), betas(end));
                    dbr = interp2d(obj.calib.dbeta, xes, betas, xcl, bcl);
                    dxb = interp2d(obj.calib.dxg(:,iaStar,:), xes, betas, xcl, bcl);
                    xraw = xraw - dxb;
                    braw = braw - dbr;
                case 'hybrid'
                    % edge_fd_hybrid_track.m two-pass per-frame block.
                    fs = smoothGauss(signal, cfg.sGrid);
                    Rbank = rieszLaplacian(fs(1:M), alphas);
                    xguess = min(max(xguess, cfg.xLo), cfg.xHi);
                    % (a) exponent pass 1 at the prediction xguess -> order ia1
                    [b1, ia1] = estExponent(fs, alphas, xguess, x, cfg.w0, cfg.w1);
                    % (b) position: Riesz cusp minimum at ia1
                    xr = (0:M-1)'/M;
                    xraw = cuspMin(Rbank(:,ia1), xr, xguess - 0.02, xguess + 0.02);
                    % (c) exponent pass 2 at the fresh position xraw
                    [braw, iaStar] = estExponent(fs, alphas, xraw, x, cfg.w0, cfg.w1);
                    % (d) debias (bilinear in xref and beta) - the hybrid PoC
                    % clips xraw (not xguess); axes as in the PoC
                    xes = (0.20:0.10:0.60)';
                    betas = alphas;
                    xcl = min(max(xraw, xes(1)), xes(end));
                    bcl = min(max(braw, betas(1)), betas(end));
                    dbr = interp2d(obj.calib.dbeta, xes, betas, xcl, bcl);
                    dxb = interp2d(obj.calib.dxg(:,iaStar,:), xes, betas, xcl, bcl);
                    xraw = xraw - dxb;
                    braw = braw - dbr;
                otherwise
                    error('fdedge:tracker', ...
                        "unknown profile '%s' (expected 'poc'|'beta'|'riesz'|'hybrid')", cfg.Profile);
            end
        end
    end
end

function tf = is2state(cfg)
tf = strcmpi(cfg.Profile, 'poc');
end
