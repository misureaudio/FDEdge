classdef test_gate < matlab.unittest.TestCase
    % Pins the two-tier replay gate: the numeric parser's tolerance rule,
    % the byte-identity fast path, and the CRLF triage path (D1/D2).
    methods (Test)
        function parserDigits(t)
            n = gate_parse('a = 0.2008  b = 1.932e-03  c = 1024  d = -0.1  e = +6.028e-04');
            t.verifyEqual(numel(n), 5);
            t.verifyEqual(n(1).v, 0.2008);
            t.verifyLessThanOrEqual(n(1).tol, 2e-4 + 1e-18);
            t.verifyEqual(n(2).v, 1.932e-3);
            t.verifyLessThanOrEqual(n(2).tol, 2e-6 + 1e-20);
            t.verifyEqual(n(3).v, 1024);
            t.verifyEqual(n(3).tol, 0);
            t.verifyEqual(n(4).v, -0.1);
            t.verifyLessThanOrEqual(n(4).tol, 2e-1 + 1e-18);
            t.verifyEqual(n(5).v, 6.028e-4);
            t.verifyLessThanOrEqual(n(5).tol, 2e-7 + 1e-21);
        end

        function identicalFilesPass(t)
            g = golden_file('results_poc.txt');
            [p, d] = gate_files(g, g);
            t.verifyTrue(p, d);
            t.verifyTrue(contains(d, 'byte-identical'), d);
        end

        function crlfIsTriaged(t)
            g = golden_file('results_poc.txt');
            % Explicit fid: the filename form fread(filename, format) is
            % broken in this MATLAB wrapper (always "Invalid file
            % identifier"); element-wise loop below is orientation-proof.
            fid = fopen(g, 'r');
            b = fread(fid, '*uint8');
            fclose(fid);
            % expand every LF to CRLF (loop: safe for any count)
            crlf = zeros(numel(b) + sum(b == 10), 1, 'uint8');
            k = 0;
            for i = 1:numel(b)
                k = k + 1; crlf(k) = b(i);
                if b(i) == 10
                    k = k + 1; crlf(k) = 13;
                end
            end
            g2 = fullfile(fileparts(g), '_crlf_test.txt');
            fid = fopen(g2, 'w'); fwrite(fid, crlf, 'uint8'); fclose(fid);
            [p, d] = gate_files(g2, g);
            delete(g2);
            t.verifyTrue(p, d);
            t.verifyTrue(contains(d, 'triaged'), d);
        end

        function numericDiffFailsStage1(t)
            g = golden_file('results_poc.txt');
            % corrupt one printed value beyond its tolerance: 2.1419 -> 2.1519
            fid = fopen(g, 'r');
            raw = fread(fid, '*uint8');
            fclose(fid);
            txt = char(raw(:))';              % force char row: regexp requires it
            % count + literal replace (4-arg regexp replace is misparsed here)
            nmatch = numel(regexp(txt, '2\.1419', 'match'));
            t.verifyEqual(nmatch, 1, 'expected the poc golden to contain 2.1419');
            txt = strrep(txt, '2.1419', '2.1519');
            g2 = fullfile(fileparts(g), '_num_test.txt');
            fid = fopen(g2, 'w'); fwrite(fid, uint8(txt), 'uint8'); fclose(fid);
            [p, d] = gate_files(g2, g);
            delete(g2);
            t.verifyFalse(p, d);
            t.verifyTrue(contains(d, 'stage1 FAIL'), d);
        end
    end
end
