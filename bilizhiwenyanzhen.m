%% =========================================================
%  GEN_ratioSamples_KRON_case33bw_h3_SNR6levels_bus8to15_1e4.m
%
%  Goal:
%    For each internal bus in 8..15, build ratio feature under 3rd harmonic:
%      ratio = Ieq_L / Ieq_R
%    Under SNR = [50 45 40 35 30 25] dB, each (snr, bus) has 10000 samples.
%
%  Key points:
%    - Keep-set computed EXACTLY like your generator:
%         pmubranches + extract_nodes_from_branches
%      plus force include boundary buses (busL, busR).
%    - Harmonic solve: build Yh via makeYh_ref_from_mpc + optional Yload
%      enforce Vh(slack)=SlackVh (default 0)
%    - Kron reduction: Yred = kron_reduce_keep_set(Yh, keep_idx, [])
%      and Ieq = Yred * Vkeep
%    - Noise: complex AWGN added to [iL,iR] such that
%         SNR(dB) = 10*log10( Psig / Pnoise )
%      where Psig = mean(|iL|^2, |iR|^2) per sample
%
%  Output:
%    outDir\ratioSamples_case33bw_h3_bus8to15_SNRxx.mat
%    variables:
%      ratioMat  : [Nsamp x 8 x 6] complex, dims = (sample, busIndex, snrIndex)
%      meta      : struct with settings, mappings, keep-set info
%
%  Dependencies (must be on path):
%    MATPOWER: define_constants, loadcase, runpf, mpoption
%    your funcs: makeYh_ref_from_mpc, kron_reduce_keep_set, pmubranches, extract_nodes_from_branches
%% =========================================================
clc; clear; close all;

%% ===================== EDIT HERE ======================
outDir = '';

internal_bus_list = 8:15;     % nodes 8-15
harm = 3;                     % only 3rd harmonic
snr_list = [50 45 40 35 30 25];
Nsamp = 10000;

% EXACTLY match your generator rules (for Nint=8 => internal_end=15)
pmuL = 6;
busL = 7;

internal_end = internal_bus_list(end); % 15
pmuR = internal_end + 2;               % 17
busR = internal_end + 1;               % 16

UseLoadImpedance = false;   % match run_hpf default
SlackVh = 0;                % match run_hpf default
alpha_reg = 0;              % match run_hpf default

% injection randomness (in Amps)
Iamp_min = 0.5;
Iamp_max = 1.5;
% =======================================================

%% ---- dependency checks ----
assert(exist('define_constants','file')==2, 'MATPOWER not on path (define_constants missing).');
assert(exist('runpf','file')==2,           'MATPOWER not on path (runpf missing).');
assert(exist('makeYh_ref_from_mpc','file')==2, 'Missing makeYh_ref_from_mpc.m');
assert(exist('kron_reduce_keep_set','file')==2, 'Missing kron_reduce_keep_set.m');
assert(exist('pmubranches','file')==2, 'Missing pmubranches.m (used in generator).');
assert(exist('extract_nodes_from_branches','file')==2, 'Missing extract_nodes_from_branches.m (used in generator).');

define_constants;

if ~exist(outDir,'dir'); mkdir(outDir); end

%% ---- load case ----
mpc0 = loadcase('case33bw');

%% ---- base PF (same as run_hpf) ----
mpopt = mpoption('verbose', 0, 'out.all', 0);
res = runpf(mpc0, mpopt);
assert(res.success, 'Base power flow did not converge.');

nb = size(res.bus,1);
bus_ids = res.bus(:, BUS_I);
bidx = @(id) find(bus_ids==id, 1, 'first');

% slack index (REF bus)
slack = find(res.bus(:, BUS_TYPE) == REF, 1, 'first');
assert(~isempty(slack), 'No REF bus found.');

% base V1 (for load impedance option)
Vm = res.bus(:, VM);
Va = res.bus(:, VA) * pi/180;
V1 = Vm .* exp(1j*Va);

% load shunt at harmonics (same as run_hpf)
Yload_mat = sparse(nb, nb);
if UseLoadImpedance
    Sload_pu = (res.bus(:, PD) + 1j*res.bus(:, QD)) / res.baseMVA; % p.u.
    Yload    = conj(Sload_pu) ./ max(abs(V1).^2, 1e-6);            % p.u.
    Yload(isnan(Yload)) = 0;
    Yload_mat = sparse(1:nb, 1:nb, Yload, nb, nb);
end

%% ---- sanity checks for nodes ----
assert(~isempty(bidx(pmuL)), 'pmuL=%g not found', pmuL);
assert(~isempty(bidx(pmuR)), 'pmuR=%g not found', pmuR);
assert(~isempty(bidx(busL)), 'busL=%g not found', busL);
assert(~isempty(bidx(busR)), 'busR=%g not found', busR);
for k = 1:numel(internal_bus_list)
    assert(~isempty(bidx(internal_bus_list(k))), 'internal bus %g not found', internal_bus_list(k));
end

pmu_nodes = [pmuL, pmuR];

fprintf('Case: case33bw | nb=%d | baseMVA=%.3g\n', nb, res.baseMVA);
fprintf('Internal buses=%s | Boundary=[%g,%g] | PMU=[%g,%g]\n', ...
    mat2str(internal_bus_list), busL, busR, pmuL, pmuR);
fprintf('Harmonic=%d | Nsamp=%d | SNR=%s dB\n\n', harm, Nsamp, mat2str(snr_list));

%% =========================================================
% KEEP-SET (EXACTLY LIKE YOUR GENERATOR)
%% =========================================================
E = mpc0.branch(:,1:2);
row = (1:size(E,1))';
E = [E, row];

pmu_branches = pmubranches(E, pmu_nodes)';  % transpose like generator

branch = mpc0.branch;
branch = [(1:size(branch,1))', branch];     % prefix row id like generator expects

nodes_dengxiao = extract_nodes_from_branches(pmu_branches, branch);

keep_bus_id = unique([nodes_dengxiao(:); busL; busR], 'stable');
keep_idx    = arrayfun(bidx, keep_bus_id);

idxL = find(keep_bus_id == busL, 1);
idxR = find(keep_bus_id == busR, 1);
assert(~isempty(idxL) && ~isempty(idxR), 'busL/busR not in keep-set (should not happen).');

fprintf('keepSize=%d | keep_bus_id=%s\n', numel(keep_bus_id), mat2str(keep_bus_id'));
fprintf('idxL=%d (busL=%g) | idxR=%d (busR=%g)\n\n', idxL, busL, idxR, busR);

%% ---- build Yh once (harm=3) ----
[Ynet_h, ~] = makeYh_ref_from_mpc(mpc0, harm);
Ynet_h = sparse(Ynet_h);
if alpha_reg > 0
    Ynet_h = Ynet_h + alpha_reg * speye(nb);
end
Yh = Ynet_h + Yload_mat;

% Kron reduce
[Yred, ~, ~, ~] = kron_reduce_keep_set(Yh, keep_idx, []);
Yred = full(Yred);

% for enforced slack solve (same as your build_Fphys)
keepMask = true(nb,1); keepMask(slack) = false;
Ykk = Yh(keepMask, keepMask);
YkS = Yh(keepMask, slack);

Vh0 = zeros(nb,1);
Vh0(slack) = SlackVh;

%% ---- allocate output: (sample, busIndex, snrIndex) ----
Nb = numel(internal_bus_list);
Ns = numel(snr_list);
ratioMat = complex(zeros(Nsamp, Nb, Ns));

%% =========================================================
% MAIN LOOP
%% =========================================================
rng(1); % reproducible

for s = 1:Ns
    snr_db = snr_list(s);
    fprintf('--- SNR = %g dB ---\n', snr_db);

    for b = 1:Nb
        busnum = internal_bus_list(b);
        busrow = bidx(busnum);

        % base current conversion at this bus (A -> p.u.)
        Vbase_kV = res.bus(busrow, BASE_KV);
        assert(isfinite(Vbase_kV) && Vbase_kV > 0, 'BASE_KV invalid at bus %g', busnum);
        Ibase_A = (res.baseMVA*1e6) / (sqrt(3) * (Vbase_kV*1e3)); % A

        % generate random injections in A: magnitude uniform + random phase
        ampA  = Iamp_min + (Iamp_max - Iamp_min) * rand(Nsamp,1);
        ph   = 2*pi*rand(Nsamp,1);
        IinjA = ampA .* exp(1j*ph);

        % convert to pu
        Iinj_pu = IinjA / Ibase_A;

        % compute noiseless (iL,iR) sample-by-sample
        iL_true = complex(zeros(Nsamp,1));
        iR_true = complex(zeros(Nsamp,1));

        for k = 1:Nsamp
            Ih = zeros(nb,1);
            Ih(busrow) = Iinj_pu(k);

            rhs = Ih(keepMask) - YkS * Vh0(slack);
            Vh = Vh0;
            Vh(keepMask) = Ykk \ rhs;

            Vkeep = Vh(keep_idx);
            Ieq   = Yred * Vkeep;

            iL_true(k) = Ieq(idxL);
            iR_true(k) = Ieq(idxR);
        end

        % add complex AWGN to iL,iR with target SNR (per-sample power scale)
        % Psig = mean( |iL|^2 and |iR|^2 ) for that sample pair
        Psig = 0.5*(abs(iL_true).^2 + abs(iR_true).^2);
        Pnoise = Psig ./ (10.^(snr_db/10));     % because SNR(dB)=10log10(Psig/Pnoise)
        sigma = sqrt(Pnoise/2);                 % complex noise: E|n|^2 = 2*sigma^2

        nL = sigma .* (randn(Nsamp,1) + 1j*randn(Nsamp,1));
        nR = sigma .* (randn(Nsamp,1) + 1j*randn(Nsamp,1));

        iL_noisy = iL_true + nL;
        iR_noisy = iR_true + nR;

        ratioMat(:, b, s) = iL_noisy ./ iR_noisy;

        fprintf('  bus %g done. (mean|ratio|=%.3g)\n', busnum, mean(abs(ratioMat(:,b,s))));
    end

    % save each SNR separately (also keep full ratioMat in workspace)
    meta = struct();
    meta.case = 'case33bw';
    meta.harm = harm;
    meta.snr_db = snr_db;
    meta.Nsamp = Nsamp;
    meta.internal_bus_list = internal_bus_list;

    meta.pmu_nodes = pmu_nodes;
    meta.busL = busL; meta.busR = busR;
    meta.keep_bus_id = keep_bus_id;
    meta.keep_idx = keep_idx;
    meta.idxL = idxL; meta.idxR = idxR;

    meta.UseLoadImpedance = UseLoadImpedance;
    meta.SlackVh = SlackVh;
    meta.alpha_reg = alpha_reg;

    meta.injection_A_range = [Iamp_min, Iamp_max];
    meta.noise_model = 'complex AWGN on (iL,iR) with per-sample Psig';

    ratioSNR = squeeze(ratioMat(:,:,s)); % [Nsamp x 8]

    outPath = fullfile(outDir, sprintf('ratioSamples_case33bw_h3_bus8to15_SNR%d.mat', snr_db));
    save(outPath, 'ratioSNR', 'meta', '-v7.3');

    fprintf('[SAVED] %s | ratioSNR size=%dx%d\n\n', outPath, size(ratioSNR,1), size(ratioSNR,2));
end

fprintf('ALL DONE. outDir=%s\n', outDir);

%% ===================== FUNCTIONS ======================

function Fphys = build_Fphys_match_runhpf_1A( ...
    mpc, res, internal_bus_id, keep_idx, idxL, idxR, ...
    harm_list, Yload_mat, slack, SlackVh, alpha_reg)

    define_constants;

    nb   = size(res.bus,1);
    Nh   = numel(harm_list);
    Nint = numel(internal_bus_id);

    % map BUS ID -> row in res.bus
    bus_ids = res.bus(:, BUS_I);

    % slack elimination mask
    keepMask = true(nb,1); keepMask(slack) = false;

    Fphys = complex(zeros(4*Nh, Nint));

    for hh = 1:Nh
        h = harm_list(hh);

        % --- build Yh (p.u.) exactly like run_hpf ---
        [Ynet_h, ~] = makeYh_ref_from_mpc(mpc, h);
        Ynet_h = sparse(Ynet_h);
        if alpha_reg > 0
            Ynet_h = Ynet_h + alpha_reg * speye(nb);
        end
        Yh = Ynet_h + Yload_mat;

        % Kron reduce to keep-set (indices)
        [Yred, ~, ~, ~] = kron_reduce_keep_set(Yh, keep_idx, []);
        Yred = full(Yred);

        % reduced blocks for enforced slack V
        Ykk = Yh(keepMask, keepMask);
        YkS = Yh(keepMask, slack);

        % slack fixed
        Vh0 = zeros(nb,1);
        Vh0(slack) = SlackVh;

        for j = 1:Nint
            busnum = internal_bus_id(j);

            busrow = find(bus_ids == busnum, 1, 'first');
            if isempty(busrow)
                error('Internal bus %d not found in res.bus BUS_I.', busnum);
            end

            Vbase_kV = res.bus(busrow, BASE_KV);
            if ~(isfinite(Vbase_kV) && Vbase_kV > 0)
                error('BASE_KV missing/invalid at bus %d.', busnum);
            end

            % 1 A -> p.u. (same as your run_hpf)
            Ibase_A = (res.baseMVA*1e6) / (sqrt(3) * (Vbase_kV*1e3)); % A
            Ipu_mag = 1.0 / Ibase_A;  % 1A in p.u. magnitude

            Ih = zeros(nb,1);
            Ih(busrow) = Ipu_mag;      % phase = 0 deg for basis

            rhs = Ih(keepMask) - YkS * Vh0(slack);
            Vh = Vh0;
            Vh(keepMask) = Ykk \ rhs;

            Vkeep = Vh(keep_idx);
            Ieq   = Yred * Vkeep;

            iL = Ieq(idxL);
            iR = Ieq(idxR);
            vL = Vkeep(idxL);
            vR = Vkeep(idxR);

            r0 = (hh-1)*4;
            Fphys(r0+1:r0+4, j) = [iL; iR; vL; vR];
        end
    end
end