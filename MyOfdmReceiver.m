function [decoded_data, results] = MyOfdmReceiver(data, varargin)
% Optional: config
% config.packet_det_method: 'LTS' or 'STS'
% config.cfo_method: 'autocorr' or 'crosscorr'
% config.cfo_coarse: true/false (enable coarse CFO using STS)
% config.num_lts: 1 or 2 (number of LTS sequences for channel est)
% config.sfo_correction: true/false
% config.phase_correction: true/false

%% Setup
% OFDM_TX;

% Default configuration
config.packet_det_method = 'STS';
config.mod_order = 2;
config.n_ofdm_syms = 500;
config.cfo_method = 'crosscorr';
config.cfo_coarse = true;
config.num_lts = 1;
config.sfo_correction = true;
config.phase_correction = true;

% Parse optional arguments
if nargin > 1
    user_config = varargin{1};
    fields = fieldnames(user_config);
    for i = 1:length(fields)
        config.(fields{i}) = user_config.(fields{i});
    end
end

results.config = config;

%% Params:

% Waveform params
N_OFDM_SYMS             = config.n_ofdm_syms;         % Number of OFDM symbols
MOD_ORDER               =  config.mod_order;          % Modulation order in power of 2 (1/2/4/6 = BSPK/QPSK/16-QAM/64-QAM)

% OFDM params
SC_IND_PILOTS           = [8 22 44 58];                           % Pilot subcarrier indices
SC_IND_DATA             = [2:7 9:21 23:27 39:43 45:57 59:64];     % Data subcarrier indices
N_SC                    = 64;                                     % Number of subcarriers
CP_LEN                  = 16;                                     % Cyclic prefix length
trel = poly2trellis(7, [171 133]);              % Define trellis
scale = sqrt(2/(3*(2^MOD_ORDER-1)));

%% Preamble
% Preamble is a concatenation of multiple copies of STS and LTS
% It is used for packet detection and CFO and channel estimation
% LTS is sufficient to be used for the above three blocks in a way similar to what is given in OFDM thesis.
% If you want to use STS in place of LTS, read the paper below:
% 'Robust Frequency and Timing Synchronization for OFDM' by Timothy M. Schmidl and Donald C. Cox

% STS
sts_f = zeros(1,64);
sts_f(1:27) = [0 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0];
sts_f(39:64) = [0 0 1+1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0];
sts_t = ifft(sqrt(13/6).*sts_f, 64);
sts_t = sts_t(1:16);

% LTS
lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
lts_t = ifft(lts_f, 64);

preamble = [repmat(sts_t, 1, 30)  lts_t(33:64) lts_t lts_t];
rx_data = data;

%% 1. PACKET DETECTION
% LTS-based packet detection
lts_corr_out = zeros(1, length(rx_data));
for m = length(lts_t)+1:length(rx_data) - length(lts_t)
    lts_corr_out(m- length(lts_t)) = abs(sum(conj(rx_data(m - length(lts_t):m - 1)) .* rx_data(m:m+length(lts_t)-1)));
end
lts_corr_out = lts_corr_out / max(abs(lts_corr_out));

% STS-based packet detection (using multiple copies)
sts_corr_out = zeros(1, length(rx_data));
for m = length(sts_t)+1:length(rx_data) - length(sts_t)
    sts_corr_out(m- length(sts_t)) = abs(sum(conj(rx_data(m - length(sts_t):m - 1)) .* rx_data(m:m+length(sts_t)-1)));
end
sts_corr_out = sts_corr_out / max(abs(sts_corr_out));

% Store correlation results
results.lts_corr = lts_corr_out;
results.sts_corr = sts_corr_out;

% Select packet detection method
if strcmp(config.packet_det_method, 'LTS')
    corr_out = lts_corr_out;
    preamble_len = length(preamble);
else % STS
    corr_out = sts_corr_out;
    preamble_len = length(preamble);
end

% Find packet start (first significant peak)
LTS_CORR_THRESH = 0.7;
peaks = find(corr_out > LTS_CORR_THRESH);
plot(corr_out);
if isempty(peaks)
    error('No packet detected');
end

% Packet starts after preamble
pkt_start = peaks(1) + preamble_len;
results.packet_start = pkt_start;

% Extract preamble (STS and LTS)
preamble_data = rx_data(peaks(1):pkt_start-1);

% Extract payload
payload_len = (N_SC + CP_LEN) * N_OFDM_SYMS; % (N_SC + CP_LEN) * N_OFDM_SYMS
% if pkt_start + payload_len - 1 > length(rx_data)
%     payload_len = length(rx_data) - pkt_start + 1;
% end
payload_data = rx_data(pkt_start:pkt_start + payload_len - 1);

%% 2. CFO ESTIMATION AND CORRECTION

% Coarse CFO estimation using STS (if enabled)
cfo_coarse = 0;
if config.cfo_coarse
    % Last k STS copies
    k = 15;
    sts_kth_start = peaks(1) + (30 - k) * length(sts_t);
    sts_30th_start = peaks(1) + 29 * length(sts_t);
    runningsum = 0;
    for i = sts_kth_start : sts_30th_start
        runningsum = runningsum + conj(rx_data(i)) * rx_data(i + length(sts_t));
    end
    cfo_coarse = angle(runningsum)/ (2 * pi * length(sts_t));
end

results.cfo_coarse = cfo_coarse;

% Extract LTS from preamble
lts_start_idx = peaks(1) + 30*length(sts_t) + 0.5*length(lts_t);
lts_symbol_1 = rx_data(lts_start_idx:lts_start_idx+length(lts_t)-1);
lts_symbol_2 = rx_data(lts_start_idx+length(lts_t):lts_start_idx+2*length(lts_t)-1);

% Apply Coarse CFO correction to LTS
lts_1_cfo = lts_symbol_1 .* exp(-1j * 2 * pi * cfo_coarse * (0:length(lts_symbol_1)-1));
lts_2_cfo = lts_symbol_2 .* exp(-1j * 2 * pi * cfo_coarse * (0:length(lts_symbol_2)-1));


% Fine CFO via autocorrelation of LTS
rho_lts_autocorr = conj(lts_1_cfo) .* lts_2_cfo;
cfo_fine_lts_autocorr = angle(sum(rho_lts_autocorr))/(2 * pi * length(lts_t));

% Fine CFO via cross-correlation with known LTS
peak_1 = sum(lts_1_cfo .* conj(lts_t));
peak_2 = sum(lts_2_cfo .* conj(lts_t));
cfo_fine_lts_crosscorr = angle(peak_2 * conj(peak_1)) / (2 * pi * length(lts_t));

if strcmp(config.cfo_method, 'autocorr')
    results.cfo_fine = cfo_fine_lts_autocorr;
else
    results.cfo_fine = cfo_fine_lts_crosscorr;
end

cfo_total = results.cfo_coarse + results.cfo_fine;
results.cfo_total = cfo_total;

% Apply Fine CFO correction to LTS
lts_1_cfo = lts_1_cfo .* exp(-1j * 2 * pi * results.cfo_fine * (0:length(lts_symbol_1)-1));
lts_2_cfo = lts_2_cfo .* exp(-1j * 2 * pi * results.cfo_fine * (0:length(lts_symbol_2)-1));

% Apply Total CFO correction to payload
cfo_correction = exp(-1j * 2 * pi * cfo_total * (0:length(payload_data)-1));
payload_corrected = payload_data .* cfo_correction;

%% 3. CP REMOVAL
% Convert payload to matrix form and remove cyclic prefix
payload_mat = reshape(payload_corrected, N_SC + CP_LEN, N_OFDM_SYMS);

% Remove CP
payload_mat_nocp = payload_mat(CP_LEN+1:end, :);
%% 4. FFT
% Convert from time domain to frequency domain
payload_fft = fft(payload_mat_nocp, N_SC, 1);

%% 5. CHANNEL ESTIMATION
% FFT of LTS
lts_1_fft = fft(lts_1_cfo, N_SC);
lts_2_fft = fft(lts_2_cfo, N_SC);

% Channel estimate from LTS
if config.num_lts == 1
    H_est = lts_1_fft ./ lts_f;
else % Average over 2 LTS copies
    H_est = (lts_1_fft + lts_2_fft) / 2 ./ lts_f;
end

results.H_est = H_est;

% Channel equalization
H_est_mat = repmat(H_est(:), 1, N_OFDM_SYMS);
payload_equalized = payload_fft ./ H_est_mat;

%% 6. SFO CORRECTION (if enabled)
if config.sfo_correction
    % Extract pilots and estimate SFO
    pilots_rx = payload_equalized(SC_IND_PILOTS, :);
    
    % Expected pilots
    pilots_expected = repmat([1 1 -1 1]', 1, N_OFDM_SYMS);
    
    % Phase error over symbols
    phase_error = zeros(1, N_OFDM_SYMS);
    for sym_idx = 1:N_OFDM_SYMS
        phase_error(sym_idx) = angle(pilots_rx(:, sym_idx)' * conj(pilots_expected(:, sym_idx)));
    end
    
    % Estimate SFO as linear phase slope
    if N_OFDM_SYMS > 1
        sfo_coeff = polyfit(1:N_OFDM_SYMS, unwrap(phase_error), 1);
        sfo_slope = sfo_coeff(1);
    else
        sfo_slope = 0;
    end
    
    % Apply SFO correction
    sfo_correction = exp(-1j * sfo_slope * (0:N_OFDM_SYMS-1));
    payload_equalized = payload_equalized .* repmat(sfo_correction, N_SC, 1);
    
    results.sfo_slope = sfo_slope;
end

%% 7. PHASE ERROR CORRECTION (if enabled)
if config.phase_correction
    pilots_rx = payload_equalized(SC_IND_PILOTS, :);
    pilots_expected = repmat([1 1 -1 1]', 1, N_OFDM_SYMS);
    
    % Per-symbol phase correction
    for sym_idx = 1:N_OFDM_SYMS
        phase_error = angle(sum(pilots_rx(:, sym_idx) .* conj(pilots_expected(:, sym_idx))));
        payload_equalized(:, sym_idx) = payload_equalized(:, sym_idx) * exp(-1j * phase_error);
    end
end

%% 8. EXTRACT DATA AND DEMODULATION
% Remove pilots and flatten to vector
rx_syms = payload_equalized(SC_IND_DATA, :);
rx_syms = reshape(rx_syms, 1, numel(rx_syms));

% Store received symbols for constellation plot
results.rx_syms = rx_syms;

switch MOD_ORDER
    case 1
        constellation=scale*[1 -1];
    case 2
        constellation=scale*[1+1i -1+1i 1-1i -1-1i];
    case 4
        constellation=scale*[1+1i 1+3i 1-1i 1-3i 3+1i 3+3i 3-1i 3-3i -1+1i -1+3i -1-1i -1-3i -3+1i -3+3i -3-1i -3-3i];
    case 6
        constellation=scale*[3+3i 3+1i 3+5i 3+7i 3-3i 3-1i 3-5i 3-7i 1+3i 1+1i 1+5i 1+7i 1-3i 1-1i 1-5i 1-7i 5+3i 5+1i 5+5i 5+7i 5-3i 5-1i 5-5i 5-7i 7+3i 7+1i 7+5i 7+7i 7-3i 7-1i 7-5i 7-7i -3+3i -3+1i -3+5i -3+7i -3-3i -3-1i -3-5i -3-7i -1+3i -1+1i -1+5i -1+7i -1-3i -1-1i -1-5i -1-7i -5+3i -5+1i -5+5i -5+7i -5-3i -5-1i -5-5i -5-7i -7+3i -7+1i -7+5i -7+7i -7-3i -7-1i -7-5i -7-7i];
    otherwise
        error('wrong choice');
end

constDiag = comm.ConstellationDiagram( ...
    ReferenceConstellation=constellation, ...
    AxesLimits=[-2 2], ...
    Title='Rx Constellation');


constDiag(rx_syms(:));

% Demodulation using demapper
demap_out = demapper(rx_syms, MOD_ORDER, 1);

%% 9. FEC DECODING
% Viterbi decoder
decoded_data = vitdec(demap_out, trel, 7, 'trunc', 'hard');

% Remove padding zeros
% decoded_data = decoded_data(1:number_of_bits);

results.decoded_data = decoded_data;
results.demap_out = demap_out;

end
