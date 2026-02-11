clc
if isappdata(0, 'snr_override')
    snr = getappdata(0, 'snr_override'); % in dB
else
    snr = 10; % in dB (default)
end
WRITE_PNG_FILES         = 0;           % Enable writing plots to PNG

%% Params:

INTERP_RATE=2;
SAMP_FREQ=20e6;
CHANNEL = 11;
simple_channel_flag = 1; % flag to use simple channel
CFO_FLAG = 0; % flag to enable CFO
DETECTION_OFFSET = 100; % to add packet detection error

TX_SPATIAL_STREAM_SHIFT=3;

% Waveform params
N_OFDM_SYMS             = 50;         % Number of OFDM symbols
MOD_ORDER               =  2;          % Modulation order in power of 2 (1/2/4/6 = BSPK/QPSK/16-QAM/64-QAM)
TX_SCALE                = 1.0;         % Scale for Tx waveform ([0:1])

% OFDM params
SC_IND_PILOTS           = [8 22 44 58];                           % Pilot subcarrier indices
SC_IND_DATA             = [2:7 9:21 23:27 39:43 45:57 59:64];     % Data subcarrier indices
N_SC                    = 64;                                     % Number of subcarriers
CP_LEN                  = 16;                                     % Cyclic prefix length
N_DATA_SYMS             = N_OFDM_SYMS * length(SC_IND_DATA);      % Number of data symbols (one per data-bearing subcarrier per OFDM symbol)

channel_coding = .5; % coding rate
trellis_end_length = 8; % bits for trellis to end


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

% We break the construction of our preamble into two pieces. First, the
% legacy portion, is used for CFO recovery and timing synchronization at
% the receiver. The processing of this portion of the preamble is SISO.
% Second, we include explicit MIMO channel training symbols.

% Legacy Preamble

% Use 30 copies of the 16-sample STS for extra AGC settling margin
% To avoid accidentally beamforming the preamble transmissions, we will
% let RFA be dominant and handle the STS and first set of LTS. We will
% append an extra LTS sequence from RFB so that we can build out the
% channel matrix at the receiver

sts_t_rep = repmat(sts_t, 1, 30);

preamble_legacy_A = [sts_t_rep, lts_t(33:64), lts_t, lts_t];
preamble_legacy_B = [circshift(sts_t_rep, [0, TX_SPATIAL_STREAM_SHIFT]), zeros(1, 160)];

% MIMO Preamble

% There are many strategies for training MIMO channels. Here, we will use
% the LTS sequence defined before and orthogonalize over time. First we
% will send the sequence on stream A and then we will send it on stream B

preamble_mimo_A = [lts_t(33:64), lts_t, zeros(1,96)];
preamble_mimo_B = [zeros(1,96), lts_t(33:64), lts_t];

preamble_A = [preamble_legacy_A, preamble_mimo_A];
preamble_B = [preamble_legacy_B, preamble_mimo_B];

% Define a half-band 2x interpolation filter response
interp_filt2 = zeros(1,43);
interp_filt2([1 3 5 7 9 11 13 15 17 19 21]) = [12 -32 72 -140 252 -422 682 -1086 1778 -3284 10364];
interp_filt2([23 25 27 29 31 33 35 37 39 41 43]) = interp_filt2(fliplr([1 3 5 7 9 11 13 15 17 19 21]));
interp_filt2(22) = 16384;
interp_filt2 = interp_filt2./max(abs(interp_filt2));

pilots_A= [1 1 -1 1].';
pilots_B= [0 0 0 0].';

%% Generate a payload of random integers
number_of_bits= (N_DATA_SYMS * MOD_ORDER - 2*trellis_end_length) * channel_coding;
tx_data = randi(2, 1, number_of_bits) - 1;

% we will use same tx data for both the TX antenna
[tx_data_a, tx_payload_vec_A,tx_syms_A]=generate_ofdm_tx(tx_data, pilots_A, MOD_ORDER, number_of_bits, N_SC, CP_LEN, SC_IND_DATA, SC_IND_PILOTS, N_OFDM_SYMS, trellis_end_length);
[tx_data_b, tx_payload_vec_B,tx_syms_B]=generate_ofdm_tx(tx_data, pilots_B, MOD_ORDER, number_of_bits, N_SC, CP_LEN, SC_IND_DATA, SC_IND_PILOTS, N_OFDM_SYMS, trellis_end_length);

[tx_vec_air_A]=generate_mimo(preamble_A,tx_payload_vec_A, INTERP_RATE, TX_SCALE,CFO_FLAG,interp_filt2);
[tx_vec_air_B]=generate_mimo(preamble_B,tx_payload_vec_B, INTERP_RATE, TX_SCALE,CFO_FLAG,interp_filt2);

tx_syms_space = [tx_syms_A; tx_syms_B];
tx_syms = reshape(tx_syms_space, 1, length(tx_syms_A)*2);


%% MIMO air channel case 2 (1X4 MIMO)
rx_vec_air_1A = tx_vec_air_A;
rx_vec_air_1B= .5 * tx_vec_air_A;
rx_vec_air_1C = tx_vec_air_A;
rx_vec_air_1D= .5 * tx_vec_air_A;


rx_vec_air_1A = [rx_vec_air_1A, zeros(1,100)];
rx_vec_air_1B = [rx_vec_air_1B, zeros(1,100)];
rx_vec_air_1C = [rx_vec_air_1C, zeros(1,100)];
rx_vec_air_1D = [rx_vec_air_1D, zeros(1,100)];

noise_power = var(rx_vec_air_1A) * 10.^(-snr/20);

rx_vec_air_1A = rx_vec_air_1A + noise_power.*complex(randn(1,length(rx_vec_air_1A)), randn(1,length(rx_vec_air_1A)));
rx_vec_air_1B = rx_vec_air_1B + noise_power.*complex(randn(1,length(rx_vec_air_1B)), randn(1,length(rx_vec_air_1B)));
rx_vec_air_1C = rx_vec_air_1C + noise_power.*complex(randn(1,length(rx_vec_air_1C)), randn(1,length(rx_vec_air_1C)));
rx_vec_air_1D = rx_vec_air_1D + noise_power.*complex(randn(1,length(rx_vec_air_1D)), randn(1,length(rx_vec_air_1D)));

rx_vec_dec_1A=generate_downsample_rx(rx_vec_air_1A,INTERP_RATE, interp_filt2);
rx_vec_dec_1B=generate_downsample_rx(rx_vec_air_1B,INTERP_RATE, interp_filt2);
rx_vec_dec_1C=generate_downsample_rx(rx_vec_air_1C,INTERP_RATE, interp_filt2);
rx_vec_dec_1D=generate_downsample_rx(rx_vec_air_1D,INTERP_RATE, interp_filt2);

%% MIMO air channel case 3 (2X2 MIMO, 1 spatial streams)
rx_vec_air_2A = tx_vec_air_A + .5 * tx_vec_air_B;
rx_vec_air_2B = tx_vec_air_B + .5 * tx_vec_air_A;

rx_vec_air_2A = [rx_vec_air_2A, zeros(1,100)];
rx_vec_air_2B = [rx_vec_air_2B, zeros(1,100)];

rx_vec_air_2A = rx_vec_air_2A + noise_power*complex(randn(1,length(rx_vec_air_2A)), randn(1,length(rx_vec_air_2A)));
rx_vec_air_2B = rx_vec_air_2B + noise_power*complex(randn(1,length(rx_vec_air_2B)), randn(1,length(rx_vec_air_2B)));

raw_rx_dec_2A=generate_downsample_rx(rx_vec_air_2A,INTERP_RATE, interp_filt2);
raw_rx_dec_2B=generate_downsample_rx(rx_vec_air_2B,INTERP_RATE, interp_filt2);

%% MIMO air channel case 4 (2X2 MIMO, 2 spatial streams)
tx_data_A_bits = randi(2, 1, number_of_bits) - 1;   % Stream A bits
tx_data_B_bits = randi(2, 1, number_of_bits) - 1;   % Stream B bits 



[tx_data_a, tx_payload_vec_A, tx_syms_A] = generate_ofdm_tx( ...
    tx_data_A_bits, pilots_A, MOD_ORDER, number_of_bits, N_SC, CP_LEN, ...
    SC_IND_DATA, SC_IND_PILOTS, N_OFDM_SYMS, trellis_end_length);

[tx_data_b, tx_payload_vec_B, tx_syms_B] = generate_ofdm_tx( ...
    tx_data_B_bits, pilots_B, MOD_ORDER, number_of_bits, N_SC, CP_LEN, ...
    SC_IND_DATA, SC_IND_PILOTS, N_OFDM_SYMS, trellis_end_length);

tx_vec_air_A = generate_mimo(preamble_A, tx_payload_vec_A, INTERP_RATE, TX_SCALE, CFO_FLAG, interp_filt2);
tx_vec_air_B = generate_mimo(preamble_B, tx_payload_vec_B, INTERP_RATE, TX_SCALE, CFO_FLAG, interp_filt2);

%% simple channel (scalar multiplication)
% [y1]   [h11 h12] [xA] + [n1]
% [y2] = [h21 h22] [xB]   [n2]
h11 = 0.3;
h12 = 0.5;
h21 = 0.7;
h22 = 1.0;

H = [0.3, 0.5; 0.7, 1.0];
[U,S,V] = svd(H);

precoded_tx_vec_air_A = V(1,1)*tx_vec_air_A + V(1,2)*tx_vec_air_B;
precoded_tx_vec_air_B = V(2,1)*tx_vec_air_A + V(2,2)*tx_vec_air_B;

rx_vec_air_2A = h11*precoded_tx_vec_air_A + h12*precoded_tx_vec_air_B;   % Rx antenna 1
rx_vec_air_2B = h21*precoded_tx_vec_air_A + h22*precoded_tx_vec_air_B;   % Rx antenna 2

noise_amp_A = sqrt(var(rx_vec_air_2A)) * 10^(-snr/20);
noise_amp_B = sqrt(var(rx_vec_air_2B)) * 10^(-snr/20);

rx_vec_air_2A = rx_vec_air_2A + noise_amp_A * complex(randn(1,length(rx_vec_air_2A)), randn(1,length(rx_vec_air_2A)));
rx_vec_air_2B = rx_vec_air_2B + noise_amp_B * complex(randn(1,length(rx_vec_air_2B)), randn(1,length(rx_vec_air_2B)));

raw_rx_dec_3A = generate_downsample_rx(rx_vec_air_2A, INTERP_RATE, interp_filt2);
raw_rx_dec_3B = generate_downsample_rx(rx_vec_air_2B, INTERP_RATE, interp_filt2);
%% Multipath channel
L = 4;                            %number of taps
pwr = exp(-(0:L-1));              %power delay profile
pwr = pwr / sum(pwr);             %normalize total power

% Independent complex Gaussian taps per link 
h11 = (randn(L,1)+1i*randn(L,1)).' .* sqrt(pwr/2);
h12 = (randn(L,1)+1i*randn(L,1)).' .* sqrt(pwr/2);
h21 = (randn(L,1)+1i*randn(L,1)).' .* sqrt(pwr/2);
h22 = (randn(L,1)+1i*randn(L,1)).' .* sqrt(pwr/2);

% Apply multipath (FIR filtering) 
rx_vec_air_2A = filter(h11, 1, tx_vec_air_A) + filter(h12, 1, tx_vec_air_B);
rx_vec_air_2B = filter(h21, 1, tx_vec_air_A) + filter(h22, 1, tx_vec_air_B);

energy = abs([rx_vec_air_2A; rx_vec_air_2B]).^2;
sigPow = mean(energy(:));         
noisePow = sigPow * 10^(-snr/10);       
n1 = sqrt(noisePow/2) * (randn(size(rx_vec_air_2A)) + 1i*randn(size(rx_vec_air_2A)));
n2 = sqrt(noisePow/2) * (randn(size(rx_vec_air_2B)) + 1i*randn(size(rx_vec_air_2B)));

rx_vec_air_2A = rx_vec_air_2A + n1;
rx_vec_air_2B = rx_vec_air_2B + n2;
if simple_channel_flag == 0
    % raw_rx_dec_3A = generate_downsample_rx(rx_vec_air_2A, INTERP_RATE, interp_filt2);
    % raw_rx_dec_3B = generate_downsample_rx(rx_vec_air_2B, INTERP_RATE, interp_filt2);
end
