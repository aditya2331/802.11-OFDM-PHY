function [tx_vec_air_A]=generate_mimo(preamble_A,tx_payload_vec_A, INTERP_RATE, TX_SCALE,CFO_FLAG, interp_filt2)


% Construct the full time-domain OFDM waveform
tx_vec_A = [preamble_A tx_payload_vec_A];

% Pad with zeros for transmission
tx_vec_padded_A = [tx_vec_A zeros(1,50)];

%% Interpolate

% Interpolate
if(INTERP_RATE == 1)
    tx_vec_air_A = tx_vec_padded_A;
elseif(INTERP_RATE == 2)
    % Zero pad then filter (same as interp or upfirdn without signal processing toolbox)
    tx_vec_2x_A = zeros(1, 2*numel(tx_vec_padded_A));
    tx_vec_2x_A(1:2:end) = tx_vec_padded_A;
    tx_vec_air_A = filter(interp_filt2, 1, tx_vec_2x_A);
    
end

% Scale the Tx vector to +/- 1
tx_vec_air_A = TX_SCALE .* tx_vec_air_A ./ max(abs(tx_vec_air_A));

TX_NUM_SAMPS = length(tx_vec_air_A);

% Perfect (ie. Rx=Tx):
% rx_vec_air = tx_vec_air;

% to enable CFO make CFO_FLAG=1

if(CFO_FLAG)
    tx_vec_air_A = tx_vec_air_A .* exp(-1i*2*pi*1e-4*[0:length(tx_vec_air_A)-1]);
end

end
