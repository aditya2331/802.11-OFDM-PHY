function [tx_data, tx_payload_vec, tx_syms]=generate_ofdm_tx(tx_data,pilots, MOD_ORDER, number_of_bits, N_SC, CP_LEN, SC_IND_DATA, SC_IND_PILOTS, N_OFDM_SYMS, trellis_end_length)

% Forward Error Correction
tx_data = double([tx_data zeros(1,trellis_end_length) ]);    % 8 bits padding
trel = poly2trellis(7, [171 133]);              % Define trellis
tx_code = convenc(tx_data,trel);            % convultional encoder

% bits to signal space mapping these are you are x_k from the class
tx_syms = mapping(tx_code', MOD_ORDER, 1);

% figure(1);
% scatter(real(tx_syms), imag(tx_syms),'filled');
% title(' Signal Space of transmitted bits');
% xlabel('I'); ylabel('Q');

% Reshape the symbol vector to a matrix with one column per OFDM symbol,
tx_syms_mat = reshape(tx_syms, length(SC_IND_DATA), N_OFDM_SYMS);


% Repeat the pilots across all OFDM symbols
pilots_mat = repmat(pilots, 1, N_OFDM_SYMS);


% Construct the IFFT input matrix
ifft_in_mat = zeros(N_SC, N_OFDM_SYMS);

% Insert the data and pilot values; other subcarriers will remain at 0
ifft_in_mat(SC_IND_DATA, :)   = tx_syms_mat;
ifft_in_mat(SC_IND_PILOTS, :) = pilots_mat;

%Perform the IFFT --> frequency to time translation
tx_payload_mat = ifft(ifft_in_mat, N_SC, 1);

% Insert the cyclic prefix
if(CP_LEN > 0)
    tx_cp = tx_payload_mat((end-CP_LEN+1 : end), :);
    tx_payload_mat = [tx_cp; tx_payload_mat];
end

% Reshape to a vector
tx_payload_vec = reshape(tx_payload_mat, 1, numel(tx_payload_mat));


end
