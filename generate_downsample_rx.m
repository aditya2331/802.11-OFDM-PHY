function raw_rx_dec_A=generate_downsample_rx(rx_vec_air_2A, INTERP_RATE, interp_filt2)

%% Decimate
if(INTERP_RATE == 1)
    raw_rx_dec_A = rx_vec_air_2A;
elseif(INTERP_RATE == 2)
    raw_rx_dec_A = filter(interp_filt2, 1, rx_vec_air_2A);
    raw_rx_dec_A = raw_rx_dec_A(1:2:end);
end


end
