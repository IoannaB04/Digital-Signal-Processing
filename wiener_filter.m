function z = wiener_filter(rank, data, data_noise, n)   
    telos = n+rank-1;
    rxx = xcorr(data_noise);      % autocorrelation of the noised waveform
    rxd = xcorr(data_noise,data); % autocorrelation of the two waveforms
    rxx = rxx(n:telos);
    rxd = rxd(n:telos);
    Rxx = toeplitz(rxx);
    w = pinv(Rxx,0.0001)*rxd;
    z = filter(w,1,data_noise);
end
