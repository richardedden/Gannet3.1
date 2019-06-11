function out = SignalFilter(in, r_flag, q_flag, MRS_struct)
% Based on: Cobas, J.C., Bernstein, M.A., Martín-Pastor, M., Tahoces, P.G.,
% 2006. A new general-purpose fully automatic baseline-correction procedure
% for 1D and 2D NMR data. J. Magn. Reson. 183, 145-51.

showPlots = 0;

[z.real, z.imag] = BaselineModeling(in, r_flag, q_flag, MRS_struct);

if showPlots == 1
    
    ii = MRS_struct.ii;
    freqRange = MRS_struct.p.sw(ii)/MRS_struct.p.LarmorFreq(ii);
    freq = (MRS_struct.p.npoints(ii) + 1 - (1:MRS_struct.p.npoints(ii))) / MRS_struct.p.npoints(ii) * freqRange + 4.68 - freqRange/2;
    
    H = figure(77);
    d.w = 0.5;
    d.h = 0.65;
    d.l = (1-d.w)/2;
    d.b = (1-d.h)/2;
    set(H,'Color', 'w', 'Units', 'Normalized', 'OuterPosition', [d.l d.b d.w d.h]);
    
    subplot(2,2,1); cla;
    plot(freq, real(in), 'k');
    set(gca, 'XDir', 'reverse', 'XLim', [0 6], 'Box', 'off');
    hold on;
    plot(freq, z.real, 'r');
    set(gca, 'XDir', 'reverse', 'XLim', [0 6], 'Box', 'off');
    
    subplot(2,2,2); cla;
    plot(freq, imag(in), 'k');
    set(gca, 'XDir', 'reverse', 'XLim', [0 6], 'Box', 'off');
    hold on;
    plot(freq, z.imag, 'r');
    set(gca, 'XDir', 'reverse', 'XLim', [0 6], 'Box', 'off');
    
    subplot(2,2,3); cla;
    plot(freq, real(in) - z.real, 'k');
    set(gca, 'XDir', 'reverse', 'XLim', [0 6], 'Box', 'off');
    
    subplot(2,2,4); cla;
    plot(freq, imag(in) - z.imag, 'k');
    set(gca, 'XDir', 'reverse', 'XLim', [0 6], 'Box', 'off');
    
    drawnow;
    %pause(0.5);
    
end

out = ifft(fftshift(complex(real(in) - z.real, imag(in) - z.imag)));

end


function [z_real, z_imag] = BaselineModeling(y, r_flag, q_flag, MRS_struct)
% Based on:
% Golotvin & Williams, 2000. Improved baseline recognition
%   and modeling of FT NMR spectra. J. Magn. Reson. 146, 122-125
% Cobas et al., 2006. A new general-purpose fully automatic
%   baseline-correction procedure for 1D and 2D NMR data. J. Magn. Reson.
%   183, 145-151

% Power spectrum of first-derivative of signal calculated by CWT
Wy = abs(cwt(real(y), 10, 'haar')).^2;

ii = MRS_struct.ii;
freqRange = MRS_struct.p.sw(ii)/MRS_struct.p.LarmorFreq(ii);
freq = (length(Wy) + 1 - (1:length(Wy))) / length(Wy) * freqRange + 4.68 - freqRange/2;
noiseLim = freq <= 12 & freq >= 10;

sigma = std(Wy(noiseLim));

w = 1:5;
k = 3;
baseline = zeros(length(Wy),1);
signal = zeros(length(Wy),1);

while 1
    if w(end) > length(Wy)
        break
    end
    h = max(Wy(w)) - min(Wy(w));
    if h < k*sigma
        baseline(w) = Wy(w);
    else
        signal(w) = Wy(w);
    end
    w = w + 1;
end

% Include lipids in baseline estimate and water, as appropriate
if r_flag
    lipidLim = freq <= 1.85 & freq >= -2;
    baseline(lipidLim) = Wy(lipidLim);
end
if q_flag
    waterLim = freq <= 5.5 & freq >= 4.25;
    baseline(waterLim) = Wy(waterLim);
end

z_real = real(y);
z_real(baseline == 0) = 0;
if r_flag
    lipids = whittaker(z_real(lipidLim), 2, 10);
end
if q_flag
    water = whittaker(z_real(waterLim), 2, 0.2);
end
z_real = whittaker(z_real, 2, 1e3);
if r_flag
    z_real(lipidLim) = lipids;
end
if q_flag
    z_real(waterLim) = water;
end

z_imag = -imag(hilbert(z_real));

end


function z = whittaker(y, d, lambda)
% Code taken from: Eilers, 2003. A perfect smoother. Anal. Chem. 75,
% 3631-3636

y = y(:);

m = length(y);
E = speye(m);
D = diff(E,d);

w = double(y ~= 0);
W = spdiags(w,0,m,m);
C = chol(W + lambda*(D'*D));
z = C\(C'\(w.*y));

end



