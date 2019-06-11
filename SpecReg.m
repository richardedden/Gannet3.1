function out = SpecReg(input, target, x)

f   = x(1);
phi = x(2);

t = 0:input.dwelltime:(length(input.data)/2-1)*input.dwelltime;
input.data = reshape(input.data, [length(input.data)/2 2]);
fid = complex(input.data(:,1), input.data(:,2));

fid = fid .* exp(1i*pi*(t'*f*2+phi/180));
fid = [real(fid); imag(fid)];

out = target - fid;

end