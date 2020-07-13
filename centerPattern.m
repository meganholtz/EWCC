function centered = centerPattern(pattern,center,useWindow)
% centerPattern shifts the CBED to be centered for 2d or 4d data
%   input:
%       pattern -- 2d or 4d scanning electron diffraction data, ordered (kx,
%               ky, x, y)
%       center -- the known center of the diffraction pattern to shift to 
%       useWindow -- logical indicating whether or not to apply hann window
%               before FFT.  Default is 1.  The window is useful to prevent
%               FFT artifacts from non-periodicity.  This is especially
%               important it the diffraction patterns have significant
%               intensity at their edges.


N = size(pattern);
f1 = (0:N(1)-1)/N(1);   % frequency grid in direction 1
f2 = (0:N(2)-1)/N(2);   % frequency grid in direction 1
sh = 1-center;          % amount to shift by
phf = exp(-2*pi*1i*sh(1)*f1')*exp(-2*pi*1i*sh(2)*f2); % phase factor

if nargin==2
    useWindow = 1;
end
if useWindow==1
    win=window2(N(1),N(2),@hann);
else
    win=ones(N(1),N(2));
end

centered = 0*pattern;
if length(N)==2
    centered =  abs( fftshift(ifft2( win.*fftshift(fft2(win.*pattern)).*phf )));
elseif length(N)==3
    for i=1:N(3)
        centered(:,:,i) = ...
            abs( fftshift(ifft2( win.*fftshift(fft2(win.*pattern(:,:,i))).*phf )));
    end
elseif length(N)==4
    for i=1:N(3)
        for j=1:N(4)
            centered(:,:,i,j) = ...
                abs( fftshift(ifft2( win.*fftshift(fft2(win.*pattern(:,:,i,j))).*phf )));
        end
    end
end
end

function w=window2(N,M,w_func)

wc=window(w_func,N);
wr=window(w_func,M);
[maskr,maskc]=meshgrid(wr,wc);

w=maskr.*maskc;

end