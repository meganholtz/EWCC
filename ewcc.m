function [ FFTdata ] = ewcc( data , useWindow, center)
%ewcc calculates the complex EWCC transform fft(log(data)) for 2d or 4d data
%   input:
%       data -- 2d or 4d scanning electron diffraction data, ordered (kx,
%               ky, x, y)
%       useWindow -- logical indicating whether or not to apply hann window
%               before FFT.  Default is 1.  The window is useful to prevent
%               FFT artifacts from non-periodicity.  This is especially
%               important it the diffraction patterns have significant
%               intensity at their edges.
%       center -- the known center of the diffraction pattern to shift to -
%               otherwise mis-centering artifacts may appear in the
%               im(ewcc)

if nargin==1
    useWindow = 1;
    center = floor([size(data,1),size(data,2)]/2)+1;
elseif nargin==2
    center = floor([size(data,1),size(data,2)]/2)+1;
end
if numel(center)==1
    center = center*[1,1];
end

[N_kx,N_ky,N_x,N_y]=size(data);
minval=min(data(:));
offsetval = (max(data(:)) - minval)*1e-3;

if useWindow
    win=window2(N_kx,N_ky,@hann);
else
    win=ones(N_kx,N_ky);
end


fftComp = @(x) fftshift(fft2(fftshift( win.*log(x-minval+offsetval)))) ;

% Convert to FFT of CBED map
FFTdata = data;
if nargin<3 % no centering needed
    for x=1:N_x
        for y=1:N_y
            FFTdata(:,:,x,y) = fftComp(FFTdata(:,:,x,y));
        end
    end
else % need to center the diffraction pattern first
    for x=1:N_x
        for y=1:N_y
            FFTdata(:,:,x,y) = fftComp(centerPattern(FFTdata(:,:,x,y),center,useWindow));
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
