function center = findCenterEWCCphase(nbed, center, precision, range, mask, doPlot)
% findCenterEWCCphase minimizes the imaginary intensity of EWCC to find
% center of the diffraction pattern
%   input:
%       nbed -- 2d scanning electron diffraction data, ordered 
%       center -- the estimated center of the diffraction pattern
%       precision -- how many orders of magnitude to check
%       range -- how many pixels to check over
%       mask -- where to check the imaginary EWCC 
%       doPlot -- plot the results, true/false

N = size(nbed);

c0 = floor(N/2)+1;

if nargin<2
    center = c0;
    precision = 2;
    range = 5;
elseif nargin<3
    precision = 2;
    range = 5;
elseif nargin<4
    range = 5;
end
if nargin<5
    mask = ones(N);
    mask(c0(1)-3:c0(1)+3, c0(2)-3:c0(2)+3) = 0;
end
if nargin<6
    doPlot = 0;
end

if doPlot
    F=figure;
end

for p = 0:precision
    
    dc = -range:0.5:range;
    dc = dc/10^p; % smaller and smaller increments
    
    
    imew = zeros(N(1),N(2),length(dc),length(dc));
    
    for i = 1:length(dc)
        for j = 1:length(dc)
            
            imew(:,:,i,j) = imag(cewpc(nbed,1, center+[dc(i),dc(j)]));
            
        end
    end
    
    imewm = imew .* repmat(mask,[1,1,size(imew,3),size(imew,4)]);
    %browseSTEM4D(bsat(imew))
    
    minmat = squeeze(sum(sum(abs(imewm),1),2));
    [~,ind] = min(minmat(:));
    [ind1,ind2] = ind2sub(size(minmat),ind);
    
    center = center+[dc(ind1),dc(ind2)];
    
    if doPlot
        
        disp([' testing over [' num2str(min(dc)) ' : ' num2str(dc(2)-dc(1)) ' : ' num2str(max(dc)) ']'])
        
        disp(['center found: '])
        disp(center)
        
        figure(F)
        subplot(precision+1,2,2*p+1); plotIM(minmat); hold on;
        plot(ind2,ind1,'rx')
        title([num2str(center(1)),', ' num2str(center(2)) ' / ' ...
            num2str((ind1)),', ' num2str((ind2))])
        
        subplot(precision+1,2,2*(p+1));
        plotIM(bsat(imew(:,:,ind1,ind2))); colorbar
        
    end
end