function [snr,noise] = mySNRmap(img, annulussize)
% SNRMAP Generates a map of signal-to-noise ratio for a square image
%   Assumes signals follow multivariate gaussian distrubution
%   sectionHW is the half width of sections used to traverse the image, and
%   must be set to the best estimate for standard deviation of the signal
%   PSF


[rows, cols] = size(img);
dataArray = zeros(rows,cols,3);
snr = zeros(size(img));
img(img<0)=0;

%If not computing in parallel, change 'parfor' to 'for'
for r = 11:length(img(1,:))-11
    for c = 11:length(img(1,:))-11
        signal = img(r,c);

        [rr,cc] = meshgrid(1:size(img));

        %creates the logical mask for the annulus
        rad_signal = sqrt((ceil(rows/2)-r)^2+(ceil(cols/2)-c)^2);
        if rad_signal>length(img(1,:))/2-11
            continue
        end
        mask = sqrt((rr-ceil(rows/2)).^2+(cc-ceil(cols/2)).^2)>= (rad_signal-annulussize/2);
        mask2 = sqrt((rr-ceil(rows/2)).^2+(cc-ceil(cols/2)).^2)<= (rad_signal+annulussize/2);
        mask_annulus = mask&mask2;
        
        l = sqrt((rr-r).^2+(cc-c).^2);
        d = sqrt((rr-ceil(rows/2)).^2+(cc-ceil(cols/2)).^2);
        angles = acos((rad_signal^2-l.^2+d.^2)./(2*rad_signal*d));
        mask_theta = angles>=22*pi/180;
        mask_bpcsp = not(rr > 34 & rr < 46 & cc > 12 & cc <23);
        mask_incsp = not(rr > 38 & rr < 49 & cc > 52 & cc <63);
        mask_bpklip = not(rr > 21 & rr < 28 & cc > 45 & cc <54);
        
        mask_fin = mask_theta' & mask_annulus;% & mask_bpcsp';% & mask_incsp';% mask_bpklip';;
        %mask_fin = mask_annulus';
        %calculates the noise value
        mu_noise = mean(img(mask_fin));
        [s,n] = sumsqr(img(mask_fin)-mu_noise);
        noise(r,c) = sqrt(s/(n-1));
        noi=noise(r,c);

        if isnan(signal/noi) || isinf(signal/noi) ||(signal/noi)>100
            snr(r,c) = 1;
        else
            snr(r,c) = (signal-mu_noise)/(noi);
        end
    end
end
