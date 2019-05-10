function [newImg,dataArray] = snrMap2all(img, sectionHW)
% SNRMAP Generates a map of signal-to-noise ratio for a square image
%   Assumes signals follow multivariate gaussian distrubution
%   sectionHW is the half width of sections used to traverse the image, and
%   must be set to the best estimate for standard deviation of the signal
%   PSF


[rows, cols] = size(img);
dataArray = zeros(rows,cols,3);
newImg = zeros(size(img));

%Creates the 2D multivariate Guassian distrubution model
%Parameters are sigma,x0,y0
f = fittype('a1*exp(-(x-sectionHW-1)^2/(2*sigma^2)-(y-sectionHW-1)^2/(2*sigma^2))','independent',{'x','y'},'dependent','z','problem','sectionHW');

%If not computing in parallel, change 'parfor' to 'for'
for r = 1:rows
    for c = 1:cols
        
        %fprintf('(%f,%f)\n',r,c)
        
        %computes the range of column and row values of current section
        c_vals = c-sectionHW:c+sectionHW;
        r_vals = r-sectionHW:r+sectionHW;
        
        %handles edge cases by truncation
        %corrected column and row ranges are stored as c_vals_corr and
        %r_vals_corr
        regc = c_vals >= 1 & c_vals<= cols;
        regr = r_vals >= 1 & r_vals <= rows;
        c_vals_corr = c_vals(regc);
        r_vals_corr = r_vals(regr);
        
        %Generates the Gaussian fit
        %Initial values of x0,y0,a1, and sigma are provided by the location
        %of the maximum value in the section, the intensity at that
        %location, and the section halfwidth respectively
        z = img(r_vals_corr,c_vals_corr);
        imagemean = mean(mean(img));
        
        if img(r,c)<imagemean
            newImg(r,c) = 1;
        else
            [x,y] = meshgrid(1:length(r_vals_corr),1:length(c_vals_corr));
            x = x(:); y = y(:); z_alt = z(:);
            [M,I] = max(z_alt);
            [I1,I2] = ind2sub(size(z),I);
            [sf,gof] = fit([x,y],z_alt,f,'StartPoint',[z(I1,I2),sectionHW],'problem',sectionHW);
            coefficients = coeffvalues(sf);
            sig_real = coefficients(2);
            
            %sets paramaters to remove outliers, these can be altered
            %Currently sets the SNR to zero if:
            % 1) the standard devation of the fit is greater than an upper bound
            % reasonable number
            if abs(sig_real)>20
                newImg(r,c) = 1;
            else
                %Signal is calculated by finding the value of the fit at the
                %current location
                signal = coefficients(1);
                
                %NOISE CALCS
                %Sets up the area for the noise to be evaluated on
                signal_c = r;
                signal_r = c;
                area_half_width = 2.355*sig_real;
                if area_half_width < 1
                    area_half_width = 1;
                end
                img_gauss = zeros(rows,cols);
                [rr,cc] = meshgrid(1:size(img_gauss));
                
                %creates the logical mask for the annulus
                rad_signal = sqrt( (rows/2-signal_r)^2+(cols/2-signal_c)^2);
                mask = sqrt( (rr-rows/2).^2+(cc-cols/2).^2)>= (rad_signal-area_half_width);
                mask2 = sqrt((rr-rows/2).^2+(cc-cols/2).^2)<= (rad_signal+area_half_width);
                mask_annulus = mask&mask2;
                
                %creates the logical mask for the wedge section
                l = sqrt((rr-signal_r).^2+(cc-signal_c).^2);
                d = sqrt((rr-rows/2).^2+(cc-cols/2).^2);
                angles = acos((rad_signal^2-l.^2+d.^2)./(2*rad_signal*d));
                %angles(imag(angles)~=0) = 0;
                maxAngle = rad_signal/rows;
                if maxAngle<.3
                    maxAngle = .3;
                end
                %disp(maxAngle)
                %mask_theta = angles<= 2*area_half_width/rad_signal;
                mask_theta = angles<= maxAngle;
                
                mask_fin = mask_theta & mask_annulus;
                %image(mask_fin,'CDataMapping','Scaled')
                %pause
                %mask_fin = mask_annulus;
                
                %generates the fit model of the signal
                for r_t = 1:rows
                    for c_t = 1:cols
                        img_gauss(r_t,c_t) = coefficients(1)*exp( -( (r_t-r)^2 / (2*sig_real^2) + (c_t-c)^2/(2*sig_real^2)));
                    end
                end
                
                %subtracts the fit model from the data
                img_subtracted = img-img_gauss;
                %image(img_subtracted,'CDataMapping','Scaled')
                %pause
                
                %calculates the noise value
                mu_noise = mean(img_subtracted(mask_fin));
                [s,n] = sumsqr(img_subtracted(mask_fin)-mu_noise);
                noise = sqrt(s/(n-1));
                
                if isnan(signal/noise) || isinf(signal/noise) ||(signal/noise)>100
                    newImg(r,c) = 1;
                else
                    newImg(r,c) = (signal-mu_noise)/(noise);
                    %newImg(r,c) = signal;
                end
                %centerArray(int8(coefficients(4)+r_vals_corr(1)),int8(coefficients(3)+c_vals_corr(1))) = centerArray(int8(coefficients(4)+r_vals_corr(1)),int8(coefficients(3)+c_vals_corr(1)))+1;
                dataArray(r,c,:) = [signal-mu_noise,noise,coefficients(2)];
            end
        end
    end
end
