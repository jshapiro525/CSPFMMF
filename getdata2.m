function [DATA,parangs,lams] = getdata2(trunc,imdim)

parangs=[];
DATA=[];
lams =[];

for i = 0:34
    

       
        name='injected';
        exten='.fits';
        filenumb=i;
        numb=num2str(filenumb);

        filename=strcat(numb,name,exten);
        
        thisparang = fitsheader(filename,'PAR_ANG');
        parangs =[parangs thisparang];

        data=fitsread(filename,'image');
        if trunc      
            dim=max(size(data));
            c=144;
            r=144;
            w=floor(imdim/2);
            data=data(r-w:r+w,c-w:c+w,:);
        end
        DATA=cat(3,DATA,data); 
        


end
        lams = [1.49461 1.50302 1.51143 1.51984 1.52825 1.53666 1.54507 1.55348 1.56189 ...
        1.5703 1.57871 1.58712 1.59554 1.60395 1.61236 1.62077 1.62918 1.63759  ...
        1.646 1.65441 1.66282 1.67123 1.67964 1.68805 1.69646 1.70488 1.71329 ...
        1.7217 1.73011 1.73852 1.74693 1.75534 1.76375 1.77216 1.78087 1.78898 1.79739];

    
[imdim, dum1, dum2] = size(data);

for i=1:imdim
    for j=1:imdim
        for k=1:length(DATA(1,1,:))
            if isnan(DATA(i,j,k))
                DATA(i,j,k)=0;
            end
        end
    end
end
end