function [rep_parangs, rep_lams] = repparlam(parangs,lams)

n=length(lams);
dur=length(parangs);

rep_lams=[];
rep_parangs=[];

for i=1:dur
    rep_lams = [rep_lams lams];
end

for i=1:dur
    rep_parangs = [rep_parangs ones(1,length(lams))*parangs(i)];
end


end