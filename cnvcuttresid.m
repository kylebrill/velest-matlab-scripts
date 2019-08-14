function cnvcuttresid(cnvfilein)
%CNVCUTTRESID cuts CNV based on gap, rms and residual from last VELEST run
%   Used after VELEST runs to cut out events with high residuals, rms, or
%   gap returned from previous velest iteration. Output is a new .CNV file
%   with "rev" preceeding filetype for use moving forward.
%
% Written by Kyle Brill
% Last edited by Kyle Brill on 17 September, 2016 at 17:00 EDT

%% Establish thresholds
minresid=0.25;
minsta=6;
maxgap=360;
maxrms=1;
numpws=5;

%% Read in outputs from previous velest run
[ot,ola,olo,od,mag,gap,rms,~,picks]=readcnv(cnvfilein);
outfilename=strrep(cnvfilein,'.CNV','.OUT');
[event,~,eqnum]=get_sta_res(outfilename);

%% Figure out which events are still worth looking at based on residuals
staname=fields(event);
staname=staname(2:end,1);% accounts for "origtime" field
staname=char(staname);
numsta=length(staname); 
keepevents=ones(length(event),1);

for in=numpws+1:length(event)
    chksta=0;
    for jn=1:numsta
        if ~isempty(event(in).(staname(jn,1:4)))
            if event(in).(staname(jn,1:4)).res >= minresid %pick residuals < 0.5
                event(in).(staname(jn,1:4))=[];
                rmpick=find(strcmp(picks(in).sta,staname(jn,1:4)));
                picks(in).sta(rmpick)=[];
                picks(in).obt(rmpick)=[];
                picks(in).pha(rmpick)=[];
                picks(in).wgt(rmpick)=[];
                
            else
                chksta=chksta+1;
            end
        end
    end
    if chksta<minsta
        keepevents(in,1)=0;
    end
end

%% Print Hypocenter info to file based on criteria
cnvfileout=strrep(cnvfilein,'.CNV','rev.CNV');
shotfileout=strrep(cnvfilein,'.CNV','sh.CNV');

fid=fopen(cnvfileout,'w');
accepted = numpws;
for nevt=1:numpws; % Write all PWS events
    [yr,mo,da,hr,mn,sc]=datevec(datenum(ot(nevt)));
    fprintf(fid,...
        '%02d%02d%02d %02d%02d %05.2f %7.4fN %8.4fW %7.2f  %5.2f   %3d     %5.2f\n',...
        yr-2000,mo,da,hr,mn,sc,...
        ola(nevt),-1*olo(nevt),od(nevt),...
        mag(nevt),gap(nevt),rms(nevt));
    
    % write out the pick info
    for npk=1:length(picks(nevt).obt)
        if npk==7 || npk==13
            fprintf(fid,'\n');
        end
        
        fprintf(fid,'%4s%1s%1d%6.2f',...
            picks(nevt).sta{npk},picks(nevt).pha{npk},...
            picks(nevt).wgt(npk),picks(nevt).obt(npk));
        
        if npk==length(picks(nevt).obt)
            fprintf(fid,'\n\n');
        end
    end
end
for nevt=numpws+1:eqnum % Write good eqs
    if gap(nevt) <= maxgap
        if rms(nevt) <= maxrms
            if keepevents(nevt)
% write out the hypocenter info
%(3i2,1x,2i2,1x,f5.2,1x,f7.4,a1,1x,f8.4,a1,1x,f7.2,2x,f5.2)
%a4,a1,i1,f6.2
%      840210 23 6 50.92 37.1653N 121.5500W    5.62   0.00     73  0.0 0.03  1.0  1.0
%      HSPMP0  1.31CADMP1  1.55HGSMP0  1.97JCBMP1  2.36CCOMP2  3.07HCAMP1  2.79
%      HGWMP1  3.19CAOMP1  3.62JSTMP1  3.83HCRMP1  3.83HFEMP0  3.98JHLMP2  4.44
%      HPLMP0  4.31JALMP0  4.53
% convert from epoch time to string to datenum to datevec
                
                [yr,mo,da,hr,mn,sc]=datevec(datenum(ot(nevt)));
                fprintf(fid,...
                    '%02d%02d%02d %02d%02d %05.2f %7.4fN %8.4fW %7.2f  %5.2f   %3d     %5.2f\n',...
                    yr-2000,mo,da,hr,mn,sc,...
                    ola(nevt),-1*olo(nevt),od(nevt),...
                        mag(nevt),gap(nevt),rms(nevt));
                
                accepted = accepted + 1; %Count number of events printed
                % write out the pick info
                for npk=1:length(picks(nevt).obt)
                    if npk==7 || npk==13
                        fprintf(fid,'\n');
                    end

                    fprintf(fid,'%4s%1s%1d%6.2f',...
                        picks(nevt).sta{npk},picks(nevt).pha{npk},...
                        picks(nevt).wgt(npk),picks(nevt).obt(npk));
                    
                    if npk==length(picks(nevt).obt)
                        fprintf(fid,'\n\n');
                    end
                end
            end
        end
    end
end
fclose(fid);

fid=fopen(shotfileout,'w');
accepteds = 0;
for nevt=eqnum+1:length(event)
    if gap(nevt) <= 360
        if rms(nevt) <= maxrms
            if keepevents(nevt)
                [yr,mo,da,hr,mn,sc]=datevec(datenum(ot(nevt)));
                fprintf(fid,...
                    '%02d%02d%02d %02d%02d %05.2f %7.4fN %8.4fW %7.2f  %5.2f   %3d     %5.2f\n',...
                    yr-2000,mo,da,hr,mn,sc,...
                    ola(nevt),-1*olo(nevt),od(nevt),...
                    mag(nevt),gap(nevt),rms(nevt));
        
                accepteds = accepteds + 1; %Count number of events printed
                % write out the pick info
                for npk=1:length(picks(nevt).obt)
                    if npk==7 || npk==13
                        fprintf(fid,'\n');
                    end

                    fprintf(fid,'%4s%1s%1d%6.2f',...
                        picks(nevt).sta{npk},picks(nevt).pha{npk},...
                        picks(nevt).wgt(npk),picks(nevt).obt(npk));
            
                    if npk==length(picks(nevt).obt)
                        fprintf(fid,'\n\n');
                    end
                end
            end
        end
    end
end
fclose(fid);

display([num2str(accepted) ' events in rev.CNV'])
display([num2str(accepteds) ' shots in sh.CNV'])