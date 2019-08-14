function loc2cnv(locfile)
%LOC2CNV - print a .CNV file from a loc file from VELEST single event mode.
%   LOC2CNV outputs a .CNV file to use as a new input for velest while
%   trying to find a minimum 1D model. Currently, cuts events based on RMS
%   of event locations from single event mode, but this can be adjusted on
%   line 47 of this function for other criteria.
%
% Written and last edited by Kyle Brill on 30/06/2016 at 13:00 EDT
eventinfo=readseloc(locfile);
cnvfile=strrep(locfile,'.loc','.CNV');
%account for deleted events during single event mode relocation
deleter=arrayfun(@(s) ~isempty(s.origin),eventinfo);
eventinfo=eventinfo(deleter);

ot=[eventinfo.origin]';
ola=[eventinfo.lat]';
olo=[eventinfo.lon]';
od=[eventinfo.depth]';
mag=[eventinfo.mag]';
gap=[eventinfo.gap]';
rms=[eventinfo.rms]';

picks(length(eventinfo),1).sta=fieldnames(eventinfo(length(eventinfo),1).picks);
picks(length(eventinfo),1).pha=cell(size(picks(length(eventinfo),1).sta));
picks(length(eventinfo),1).pha(:)={'P'};
picks(length(eventinfo),1).wgt=zeros(length(picks(length(eventinfo),1).pha),1);
picks(length(eventinfo),1).obt=picks(length(eventinfo),1).wgt;
for wnt=1:length(picks(length(eventinfo),1).pha)
    picks(length(eventinfo),1).wgt(wnt,1)=eventinfo(length(eventinfo),1).picks.(picks(length(eventinfo),1).sta{wnt,1}).prmk;
    picks(length(eventinfo),1).obt(wnt,1)=eventinfo(length(eventinfo),1).picks.(picks(length(eventinfo),1).sta{wnt,1}).tpobs;
end

for in=1:length(eventinfo)-1
    picks(in,1).sta=fieldnames(eventinfo(in,1).picks);
    picks(in,1).pha=cell(size(picks(in,1).sta));
    picks(in,1).pha(:)={'P'};
    picks(in,1).wgt=zeros(length(picks(in,1).pha),1);
    picks(in,1).obt=picks(in,1).wgt;
    for jn=1:length(picks(in,1).pha)
       picks(in,1).wgt(jn,1)=eventinfo(in,1).picks.(picks(in,1).sta{jn,1}).prmk;
       picks(in,1).obt(jn,1)=eventinfo(in,1).picks.(picks(in,1).sta{jn,1}).tpobs;
    end
end

fid=fopen(cnvfile,'w');
accepted = 0;
for nevt=1:length(ot)
    if rms(nevt,1) < 1 %&& gap(nevt,1) < 180  % This can be adjusted to other criteria
        % write out the hypocenter info
        %(3i2,1x,2i2,1x,f5.2,1x,f7.4,a1,1x,f8.4,a1,1x,f7.2,2x,f5.2)
        %a4,a1,i1,f6.2
        %      840210 23 6 50.92 37.1653N 121.5500W    5.62   0.00     73  0.0 0.03  1.0  1.0
        %      HSPMP0  1.31CADMP1  1.55HGSMP0  1.97JCBMP1  2.36CCOMP2  3.07HCAMP1  2.79
        %      HGWMP1  3.19CAOMP1  3.62JSTMP1  3.83HCRMP1  3.83HFEMP0  3.98JHLMP2  4.44
        %      HPLMP0  4.31JALMP0  4.53
        % convert from epoch time to string to datenum to datevec
        
        [yr,mo,da,hr,mn,sc]=datevec(datenum(ot(nevt)));
        fprintf(fid,'%02d%02d%02d %02d%02d %05.2f %7.4fN %8.4fW %7.2f  %5.2f   %3d      %4.2f\n',yr-2000,mo,da,hr,mn,sc,ola(nevt),-1*olo(nevt),od(nevt),mag(nevt),gap(nevt),rms(nevt));
        accepted = accepted + 1; %Count number of events printed
        % write out the pick info
        for npk=1:length(picks(nevt).obt)
            if npk==7 || npk==13 || npk ==19
                fprintf(fid,'\n');
            end
            %                         sprintf('%4s%1s%1d%6.2f',sta,phase,wt(npk),diffsec(npk))
            %                         fprintf(fid,'%4s%1s%1d%6.2f',sta,phase,wt(npk),diffsec(npk));
            %                         sprintf('%4s%1s%1d%6.2f',sta,phase,wt,diffsec)
            fprintf(fid,'%4s%1s%1d%6.2f',picks(nevt).sta{npk},picks(nevt).pha{npk},picks(nevt).wgt(npk),picks(nevt).obt(npk));
            
            if npk==length(picks(nevt).obt)
                fprintf(fid,'\n\n');
            end
        end
    end
end
fclose(fid);

display([num2str(accepted) ' events in new .CNV'])