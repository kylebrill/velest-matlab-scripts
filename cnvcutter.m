function cnvcutter(cnvfilein)
%CNVCUTTER - Used after VELEST runs to trim events based on GAP and outputs
%a new .CNV for use moving forward.

[ot,ola,olo,od,mag,gap,rms,~,picks]=readcnv(cnvfilein);

%% Print Hypocenter info to file based on criteria
cnvfileout=strrep(cnvfilein,'.CNV','rev.CNV');
fid=fopen(cnvfileout,'w');
accepted = 0;
for nevt=1:length(ot)
    if gap(nevt) <= 180
        if rms(nevt) <= 1
        % write out the hypocenter info
        %(3i2,1x,2i2,1x,f5.2,1x,f7.4,a1,1x,f8.4,a1,1x,f7.2,2x,f5.2)
        %a4,a1,i1,f6.2
        %      840210 23 6 50.92 37.1653N 121.5500W    5.62   0.00     73  0.0 0.03  1.0  1.0
        %      HSPMP0  1.31CADMP1  1.55HGSMP0  1.97JCBMP1  2.36CCOMP2  3.07HCAMP1  2.79
        %      HGWMP1  3.19CAOMP1  3.62JSTMP1  3.83HCRMP1  3.83HFEMP0  3.98JHLMP2  4.44
        %      HPLMP0  4.31JALMP0  4.53
        % convert from epoch time to string to datenum to datevec
        
        [yr,mo,da,hr,mn,sc]=datevec(datenum(ot(nevt)));
        fprintf(fid,'%02d%02d%02d %02d%02d %05.2f %7.4fN %8.4fW %7.2f  %5.2f   %3d     %5.2f\n',yr-2000,mo,da,hr,mn,sc,ola(nevt),-1*olo(nevt),od(nevt),mag(nevt),gap(nevt),rms(nevt));
        accepted = accepted + 1; %Count number of events printed
        % write out the pick info
        for npk=1:length(picks(nevt).obt)
            if npk==7 || npk==13
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
end
fclose(fid);
display([num2str(accepted) ' events in new .CNV'])