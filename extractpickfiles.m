% write pick file from antelope db
clear, close all
db = dbopen('VBlocal','r');
% dbarrival=dblookup_table(db,'arrival');
dbor =dblookup_table(db,'origin');
% dbevent = dblookup_table(db,'event');
% dbor = dbjoin(dbor,dbevent);
% dbfree(dbevent)

dbor_size = dbquery(dbor,'dbRECORD_COUNT');
disp(['sorting through ',num2str(dbor_size),' events'])

% origerr
dborerr = dblookup_table(db,'origerr');
% nrecords = dbquery(dborerr,'dbRECORD_COUNT')
orids=dbgetv(dbor,'orid');
% dbfree(dborerr);


% origin x assoc x arrival
dbassoc1 = dblookup_table(db,'assoc');
% phase=dbgetv(dbassoc1,'phase')
dbassoc = dbjoin(dbor,dbassoc1,'orid');

dbarr1 = dblookup_table(db,'arrival');
dbarr = dbjoin(dbassoc,dbarr1);
dbarr = dbsort(dbarr,'orid');
dbarr_size = dbquery(dbarr,'dbRECORD_COUNT')

[oo,ola,olo,ot,od]=dbgetv(dbor,'orid','lat','lon','time','depth');
[ao,at,as,ap,ach,aperr]=dbgetv(dbarr,'origin.orid','arrival.time','arrival.sta','arrival.iphase','arrival.chan','arrival.deltim');

%%
fid=fopen('testhypo.txt','w')
a=find(oo>0);
%  oridstmp=orids(a)
for nevt=1:length(orids)
    if oo(nevt)>0
        tmp=find(oo(nevt)==ao);
        
        % write out the hypocenter info
        %(3i2,1x,2i2,1x,f5.2,1x,f7.4,a1,1x,f8.4,a1,1x,f7.2,2x,f5.2)
        %a4,a1,i1,f6.2
        %      840210 23 6 50.92 37.1653N 121.5500W    5.62   0.00     73  0.0 0.03  1.0  1.0
        %      HSPMP0  1.31CADMP1  1.55HGSMP0  1.97JCBMP1  2.36CCOMP2  3.07HCAMP1  2.79
        %      HGWMP1  3.19CAOMP1  3.62JSTMP1  3.83HCRMP1  3.83HFEMP0  3.98JHLMP2  4.44
        %      HPLMP0  4.31JALMP0  4.53
        % convert from epoch time to string to datenum to datevec
        [yr,mo,da,hr,mn,sc]=datevec(datenum(strtime(ot(nevt))));
        odatenum=datenum(strtime(ot(nevt)));
        fprintf(fid,'%02d%02d%02d %02d%02d %05.2f %7.4fN %8.4fW %7.2f   0.00   360  9.9 9.99  9.9  9.9\n',yr-2000,mo,da,hr,mn,sc,ola(nevt),-1*olo(nevt),od(nevt));
        
        % write out the pick info
        for npk=1:length(tmp)
            wt=[]; diffsec=[];
            % adjust pick time
            [yrp,mop,dap,hrp,mnp,scp]=datevec(datenum(strtime(at(tmp(npk)))));
            pdatenum=datenum(strtime(at(tmp(npk))));
            diffsec(npk)=86400*(pdatenum-odatenum);
            % adjust pick weight
            if aperr(tmp(npk))<=0.1
                wt(npk)=0;
            elseif aperr(tmp(npk))<=.2
                wt(npk)=1;
            elseif aperr(tmp(npk))<=.3
                wt(npk)=2;
            elseif aperr(tmp(npk))<=.4
                wt(npk)=3;
            else
                wt(npk)=4;
            end
            phase=char(ap(tmp(npk)));
            tmpsta=char(as(tmp(npk)));
            if length(tmpsta)>4, sta(1:4)=tmpsta(1:4);end
            if length(tmpsta)==3, sta(1:4)='____';sta(1:3)=tmpsta(1:3);end
            if length(tmpsta)==4, sta(1:4)=tmpsta(1:4);end
            if npk==7 || npk==13
                fprintf(fid,'\n');
            end
            sprintf('%4s%1s%1d%6.2f',sta,phase,wt(npk),diffsec(npk))
            fprintf(fid,'%4s%1s%1d%6.2f',sta,phase,wt(npk),diffsec(npk));
            if npk==length(tmp)
                fprintf(fid,'\n');
            end
        end
        clear sta phase wt diffsec
    else
        %bogus orid
    end
end
fclose(fid);