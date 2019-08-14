% use to run 1000 random velest models
tic
setenv('DYLD_LIBRARY_PATH', '/home/campus23/kabrill/bin/');

%Premake all velocity models and write to file
%outfile0='mnum';
outdir='mnum/';

% fix bottom layer to submoho vel at 8 km/s?
% next to bottom between 6.5 and 7.5
nmodels=1000;
% l8=6.5 + (7.5-6.5)*rand(nmodels,1);
% l1=1.5 + (4.0-1.5)*rand(nmodels,1);%2.5
%KAB-EDIT fix bottom layers to Molina and Tenorio 2000 referenced in Franco
%et. al 2009
l6=2.0 + (6.50-2.0)*rand(nmodels,1);
l1=0.60 + (4.00-0.6)*rand(nmodels,1);

m=zeros(9,nmodels);
%array of depths
%depths=[-3.0,0.0,2.0,4.0,6.0,12.0,18.0,34.0];
depths=[-2,0.0,1.5,3.0,5.0,7.0,9.0,17.0,37.0];
% depths and velocities from Molina and Tenorio 2000 referenced in Franco
% et al. 2009
fprintf(1,'Preparing Model Number     ');
for n=1:nmodels;
    
    %Make a directory for this model and all velest output files
    outdir1=strrep(outdir,'num',num2str(n));
    
    modliststr=cat(2,'velocitymods/',outdir1);
    if isempty(dir(modliststr))%only make folder if it doesn't exist yet
        mkdir(modliststr)
    end
    
    copyfile('velest.cmn',cat(2,modliststr,'velest.cmn'));
    
    inputname=cat(2,modliststr,'fgpws-rand01.mod');
    
    %actualfile=dir(inputname);
    %actualfile=inputname(15:23);
    fid=fopen(inputname,'w');
    
    %fid=fopen(strrep(outfile0,'num',num2str(n)),'w');
    m(1,n)=l1(n);
    m(6,n)=l6(n);
    %now select the 3 km layer (layer 4)
    m(4,n)=m(1,n)+(m(6,n)-m(1,n))*rand(1);
    
    % now the shallowest layers
    m(2:3,n)=sort(m(1,n)+(m(4,n)-m(1,n))*rand(2,1));
    
    % and lyer 5
    m(5,n)=m(4,n)+(m(6,n)-m(4,n))*rand(1);
    
    %     m(9,n)=8.04;
    m(7,n)=6.55;
    m(8,n)=6.75;
    m(9,n)=7.95;
    
    %Print to file
    
    %        dlmwrite(inputname,'PACAYA1D-model (mod1.1 FL090115)   Ref. station PS09','delimiter','','newline','unix')
    %        dlmwrite(inputname,'8        vel,depth,vdamp,phase (f5.2,5x,f7.2,2x,f7.3,3x,a1)','delimiter','','-append','newline','unix');
    fprintf(fid,' %s\n%s\n','FUEGO1D-model (mod1.1 EK280993)   Ref. station NW1_','9        vel,depth,vdamp,phase (f5.2,5x,f7.2,2x,f7.3,3x,a1)');
    fprintf(fid,' %4.2f     %7.2f   %06.2f            %s\n',m(1,n),depths(1),001.00,'P-VELOCITY MODEL');
    %     for i=2:8
    for i=2:6
        fprintf(fid,' %4.2f     %7.2f   %06.2f\n',m(i,n),depths(i),001.00);
    end
    
    fprintf(fid,' %4.2f     %7.2f   %06.2f\n',m(7,n),depths(7),999.00);
    fprintf(fid,' %4.2f     %7.2f   %06.2f\n',m(8,n),depths(8),999.00);
    fprintf(fid,' %4.2f     %7.2f   %06.2f\n',m(9,n),depths(9),999.00);
    
    fclose(fid);
    % Print Status
    fprintf(1,'\b\b\b\b%4.4d',n); pause(.000001)
end
fprintf('\n')

parfor n=1:nmodels;
    
    %Make a directory for this model and all velest output files
    outdir1=strrep(outdir,'num',num2str(n));
    
    modliststr=cat(2,'velocitymods/',outdir1);
    if isempty(dir(modliststr))%only make folder if it doesn't exist yet
        mkdir(modliststr)
    end
    
    %Navigate to this directory
    cd(modliststr);
    
    %Run velest
    system('/home/campus23/kabrill/bin/velest >/dev/null');
    
    %!velest
    cd('../../');
    
    fprintf([num2str(n) ' of ' num2str(nmodels) ' models complete.\n']);
    % system(cat(2,'cd ',' ./', modliststr,'/'))
    
    %[status cmdout]=system(cat(2,'cp ~/Desktop/src/velest/pacaya_tests/velest.cmn ','~/Desktop/src/velest/pacaya_tests/',modliststr,'./','velest.cmn'))
    %[status cmdout]=system(cat(2,'~/bin/velest ', '~/Desktop/src/velest/pacaya_tests/',modliststr,'/velest.cmn'))
    
    
end
toc
%lazy plot of the approx middle of each layer
%plot(m,[-1, 1, 3, 5, 10, 15, 25, 40]); axis ij