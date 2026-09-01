%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regrid CM1 output files, taking only what is needed for the hail
%trajectory model. Originally by MK ca. 2019; modified 9/11/25 by MK with
%ideas from Lydia Spychalla and Yuzhu Lin. For ECSS Shear study by K.
%Lombardo et al.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

%-----  User Specified Parameters  -----------------------------------------
altitude_select = 6.5;    %select altitude at which to track updraft (km); usually have used 6-7 km
compsize = 50;            %size of half of composite grid (elements; half of domain used here)
                          %so, e.g., 50 means grid is ~100 x 100 elements

start_file_num = 13;      %select your starting CM1 output file number
end_file_num = 37;        %select your ending CM1 output file number

w_thresh = 15;            %select updraft threshold for searching for centroid. Use 15/20 m/s or so
y_ind_cut = 120;          %index for cutting out left mover. Guesstimated from looking at images.
%------------------------------------------ end user specified params---------

%To make sure you have the right sizes allocated, read in a test file to
%get the grid parameters. Modified for CM1 ca. 2025 (on Apple)
testfile = 'cm1out_000012.nc';
xgrid = double(ncread(testfile,'xh')); delx = mode(diff(xgrid)); %specify grid spacing, or could read it in explicitly
ygrid = double(ncread(testfile,'yh')); dely = mode(diff(ygrid));
zgrid = double(ncread(testfile,'zh')); delz = mode(diff(ygrid)); %specify grid spacing, or could read it in explicitly

%Computes the vertical grid index closest to your selected altitude for
%updraft tracking
find_indz = abs(zgrid-altitude_select);
zpick =find(find_indz == min(find_indz));  

%this loops over all the cm1 output files in the directory -- was hardcoded
%before I learned how to automatically do it.
filen_list = [start_file_num:end_file_num];
num_f = numel(filen_list);
        %Could do this for automatically going through all files in directory
        %files = dir('cm1out_0000*.nc');
        %numfiles = numel(files);
findex = 0;
w_area = zeros(num_f,1); %allocate array for updraft area exceeding threshold specified above (km^2)
x_loc = zeros(num_f,1);  %allocate array for updraft centroid x location (km)
y_loc = x_loc;           %allocate array for updraft centroid y location (km)
w90 = x_loc;             %allocate array for updraft core 90th pecentile (m/s)

%Now loop over the files to extract data
for xx=1:num_f
    findex = findex +1;                    %file index counter
    filen = ['cm1out_0000',num2str(filen_list(xx)),'.nc']; %construct filename
    
        %THESE ARE THE QUANTITIES WE USE FOR THE HAIL MODEL
     w = double(ncread(filen,'winterp'));  %w wind field (m/s)
     u = double(ncread(filen,'uinterp'));  %u wind field (m/s)
     v = double(ncread(filen,'vinterp'));  %v wind field (m/s)
     dbz = double(ncread(filen,'dbz'));    %simulated reflectivity (dBz)
     ncr = double(ncread(filen,'ncr'));    %num mixing ratio rain (kg^-1)
     qr = double(ncread(filen,'qr'));      %rain mass mixing ratio (kg/kg)
     ncs = double(ncread(filen,'ncs'));    %num mixing ratio snow (/kg)
     qs = double(ncread(filen,'qs'));      %snow mass mixing ratio (kg/kg)
     nci = double(ncread(filen,'nci'));    %ice num mixing ratio (/kg)
     qi = double(ncread(filen,'qi'));      %ice mass mixing ratio (kg/kg)
     ncg = double(ncread(filen,'ncg'));    %graupel number mixing ratio (/kg)
     qg = double(ncread(filen,'qg'));      %graupel mass mixing ratio (kg/kg)
     qc = double(ncread(filen,'qc'));      %cloud mass mixing ratio (kg/kg)
     qv = double(ncread(filen,'qv'));      %vapor mixing ratio (kg/kg)
     prs = double(ncread(filen,'prs'));    %pressure field (Pa)
     th = double(ncread(filen,'th'));      %potential temperature field (K)
     rho_dry = double(ncread(filen,'rho'));%dry air density (kg/m3)

      %extract horizontal slice through domain at selected altitude
      wslice = w(:,:,zpick);   %extracts 2D horizontal cross section of w at desired altitude
      si = size(wslice);       %gets the size of the 3D arrays
        
             %original code looked for max w
             %tempp=find(wslice == max(max(wslice)));  %finds the max w value; note, 
             %                                         %could do something smarter like my students did, 
             %                                         %like the geometric center. Max
                                                       %is a bit volatile in location
      %The next few lines identifies the updraft core centroid, with indices a, b. 
      %It does so by first finding the core (w>threshold specified above); 
      %calculating the 90th percentile speed within that core; 
      %identifying all grid points with w exceeding that 90th percentile;
      %determining the median x, y coords of that subset (called a, b).
      
      poo = find(wslice>=w_thresh);        %find indices where updraft at this altitude exceeds specified
                                            %threshold
      
      poop=quantile(wslice(poo),[.9 .95]); %find 90th and 95th percentile speed of updraft contained 
                                           % within specified threshold contour
      poop2=find(wslice>poop(1));          %here, find indices in cross section where w exceeds this 90th percentile
      
      [ii, jj] = ind2sub(si,poop2);        %Return the indices of where those points are located
      
       cut_lm = find(jj <= y_ind_cut);           %to cut out strong updraft from left mover 
       ii = ii(cut_lm); jj = jj(cut_lm);         %remove any points found northward of line specified above

      a = floor(median(ii)); b = floor(median(jj));      %calculate the median index of that subset -- used as updraft CENTROID
                                                         %added "floor" because could be a xx.5 value
 
        %THESE ARE THE QUANTITIES THAT ARE WRITTEN OUT AND USED AS INPUT TO THE
        %HAIL MODEL
         comp_w = w(a-compsize:a+compsize,b-compsize:b+compsize,:);
         comp_u = u(a-compsize:a+compsize,b-compsize:b+compsize,:);
         comp_v = v(a-compsize:a+compsize,b-compsize:b+compsize,:);
         comp_qc = qc(a-compsize:a+compsize,b-compsize:b+compsize,:);
         comp_qv = qv(a-compsize:a+compsize,b-compsize:b+compsize,:);
         comp_nr = ncr(a-compsize:a+compsize,b-compsize:b+compsize,:);
         comp_qr = qr(a-compsize:a+compsize,b-compsize:b+compsize,:);
         comp_ni = nci(a-compsize:a+compsize,b-compsize:b+compsize,:);
         comp_qi = qi(a-compsize:a+compsize,b-compsize:b+compsize,:);
         comp_ns = ncs(a-compsize:a+compsize,b-compsize:b+compsize,:);
         comp_qs = qs(a-compsize:a+compsize,b-compsize:b+compsize,:);
         comp_ng = ncg(a-compsize:a+compsize,b-compsize:b+compsize,:);
         comp_qg = qg(a-compsize:a+compsize,b-compsize:b+compsize,:);
         comp_dbz = dbz(a-compsize:a+compsize,b-compsize:b+compsize,:);
         comp_prs = prs(a-compsize:a+compsize,b-compsize:b+compsize,:);
         comp_theta = th(a-compsize:a+compsize,b-compsize:b+compsize,:);
         comp_rhod = rho_dry(a-compsize:a+compsize,b-compsize:b+compsize,:);

         %a few more quantities of interest for updraft properties
         w_area(findex) = numel(poo) .* delx .* dely;   %calculate updraft area (km^2)
         x_loc(findex) = xgrid(a);                      %save off x location (km) of updraft centroid
         y_loc(findex) = ygrid(b);                      %save off y location (km) of updraft centroid
         w90(findex) = poop(1);                         %save off 90th percentile core speed m/s
        
        savnam = ['CAPE_2500_NOMOD_',num2str(filen_list(xx)),'.mat'];
        save(savnam,'comp_u','comp_v','comp_w','comp_qc','comp_qv','comp_nr','comp_qr',...
            'comp_ni','comp_qi','comp_ns','comp_qs','comp_ng','comp_qg','comp_dbz','comp_prs','comp_theta',...
            'comp_rhod','xgrid','ygrid','zgrid');

        xx
end   


%Can compute storm motion in a number of ways, lets do an individual one
%for each timestep here. ASSUMES OUTPUT EVERY 5 MINUTES!

 est_u = gradient(x_loc).*1000./(5*60);
 est_v = gradient(y_loc).*1000./(5*60);
save('timeseries_quants.mat','w_area','x_loc','y_loc','w90','est_u','est_v')
