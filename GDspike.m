function [ grp_delay , signal] = GDspike( FOUT,sampling_rate, smoothening_factor, winScaleFactor, thres )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%GDSPIKE Algorithm: Estimates the spike positions from a calcium fluorescence signal
% Author : Jilt Sebastian & Mari Ganesh Kumar
% Last Updated on: Feb 21, 2017
%
% INPUT Arguments:
% FOUT : Calcium fluorescence signal
% sampling rate : Sampling rate of calcium signal
% smoothening factor : For moving average preprocessing: Default valueis 1  (no smoothing)
% winScaleFactor : Parameter to select the length of root cepstral signal on which group delay is computed
% thres: Threshold for spike train estimation from the spiking information
% 
% 
% OUTPUT Arguments
%
% grp_delay: Output the group delay domain representation
% signal: Spike train
%
% Contact information:
%jiltsebastian@gmail.com, mariganeshkumar@gmail.com
%This work is accepted for poster presentation  at ICASSP 2017, New%Orleans, USA
%"GDspike: An Accurate Spike Estimation Algorithm from Noisy Calcium Fluorescence Signals" 
%by Jilt Sebastian , Mari Ganesh Kumar , Y. S.Sreekar and Hema A. Murthy from IIT Madras India and  Rajeev Vijay Rikhye
%and Mriganka Sur from Suur Lab, MIT, USA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % winScaleFactor = 5:5:45;
   
                    %Part1: Calculating the GD signal
                    S = smooth(FOUT,smoothening_factor,'moving');%MA smoothing-preprocessing
                    grp_delay = ones(length(S),1);
                    gd_sum = ones(length(S),1);

                    for wsfIndex = 1:length(winScaleFactor)
                        tempDir = sprintf('temp_%d',wsfIndex);
                        warning('off','all')
                        mkdir(tempDir); cd(tempDir);
                        energy_file_name = strcat(sprintf('neuron'),'.en');
                        dlmwrite(energy_file_name,S,'\n');
                        spec_file_name = energy_file_name(1:end-2);
                        spec_file_name =strcat(spec_file_name,'spec');

                        % Invoking the binary   : contains control file for group delay computation--see the source file fe-words.base                     
                        copyfile('../fe-words.base','fe-words.base');
                        ctrl_file = 'fe-words.base';
                        temp_ctrl_file = strcat('temp.base');

                        % Changing the winscalefactor parameter in config file
                        a = importdata(ctrl_file);
                        a = struct2cell(a);
                        a{1}(3) = winScaleFactor(wsfIndex);
                        fid0 = fopen(temp_ctrl_file,'w');
                        for i = 1:length(a{1})
                            fprintf(fid0,'%s %s %f\n',char(a{2}(i,1)),char(a{2}(i,2)),a{1}(i));
                        end
                        copyfile(temp_ctrl_file,ctrl_file);                       
                        fclose(fid0);                        
                        dummy1 = 'b';
                        dummy2 = 'c';
                        dummy3 = 'd';
                        dummy4 = 'e';
                        dump = 'dump.txt';
                       %Invoking the binary file WordSegmentsWithSilenceRemoval --see the source file  WordSegmentsWithSilenceRemoval.c
                       %Input-fluorescence signal Output-GD representation
                        system(sprintf('../WordSegmentWithSilenceRemoval %s %s %s %s %s %s %s > %s 2>&1',ctrl_file,energy_file_name,spec_file_name,dummy1,dummy2,dummy3,dummy4,dump));
                          
                        delete(energy_file_name);
                        temp = load(spec_file_name);%spec_file with GD domain signal
                        temp = temp(:,5);
                        temp(length(S)+1:end) = [];
                        grp_delay = grp_delay.*temp;
                        temp = temp - mean(temp);
                        gd_sum = gd_sum + cumsum(temp);
                        cd ..; 
                    end

                    grp_delay = diff(gd_sum);
                    grp_delay = grp_delay/max(grp_delay); %Normalization step
                    grp_delay=[grp_delay(1:end-20);grp_delay(end-20).*ones(20,1)];
                    assignin('base','grp_delay',grp_delay);
                   
                    %======================================================================
                    % Part2: Reading the contents of group delay file, and getting the
                    % spike positions
                    threshold = thres(indexD);

                    stroke_loc = zeros(1,length(grp_delay));
                    % Go to each minima, and calculate height till next maxima. Keep a
                    % threshold on this to decide if stroke!


                    t = 1:length(grp_delay);
                    [ymax,imax,ymin,imin] = extrema(grp_delay);

                    % sort the minimas and maximas;
%                     imin;
%                     ymin;
                    temp_min = sortrows([imin ymin]);
                    imin = temp_min(:,1)';
                    ymin = temp_min(:,2)';
                    clear temp_min;

                    temp_max = sortrows([imax ymax]);
                    imax = temp_max(:,1)';
                    ymax = temp_max(:,2)';
                    clear temp_max;



                    if (imin(1) < imax(1) )  % fine, just truncate the maximum
                        imin(1) = []; ymin(1) = [];

                        if (length(imin) > length(imax) )
                            imin(length(imax)+1:end) = [];
                            ymin(length(imax)+1:end) = [];
                        elseif (length(imin) < length(imax) )
                            imax(length(imin)+1:end) = [];
                            ymax(length(imin)+1:end) = [];
                        end
                    else                                                    

                            if (length(imin) > length(imax) )
                            disp('this shouldnt have come');
                            imin(length(imax)+1:end) = [];
                            ymin(length(imax)+1:end) = [];
                        elseif (length(imin) < length(imax) )

                            imax(length(imin)+1:end) = [];
                            ymax(length(imin)+1:end) = [];
                        end
                    end


                    assignin('base','ymax',ymax);
                    assignin('base','imax',imax);
                    assignin('base','ymin',ymin);
                    assignin('base','imin',imin);
                    assignin('base','grp_delay',grp_delay);


                     %==================================================================
                     % Algorithm1  for stroke location
                     index_stroke = 1;
                     peak_valley_heights = ymax - ymin;
                     peak_valley_heights = peak_valley_heights(1:length(peak_valley_heights));

                     for index = 1:1:length(peak_valley_heights)
                         if (peak_valley_heights(index) > threshold)
                             %stroke_loc(index_stroke) = ceil((imin(index) + imax(index))/2);
                             stroke_loc(index_stroke) = imin(index) ;
                             index_stroke = index_stroke + 2;
                         end
                     end
                     %==================================================================


                    stroke_loc(stroke_loc==0) = [];

                    assignin('base','stroke_loc',stroke_loc);
                    assignin('base','peaks',peaks);

                    %======================================================================


                    % Part 3: Traingulation Step
                    %======================================================================

                    X=FOUT;
                    Fs=sampling_rate;
                    time=length(S)/Fs;% Converting into seconds
                    isort = sort([imax imin]);
                    exts=0;
                    exte=0;
                    iexts=1;
                    iexte=1;
                    ind=1;


                    for ind= 1:length(isort)-1
                        iexts=isort(ind);
                        iexte=isort(ind+1);
                        exts = grp_delay(iexts);
                        exte = grp_delay(iexte);
                        valley_lenght=abs(exts-exte);
                        signal(iexts)=0;
                        signal(iexte)=0;
                        intermediate_samples=iexte-iexts-2;
                        increment_value=valley_lenght/(floor((intermediate_samples+1)/2));
                        sum=0;
                        if valley_lenght>0.1*max(grp_delay)
                            for j=iexts+1:iexts+floor((intermediate_samples+1)/2)
                                sum=sum+increment_value;
                                signal(j)=sum;
                            end 
                            sum=0;
                            for j=iexte-1:-1:iexts+floor((intermediate_samples+1)/2)+1
                                sum=sum+increment_value;
                                signal(j)=sum;
                            end
                        else
                            signal(iexts:iexte)=0;
                        end
                    end
                    signal(isort(length(isort)):time*Fs)=0;
             
          
end

