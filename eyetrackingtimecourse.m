function eyetrackingtimecourse(xlsx)
%{ 
The function of this script is to take as input a .xlsx spreadsheet output 
- converted in Excel to .xlsx from the Eyelink .xls 'fixation report' - and to calculate a viewing timecourse to various ROIs. 

Any questions/comments/necessary edits can be sent to Lauren 'I'm so great at coding' Hopkins at lsh5968@gmail.com
However, be aware if you don't read the documentation and email me stupid things that the docu explains I will be pissed. 
This thing is so heavily commented my grandma can probably use it.

Usage: eyetrackingtimecourse('your-eyetracking-filename.xlsx')
e.g. eyetrackingtimecourse('s18_fixation_data.xlsx')

There are basically 3 VERY important prerequisities to using this script:
PREREQ 1. The fixation report that YOU MAKE YOURSELF in DataViewer MUST include certain variables which are listed below but which can be in any order:
    1) TRIAL_START_TIME
    2) IP_START_TIME
    3) CURRENT_FIX_INTEREST_AREA
    4) CURRENT_FIX_START
    5) CURRENT_FIX_END
    6) CURRENT_FIX_DURATION
    7) CURRENT_FIX_INTEREST_AREAS

PREREQ 2: After making your fixation report, DataViewer will output it as
an .xls file. You MUST open it and re-save/convert it to an .xlsx file. The
reason why is specified in the code below but long story short, Macs &
Linuxes and Windows computers without Office installed cannot read .xls
files with the matlab commands used to import mixed variable data.

PREREQ 3: VERY VERY VERY IMPORTANT - REFERS TO DATA SETUP IN DATAVIEWER.
This current version assumes that, if you have a background
ROI (and/or an ROI that encompasses other ROIs) that that ROI is the FINAL
one you made. So for example, if you have 6 ROIs including the bg, the bg
will be ROI #6 (as opposed to ROI #1, #2, etc.) and will be output this way in the
CURRENT_FIX_INTEREST_AREAS variable of the fixation report. It is NECESSARY
for you to set up your data this way in order to get an accurate output
here otherwise all fixations in the script output will be bg fixations.

One potential issue referring to ROIs is that, if there is a participant
who NEVER fixates in your last ROI(s), it/they will not be included in the output and the 'off screen ROI' will shift.
Optimally, if you have a bg ROI as your last ROI and it encompasses the
whole screen, this will rarely be a problem. But it might happen 1/100
times so just be aware. What this means is, if you have 6 ROIs (bg = ROI6)
and the subject never fixates the bg, you will still get accurate binning
for ROI 1-5 but you will get no output for the bg and the off-screen ROI
will be saved as ROI6 instead of the bg. Similarly, if the subject never
fixates ROI 5 OR 6, ROIs 1-4 will be fine but there will be no output for
ROI 5 & 6 and the off screen ROI will be stored as ROI5. However, if
someone doesn't fixate to the middle ROI(s) (e.g. no fixations to ROI3) but
DOES fixate to the last one, the output will be normal. AS LONG AS THERE IS
A FIXATION TO THE LAST ROI THIS PROBLEM WILL NOT OCCUR. Just be aware of
your data so you don't end up accidently recording off screen fixations as
ROI fixations.

This script DOES NOT break up trials by 'trial type' (e.g. familiar vs.
novel face trials or CS+ vs. CS- trials) because trial types are specified
by uniquely named variables defined by the user in their individual
experiment builder files (read. there are no precompiled SR Research
variables to represent trial types).
So unless everyone establishes a standard for
naming trial types - which is absurd given that every experiment will be
different in this regard - you'll have to do that yourself in Excel. I
can't do everything for you.

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%SET STATIC VARIABLES HERE TO CORRESPOND TO YOUR EXPERIMENT%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%For this to work, first just change the following 3 values to correspond to the kind of timecourse you want.
%%binsize = size of bins in ms (e.g. 250 means bin1 will be 0-250ms, bin2 will be 250-500ms, etc.
%%trialseconds = Length of your interest period in seconds (e.g. 16 for a 16-second trial)
%%filename = What you want your output file to be named. Recommended having subject identifier in there somewhere
binsize = 250;  %In ms
trialseconds = 16; %In s; When adapting this for your own purposes just replace the first value with your trial time in seconds
filename = 's18' %What you want your output file to be named
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Now STAY OUT OF THE REST OF THIS %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

triallength = trialseconds * 1000; %Converts trialseconds to ms for matchup with spreadsheet timestamps
timepoints = [0:binsize:triallength];
numbins = length(timepoints);

%%THE INPUT FILE MUST BE AN .XLSX FORMATTED FILE. A .XLS file will NOT
%%work for the reason specified in the next paragraph. There is truly no way to get around this because the .xls file has
%%both number and text data and there is no standard for how many variables
%%each worksheet will contain and in what order. Literally just save it as
%%.xlsx. That's the easiest I can make it for you.

%%IMPORTANT NOTE: The xlsxinfo/xlsread/xlswrite commands require
%%(apparently) a working version of Excel but like...specifically a Windows
%%version? Any xls... command will work in Windows and the xlsread will
%%work on other computers with Excel variants (e.g. My mac has Mac versions
%%of Office and it works on those) but if there is no communication with
%%the Excel server these commands won't work. Again, because of the
%%structure of the spreadsheet output and no universality for output file
%%content/format/organization, there is no way around this.
%%Therefore, if you get a loading error, resave your file as a
%%NON-READ-ONLY .xlsx formatted spreadsheet.

%%Raw gives you the whole spreadsheet but in a cell array. 
if isempty(xlsfinfo(xlsx))
    error('FATAL ERROR: Your input file is not an .XLSX format. It is likely still a regular .XLS file. Save it as .XLSX in Excel and try again.');    
else
    [~, ~, raw] = xlsread(xlsx);    
end

%Calls function for data extraction
[rawtrialnumber ipstarttime, trialstarttime, fixstarttime, fixdurationtime, rois] = dataextraction(raw);

%Establish trial anchor point
stimonset = ipstarttime - trialstarttime;
%How many ROIs are there?
%THIS MIGHT BE A PROBLEM IF THEY NEVER HIT THE MAX NUM ROI - LOOK INTO IT LATER
numrois = max(rois);

%%Set up trial numbers
%%Remember - Sometimes trial data will not start at 1 so we have to make it
%%so it does for our loops. I'm going to do this a really weird way.
[trialnumber] = organizetrialnumbers(rawtrialnumber);

%This is a problem if the last trial is a bad one (read. no fixations on
%last trial) but there is really no getting around this unless you explicitly specific trial number. As long as you
%know your data, if your 'numtrials' variable is less than how many
%trials you KNOW you have, you'll know the end trials are bad. You should
%probably know that anyway.
numtrials = max(trialnumber);
%Start triallabel at 0 so it will be the same size as the output INCLUDING
%the bin header so we can eventually concatenate everything for output (so
%trial labels will be e.g. 0-13 instead of 1-13 but the 0 won't mean
%anything; will only be where bin labels are)
triallabel = [0:1:numtrials]';

%Establish when each fixation begins relative to the trial onset
%IMPORTANT: Some of the values will be negative which means they were
%already fixated there prior to the trial start (i.e. in the center for
%boundry trigger). These are reset to '0'
fixonset = fixstarttime - stimonset;
for j = 1:length(fixonset)
    if fixonset(j,1) < 0
        adjustedfixonset(j,1) = 0;
        adjustedfixduration(j,1) = fixdurationtime(j,1) - abs(fixonset(j,1));
    else
        adjustedfixonset(j,1) = fixonset(j,1);
        adjustedfixduration(j,1) = fixdurationtime(j,1);
    end
end

fixendtime = adjustedfixonset + adjustedfixduration;
%Aggregate data input matrix
dataarray = [trialnumber fixonset adjustedfixonset fixendtime adjustedfixduration rois];

%%The output matrix will be of size trialXbinXroi
outputmatrix = zeros(numtrials, numbins, numrois);
binneddata = outputmatrix;

trialcounter = 1;
bincounter = 1;
%Bin it
for h = 1:numtrials
    %Individual trial extraction
    [individualtrialarray] = trialextraction(dataarray, timepoints, h);
    %%Let's establish the format right now: COL1 = trialnumber
    %%COL2 = Actual time onset; COL3 = Relative Onset
    %%COL4 = Relative Offset; COL5 = Duration; COL6 = ROI

    %First check to see if it's a 'padded' trial (e.g. bad trial w/ no fixations). If it is, fill the trial row with NaNs
    j = 1;
    if all(isnan(individualtrialarray))
        outputmatrix(h,:,:) = nan;
        binneddata(h,:,:) = nan;
    else
        for k = 1:length(individualtrialarray)
            if k == length(individualtrialarray) %FIRST make sure you're NOT at the last element otherwise the program will throw an 'out of bounds' error
                [outputmatrix, j] = lastelementcalc(individualtrialarray, binsize, timepoints, outputmatrix, k, j, h); %If it is last element, call special function 
            else
                [outputmatrix, j] = bincheck(individualtrialarray, binsize, timepoints, outputmatrix, k, j, h); %If it isn't, call regular binning function
           end
        end
        
    totaltimes = sum(outputmatrix,3);
        %Calculate proportion viewing instead of raw
        for a = 1:length(timepoints)
            for b = 1:numrois
                binneddata(h, a, b) = outputmatrix(h,a,b) / totaltimes(h,a);
            end
        end
    end
end

%%Output section
%{
Like I said above, the spreadsheet writing 'xls...' commands will not
work on certain computers - e.g. Linux, any computer without Excel installed - (read. will not universally work)
That's crap so to improve portability I have established a less organized yet all-purpose output of the regular delimited txt output variety

OUTPUT FORMAT: The script will output a single tab-delimited file with the name you
specified in the variable 'filename' at the beginning of the script.
The output file will have all the ROIs in the same sheet, one after
another - with one label above each ROI output (e.g. s18ROI:1 for ROI 1). 
The rows are trials (numbered 0-trialnumber; 0 is there to line up with the bin labels and means nothing) 
and the columns are bins (numbered 0-trial length in increments of 'binsize' so e.g. 0-16000 in 250 increments for a 16s trial with 250ms bins). 
So if you open the file in Excel, the binned
timecourse for ROI 1 will start in row 2, col 1 and will go down in rows
however many trials you have (e.g. for a 27 trial experiment, ROI 1 will
be rows 2-28). Then there will be a blank row followed by the output for
ROI 2, etc. for however many ROIs you have. Remember: You will always have
ROIs+1 matrices because the last matrix represents fixations occurring in
NONE of the ROIs (e.g. corners or anything not taken into account by bg)
%}

%%Note: If your output file already exists it will be deleted and then
%%written again. Because we're appending all our ROIs to the same file, not
%%deleting an already existing file will just append data to it, which is
%%bad if you're looking to overwrite files.
%%SO REMEMBER TO RENAME EVERY NEW FILE OR THEY WILL JUST BE OVERWRITTEN
if exist(filename, 'file') == 2
    delete(filename);
end

%%RIGHT NOW THIS DOESN'T WORK - OVERWRITES SO ONLY FINAL SHEET/ROI IS
%%WRITTEN - CHECK https://www.mathworks.com/matlabcentral/fileexchange/38591-xlwrite--generate-xls-x--files-without-excel-on-mac-linux-win
%%TO SEE IF RESOLVED OTHERWISE JUST PUT IT IN SAME SHEET w/ dlmwrite
% for a=1:numrois
%     roioutput = squeeze(binneddata(:,:,a));
%     %Add bin labels to beginning columns
%     output = [timepoints; roioutput];
%     %Then this adds the trial labels to beginning rows (with 0 where bin
%     %value labels are)
%     output = [triallabel, output];
%     xlwrite(filename, output, a);
% end

for a=1:numrois
    label = strcat(filename, ' ROI: ', num2str(a));
    label = repmat(label,size(triallabel,1), size(triallabel,2));
    roioutput = squeeze(binneddata(:,:,a));
    %Add bin labels to beginning columns
    output = [timepoints; roioutput];
    %Then this adds the trial labels to beginning rows (with 0 where bin
    %value labels are)
    output = [triallabel, output];
    %Finally, this adds the ROI and filename labels
    fid = fopen(filename,'a');
    fprintf(fid, strcat(filename, 'ROI: ', num2str(a)));
    fclose(fid);
    dlmwrite(filename, output, 'delimiter','\t', '-append', 'roffset',1);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% I) ALGORITHM FUNCTIONS %%%%%%%%%%%%%%%%%%%

%{
1) Checks to see what bin you're in. If a fixation occurs in bins OTHER
THAN the the same as or the one subsequent to the last fixation's bin, move
on to the next bin and check to see if this fixation is there instead. If
the fixation is in the same or next bin (default behavior) call the math
algorithm
%}
function [restart, j] = bincheck(individualtrialarray, binsize, timepoints, outputmatrix, k, j, h)
    if individualtrialarray(k,3) >= timepoints(1,j) && individualtrialarray(k,3) < timepoints(1,j) + binsize
        [restart, j] = binalgorithm(individualtrialarray, binsize, timepoints, outputmatrix, k, j, h);
    else
        j=j+1;
        [restart, j] = bincheck(individualtrialarray,binsize, timepoints, outputmatrix, k, j, h);
    end
end

%{
2. The actual math algorithm. Computes time spent fixated in each ROI for
each bin specified in size by you in the beginning section.
%}
function [outputmatrix, j] = binalgorithm(individualtrialarray, binsize, timepoints, outputmatrix, k, j, h)
                   %For situations in which the fixation duration is more than
                   %the binsize (e.g. > 250ms) so it spans multiple bins
                   if individualtrialarray(k,5) > binsize
                       fixtimeleft = individualtrialarray(k,5);
                       binvalue = (timepoints(1,j) + binsize) - individualtrialarray(k,3);
                       while fixtimeleft >= binsize
                           outputmatrix(h,j, individualtrialarray(k,6)) = binvalue + outputmatrix(h,j, individualtrialarray(k,6));
                           fixtimeleft = fixtimeleft - binvalue;
                           binvalue = binsize;
                           if j == length(timepoints)
                               if fixtimeleft >= binvalue
                                   fixtimeleft = binvalue;
                                   break;
                               else
                                   break;
                               end
                           else
                               j = j + 1;    
                           end
                       end
                       %Then when it gets below 250
                       outputmatrix(h, j, individualtrialarray(k,6)) = fixtimeleft;
                       %Stay in same bin, though, in case next fix is here
                       %too. New fixation calls to 'checkbin' function check this 
                   else %For fixations lengths less than binsize (e.g. < 250ms) in duration
                       if individualtrialarray(k,3) + individualtrialarray(k,5) > timepoints(1,j) + binsize
                           %IMPORTANT: Make sure you're only putting in as many ms
                           %that fit into that bin. So a 200ms fixation starting at
                           %400ms ONLY gets 100ms in bin 2 - not all 400 - and then
                           %the remaining 100ms in bin 3
                           binvalue = (timepoints(1,j)+binsize) - individualtrialarray(k,3);
                           outputmatrix(h,j, individualtrialarray(k,6)) = binvalue + outputmatrix(h,j, individualtrialarray(k,6));
                           nextbinvalue = individualtrialarray(k,5) - binvalue;                       
                           j = j+1;
                           outputmatrix(h,j, individualtrialarray(k,6)) = nextbinvalue;
                       else
                           %And the +outputmatrix part here is in case there are
                           %multiple fixations to same ROI in same bin
                           outputmatrix(h,j, individualtrialarray(k,6)) = individualtrialarray(k,5) + outputmatrix(h,j,individualtrialarray(k,6));
                           %Then check to see if next fixation is in same bin
                           if timepoints(1,j) ~= timepoints(end)
                               if individualtrialarray(k+1,3) > timepoints(1,j+1)
                                   j = j+1;
                               else
                                  j; 
                               end
                           end
                       end
                   end
end

%{
3. In order to make sure my for loops don't error - calculating viewing
time to the last fixation (i.e. last matrix element) gets its own code 
which is basically the same as the non-last element code just without the j+1 lines.
%}
function [outputmatrix, j] = lastelementcalc(individualtrialarray, binsize, timepoints, outputmatrix, k, j, h)
    last = individualtrialarray(end,:);
    lastbin = last(3);
    if lastbin < timepoints(end) && lastbin > timepoints(end-1);
        outputmatrix(h, j, last(1,6)) = timepoints(end)- last(1,3);
    else
        [outputmatrix, j] = bincheck(individualtrialarray, binsize, timepoints, outputmatrix, k, j, h);
    end
end

%%%%%%%%% II) DATA EXTRACTION FUNCTIONS %%%%%%%%%%%%%%%%

%{
1. Takes everything out of the SR Fixation Report, extracts the stuff we
need, and makes it appropriate numeric format so we can put it into one
giant matrix - the algorithm will run on these numeric vectors
%}
function [trialnumber ipstarttime, trialstarttime, fixstarttime, fixdurationtime, rois] = dataextraction(alldata)
    %Create temp storage variables to find appropriate columns we will need to
    %extract data from. This is so when outputting the spreadsheet files you
    %don't have to have the variables in the same column every time. This
    %is portability 101 and you're all welcome otherwise you'd be screwed.
    trialnumtmp = zeros(1,size(alldata,2));
    iptmp = zeros(1,size(alldata,2));
    trialtmp = zeros(1,size(alldata,2));
    roitmp = zeros(1,size(alldata,2));
    fixstarttmp = zeros(1,size(alldata,2));
    fixdurationtmp = zeros(1,size(alldata,2));

    %This finds the necessary data columns (e.g. 'IP_START_TIME' and 
    %'TRIAL_START_TIME', etc.) - whereever they are - and extracts those data 
    %into their own numeric or char variables
    for i = 1:size(alldata,2);
        trialnumtmp(1,i) = strncmpi(alldata{1,i}, 'TRIAL_INDEX', 11);
        iptmp(1,i) = strncmpi(alldata{1,i}, 'IP_START_TIME', 13);
        trialtmp(1,i) = strncmpi(alldata{1,i}, 'TRIAL_START_TIME', 16);
        fixstarttmp(1,i) = strncmpi(alldata{1,i}, 'CURRENT_FIX_START', 17);
        fixdurationtmp(1,i) = strncmpi(alldata{1,i}, 'CURRENT_FIX_DURATION', 20);
        roitmp(1,i) = strncmpi(alldata{1,i}, 'CURRENT_FIX_INTEREST_AREAS', 26);
    end
    
    %This is for user friendliness as I'm assuming you guys are gonna mess
    %this up, inevitably. This makes sure each of the variables necessary
    %for the program to run are part of your spreadsheet input. If they are
    %not - it will tell you which variables you're missing. But if you just
    %read the documentation you won't have this problem and you're an ass
    %if you do have this problem.
    trialnumber = errorcheck(alldata,trialnumtmp);
    ipstarttime = errorcheck(alldata, iptmp);
    trialstarttime = errorcheck(alldata, trialtmp);
    fixstarttime = errorcheck(alldata, fixstarttmp);
    fixdurationtime = errorcheck(alldata, fixdurationtmp);
    rois = errorcheck(alldata, roitmp);
end

%{
2. Extracts individual trial data to be run through the algorithm
%IMPORTANT: If there are any bad trials where there are no fixations, this
%function will make a 'dummy' matrix of just one row full of NaNs in order
%to pad out the output matrix so there are no errors and people can easily
%see bad trials.
%}
function trialoutput = trialextraction(alldata, timepoints, h)
    trialnumcounter = 1;
    for a = 1:length(alldata)
        if alldata(a,1) == h %Because trials is in column 1 of dataarray
           A = exist('tmparray');
           if A == 0
               tmparray(1,:) = alldata(a,:);
               trialnumcounter = trialnumcounter + 1;
           else
               tmparray(trialnumcounter, :) = alldata(a,:);
               trialnumcounter = trialnumcounter + 1;
           end    
        elseif alldata(a,1) > h %For speed
            break;
        end
    end
    %The 'bad trial' check and pad
    A = exist('tmparray');
    if A == 0
        tmparray(1,:) = nan(1, size(alldata,2));
    end
    %Some trial output will have fixations after the end time you specified (e.g.
    %at 16.5s for a 16s desired trial) so this removes those so the loops
    %don't error
    for a = size(tmparray,1):-1:1
        if tmparray(a, 3) > timepoints(end) 
            tmparray(a,:) = [];
        else
            break;
        end
    end
    trialoutput = tmparray;
end

%{
3. Purpose of this function is to check to see if all of the relevant fields
necessary for this script (as specified in the first comment at the top of this file) to run are present. If they are not it will
give an error message specifying which variable is missing.
Unfortunately the user will need to know what field the variable here
corresponds to but I don't really think that's too much to ask so deal with it.
%}

%% But if you're really that clueless, here's the quikguide:
%{
FIELD TO VARIABLE MAPPING:
trialnumtmp = TRIAL_INDEX
iptmp = IP_START_TIME
trialtmp = TRIAL_START_TIME
fixstarttmp = CURRENT_FIX_START
fixdurationtmp = CURRENT_FIX_DURATION
roitmp = CURRENT_FIX_INTEREST_AREAS
%}
function x = errorcheck(alldata, fields)
    x=0;
    %This is a stupid-ass bit of code necessary to make the roi extraction
    %work. Basically everything until the next blank line makes it so the
    %input variable to this function and the roi variable name are the same
    %matrix size so we can compare them and do all the extra ROI column
    %processing. Just don't mess with this unless necessary.
    varname = blanks(14);
    comparitor = inputname(2);
    roicomparitor = 'roitmp        ';   
    for xx = 1:length(comparitor)
       varname(1,xx) = comparitor(1,xx);
    end
    
    for i = 1:size(alldata,2);            
        %If the field exists, extract data and convert appropriately
        if fields(1,i) == 1
            x = alldata(2:end,i);
            if varname == roicomparitor %roitmp is a special case (output is a string (e.g. [4] or [1,4] instead of just 4) instead of a number)
                %like all the other variables so it gets its own processing instructions.
                %IMPORTANT: Since the IA output gives EVERY IA
                %you're in, including bg, if you're ever in a non-bg IA you get 2 IA
                %values (e.g. [1, 6] if you're in IA1 and bg is IA6). As long as you make it so the bg is your final ROI you 
                %just need to extract the first value here.
                roichar = char(x); roichar = roichar(:,3);
                %IMPORTANT: There are rare situations in which a fixation occurs
                %not anywhere (even bg) and those just come up as '[]' so you have
                %to convert them to 0 to be able to use in the analysis
                for k = 1:length(roichar)
                    if roichar(k,:) == ']'
                        roichar(k,:) = '0';
                    end
                    roiconvert(k,1) = str2num(roichar(k,:));
                end
                %Will set the 'out of bounds ROI to whatever the max value
                %was+1 so it can have its own matrix to be stored in. (e.g.
                %if you have 6 ROIs (bg = 6) the 'out of bounds' ROI will
                %be set to 7)
                new = max(roiconvert) + 1;
                for k = 1:length(x)
                    if roiconvert(k,1) == 0
                        roiconvert(k,1) = new;
                    end
                end
                x = roiconvert;
            else
                x = cell2mat(x);
            end
        end
    end    
    if x == 0   %Otherwise variable is not present; send error
        error('FATAL ERROR: Your input does not contain the NECESSARY variable: %s', char(varname));
    end
end

%{
4. This function takes the raw trial numbers of any range (e.g. 141-200 for
a multiple phase experiment) and shifts them to start at 1.
Like I said, I'm doing this a weird way but sue me, it's the first thing I thought of.
IMPORTANT NOTE: There is literally no way to program it so the script
will know if the first or last trial is bad or not unless you hard code
either the first or last trial number - which is exactly what we're
trying to avoid. Thus, the output here assumes both the first and last
trials have at least one fixation in them. Know your data: If this is not
the case, the ultimate excel output will still be correct but you will
have to numerically shift it appropriately (e.g. make trial "1" trial 2, etc.)
%}
function trialnumber = organizetrialnumbers(rawtrialnumber)
    tmptrials = rawtrialnumber;
    %Create temp storage matrix with raw trial numbers ONLY for first fix in that trial.
    %Set all other values to 0.
    for i = 2:length(rawtrialnumber)
        if rawtrialnumber(i) ~= rawtrialnumber(i-1)
            ;
        else
            tmptrials(i) = 0;           
        end
    end
    %Make matrices of just trial numbers, skipping numbers where no
    %fixations occurred. 
    newtrialarray = tmptrials(tmptrials~=0);
    newtrialindexstarts = find(tmptrials);
    trialnumber = ones(length(rawtrialnumber), 1);
    difference = ones(length(newtrialindexstarts), 1);
    %Then figure out which trials were skipped because of no fixations.
    for i = 2:length(newtrialindexstarts)
       difference(i) = difference(i-1) + (newtrialarray(i) - newtrialarray(i-1));
    end
    
    trialnumber(newtrialindexstarts) = difference; %Can't believe this worked tbh, Jesus.
    %Finally, assigns new converted trial numbers to all fixation rows
    for i = 1:length(difference)
        j = newtrialindexstarts(i);
        if i == length(difference)
            trialnumber(j:end) = difference(i);
        else
            k = newtrialindexstarts(i+1)-1;
            trialnumber(j:k) = difference(i);
        end
    end
end