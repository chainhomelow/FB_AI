function pca_nuisance_regressors_final(data)

%File IO
input_dir = '/media/Lapis/DougDiss_MVCS/RESTING';
acq_dir = '/media/Lapis/DougDiss_MVCS/ACQ';
ext_dir = '/media/Lapis/DougDiss_MVCS/EXT';
sub = dir(input_dir);      %Should be same in RESTING as in ACQ and EXT
sub = sub(3:end); 
subs = {sub.name};
num_subs = length(subs);

%This will put names of all nuisance BOLD files for ALL scans (INCLUDING ACQ AND EXT SCANS) - 7 per subject - into one giant
%cell vector to be loaded in next loop. The files will be organized by row.
%Each row is same sub. Each col is 1-5 a RS scan, 6 - their ACQ scan, and 7
%- their EXT scan
for subjects = 1:length(subs)
    %Remember, all files are in /RESTING/../../ALIGN_ALL now
    path = char(strcat(input_dir, '/', subs(subjects), '/', subs(subjects), '.REST.results/ALIGN_ALL'));
    r1_filenames = dir(fullfile(path, '*r01.nuisance.BOLD.zscore.nii'));
    r2_filenames = dir(fullfile(path, '*r02.nuisance.BOLD.zscore.nii'));
    r3_filenames = dir(fullfile(path, '*r03.nuisance.BOLD.zscore.nii'));
    r4_filenames = dir(fullfile(path, '*r04.nuisance.BOLD.zscore.nii'));
    r5_filenames = dir(fullfile(path, '*r05.nuisance.BOLD.zscore.nii'));
    acq_filenames = dir(fullfile(path, '*ACQ.nuisance.BOLD.zscore.nii'));
    ext_filenames = dir(fullfile(path, '*EXT.nuisance.BOLD.zscore.nii'));
    num_runs = 7; %Fuck it
    
    if subjects == 1
        files = cell(num_subs, num_runs);
            files{subjects,1} = fullfile(path, r1_filenames.name);
            files{subjects,2} = fullfile(path, r2_filenames.name);
            files{subjects,3} = fullfile(path, r3_filenames.name);
            files{subjects,4} = fullfile(path, r4_filenames.name);
            files{subjects,5} = fullfile(path, r5_filenames.name);
            files{subjects,6} = fullfile(path, acq_filenames.name);
            files{subjects,7} = fullfile(path, ext_filenames.name);
    else
            files{subjects,1} = fullfile(path, r1_filenames.name);
            files{subjects,2} = fullfile(path, r2_filenames.name);
            files{subjects,3} = fullfile(path, r3_filenames.name);
            files{subjects,4} = fullfile(path, r4_filenames.name);
            files{subjects,5} = fullfile(path, r5_filenames.name);
            files{subjects,6} = fullfile(path, acq_filenames.name);
            files{subjects,7} = fullfile(path, ext_filenames.name);
    end
end

%This loads each file one at a time from the cell vector 'files' which just
%contains all 1x7 (or whatever) nuisance BOLD timeseries and runs the
%analysis code on them
min_confirmed_eigs = zeros(length(files),1);
%Variable 'confirmed_eigs' is number of retained Eigenvalues after Monte
%Carlo sim
confirmed_eigs = zeros(num_runs,1);

%First find number of Eigenvalues to keep
for subjects = 1:length(files)
    for scans = 1:size(files,2)
        im = load_nii(files{subjects,scans});
        imag = im.img;
        imagdouble = double(imag);

        %This will output a voxelXtime (i.e. voxelX236) matrix which we MUST transpose for computation purposes
        %See here: http://stats.stackexchange.com/questions/197071/removing-nuisance-pca-components-from-the-fmri-data/197141#197141
        %WHOOPS ONLY RS SCANS ARE 236 - ACQ & EXT ARE 286.
        %%Masking
        temp = imagdouble(:,:,:,1); %Just need spatial info
        all_ind = find(or(temp == 0, temp > 0));    %For reconstruction
        ind_numbers = find(temp);                   %Actually gets numerical indices
        ind = (temp ~= 0); 
        length_ind = length(find(ind ~= 0));
        masked_scan = zeros(length_ind, size(imagdouble,4));
        %loop over files
        for TR = 1:size(imagdouble,4)
            temp_scan = squeeze(imagdouble(:,:,:,TR));
            temp_scan(isnan(temp_scan)) = 0;
            temp_scan = temp_scan(ind);
            masked_scan(:, TR) = temp_scan(:);
        end
        
        masked_scan = masked_scan';
        eigens = size(masked_scan,1);
        [COEFF, SCORE, LATENT] = pca(masked_scan, 'Centered', false);       %DO NOT CENTER: Already Zscored
        
        %%RUN MONTE CARLO SIMS TO FIND NUMBER OF EIGENVECTOR TO RETAIN
        nShuffle = 1000;        %MC iterations
        alpha = 0.01;           %99th CI
        xShuffle = masked_scan;           %Orig data
        latentShuffle = zeros(length(LATENT), nShuffle);    %e.g. 236x1000
        for iShuffle = 1:nShuffle
            for dim = 1:size(masked_scan,2)
                    xShuffle(:, dim) = masked_scan(randperm(size(masked_scan,1)), dim);
            end
            [~, ~ ,latentShuffle(:,iShuffle)] = pca(xShuffle);
        end
        latentHigh = quantile(latentShuffle', 1-alpha)';
        latentLow = quantile(latentShuffle', alpha)';
        
        %%Now compare nuisance PC Eigenvalue (LATENT) to 99% CI of norm dist
        %%Eigenvalue
        keepers = zeros(size(LATENT,1), 1);

        for l = 1:eigens
            if LATENT(l, 1) > latentHigh(l,1)
                keepers(l, 1) = LATENT(l,1);
            else
                keepers(l, 1) = 0;
            end
        end

        %Remove 0s
        keepers(keepers == 0) = [];
        num_good_eigs = length(keepers)
        confirmed_eigs(scans, 1) = num_good_eigs;
        
        %%%Stick all SCORES into one 3D matrix so as to not recompute
        if scans == 1
           all_SCORES = zeros(size(SCORE, 1), size(SCORE, 2), size(files,2));   %e.g. 240xPCsx7
           all_SCORES(:,:,scans) = SCORE;
        else
            all_SCORES(:,:,scans) = SCORE;
        end
    end
    min_confirmed_eigs(subjects) = min(confirmed_eigs);
    
    %%Now actually take out the PCs by reloading the variables and taking
    %%out the minimum required number
    for scans = 1:size(files,2)   
        %Get out heade info for saving as NII later
        voxels = im.hdr.dime.pixdim(2:4);
        origin = im.hdr.hist.originator(1:3);       %I honestly don't know if it's 1:3 or 2:4 so make sure to visually check
    
        
    %%%%%%%%%%%%%   SEE IF YOU CAN RECONSTRUCT JUST FOR FUN     %%%%%%%%%%%
    [~,reconstructed] = pcares(masked_image, min_confirmed_eigs(subject));
    rebuild = zeros(size(imagdouble,1), size(imagdouble,2), size(imagdouble,3), size(imagdouble,4));
    for n = 1:size(imagdouble,4)
        for o = 1:size(ind_numbers);
            rebuild(ind_numbers) = reconstructed(n,o);
        end
    end
    niiFILE = make_nii(rebuild, voxels, origin);
    save_nii(niiFILE, 'test_writeout');
    
        %%%Get out the Principal Components from SCORE
        PCs = squeeze(all_SCORES(:, 1:min_confirmed_eigs, scans));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%           I/O             %%%%%%%%%%%%%%%%%%
       if scans == 6
           dlmwrite(fullfile(path,'ACQ.nuisance.PCs.txt'), PCs, 'delimiter','\t'); 
       elseif scans == 7
           dlmwrite(fullfile(path,'EXT.nuisance.PCs.txt'), PCs, 'delimiter','\t');
       else
           dlmwrite(fullfile(path,'ACQ.nuisance.PCs.r0', scans, '.txt'), PCs, 'delimiter','\t');
       end
    end
       %%%%%%AND HERE'S THE KICKER - your ultimate output will be SIGNIFICANTLY SMALLER than your original
       %%%%%%NII - like 236x8 or something. I'm pretty sure you just output
       %%%%%%it as the equivalent of a 1D file or some shit and use it
       %%%%%%similar to a motion parameter file into 3dREMLfit or
       %%%%%%3dTproject or some shit
        
       %%E.G. Next command: 1d_tool.py -infile pc.txt -write pc1.1D
       %%Then: 3dTproject -input 1R103T1A3d.nii -ort pc1.1D -prefix new_out
end
end
        