function get_ksphy_waveforms(sess_name)%% loop over tetrodes and load up relevant file

% For each bundle
    % 1. use sp = loadKSdir to get the spike times and cluster ids of all mua and good clusters
        % convert sp.st into an index into the mda file (might be some
        % annoying offsets here?)
    % 2. use the cluster_info.tsv file to figure out which tetrode each cluster is
    % on
    % 3. For each tetrode
        % 1. load the relevant mda file (make sure its in uv)
        % 2. grab timepoints around each of N spikes
        % 3. average each clusters spikes together
        % 4. Append the spike times, cluster labels (relabeled to be
        % continuous + separate struct of phy labels), and average
        % waveforms and sem of waveforms - i could also try to save all
        % waveforms but i think these files will get too big - can make it
        % an optional argument

% get the syncing parameters
% sync timestamps of timing variable

if nargin  < 1
    sess_name   = 'data_sdc_20190905_170428_fromSD';
end

if ~ispc 
    brody_dir   = '/Volumes/brody';
else
    brody_dir   = 'Y:\';
end

if ~exist(brody_dir),
    error(sprintf('can''t find brody directory: %s',brody_dir));
end

ratname     = 'H191';
phys_dir    = fullfile(brody_dir,'/RATTER/PhysData');
raw_dir     = fullfile(phys_dir, 'Raw');
sorted_dir  = fullfile(phys_dir, 'Sorted');
sess_dir    = fullfile(sorted_dir, 'Ahmed/SpikeGadgets/', ratname, sess_name);
mda_dir     = fullfile(raw_dir, 'Ahmed/SpikeGadgets/', [sess_name '.mda']);
bundle_dirfun   = @(bb) fullfile(sess_dir,sprintf('%s_bundle%i',sess_name,bb));
mda_filefun     = @(tt) fullfile(mda_dir,sprintf('%s.nt%i.mda',sess_name,tt) );

uv_per_bit  = 1;
warning('uv per bit set arbitrarily')
nchperb = 32;
nbundles = 4;
nbundles = 1; warning('you''ve got nbundles==1')
wave_x      = -6:25;
event_clus  = [];
event_ind   = [];
event_waves_uv = [];
nwaves   = 2000;


for bb = 1:nbundles;
    % load up spike times and cluster ids from phy
    bundle_dir      = bundle_dirfun(bb);
    clus_info_path  = fullfile(bundle_dir,'cluster_info.tsv');
    sp              = loadKSdir(bundle_dir);
    is_mua          = sp.cgs == 1;
    is_single       = sp.cgs == 2;
    fs              = sp.sample_rate;
    spike_times     = sp.st; % spike time in seconds 
    phy_clus_id     = sp.clu; % cluster ids as shown in phy
    sp.nspk         = nan(size(sp.cgs));
    % create a filter for the waveforms 
    [n,f0,a0,w] = firpmord([0 1000 6000 6500]/(sp.sample_rate/2), [0 1 0.1], [0.01 0.06 0.01]);
    spk_filt    = firpm(n,f0,a0,w,{20});

    % Find out which tetrode each cluster is on within the bundle using the
    % cluster info file
    fid = fopen(clus_info_path);
    C = textscan(fid, '%s%s%s%s%s%s%s%s%s%s%s%s');
    assert(strcmp(C{1}(1), 'id')); % these ids match the phy gui
    clu_info_id = cellfun(@str2num,C{1}(2:end));
    assert(strcmp(C{6}(1), 'ch')); % channel w/ strongest template
    clu_info_ch = cellfun(@str2num, C{6}(2:end));
    ncids = length(sp.cids);
    ch1 = nan(size(sp.cids)); % 1 indexed channel id (note: it's 0 indexed in phy)
    tt1 = nan(size(sp.cids)); % 1 indexed tetrode id
    
    for cc = 1:ncids
        clu_ix      = sp.cids(cc); % get phy cluster id
        info_ix     = find(clu_info_id == clu_ix); % find index in info file
        ch1(cc)     = clu_info_ch(info_ix) + 1; % best channel for cluster
        tt1(cc)     = ceil(ch1(cc) / 4) + (bb-1)*nchperb/4; % tetrode with cluster
    end
    
    this_event_waves = nan(nwaves, ncids, length(wave_x), 4);

    %% Load each tetrode and compute waveforms 
    active_tts = unique(tt1);
    for tt = 1:length(active_tts)
        this_tt     = active_tts(tt);
        active_clu  = sp.cids(tt1 == this_tt);
        this_mda    = mda_filefun(this_tt);
        fprintf('loading mda file for tetrode %i...',this_tt);
        tic; 
        dat    = readmda(this_mda); 
        toc; 
        fprintf('filtering');
        tic;
        dat_filt    = filtfilt(spk_filt,1,uv_per_bit*dat')';
        toc
        % loop over clusters on this tetrode and grab waveforms
        for cc = 1:length(active_clu)
            spk_ix_keep = nan(1,nwaves);
            this_clu    = active_clu(cc);
            cx = sp.cids == this_clu;
            this_spk_ix = round(sp.st(sp.clu==this_clu)*sp.sample_rate);
            this_nspk   = length(this_spk_ix);
            sp.nspk(cx) = this_nspk;
            rp_spk_ix   = this_spk_ix(randperm(this_nspk));
            ix          = 1:min([nwaves this_nspk]);
            spk_ix_keep(ix) = sort(rp_spk_ix(ix));
            % wave_inds   = spk_ix_keep(:) + repmat(wave_x,nwaves,1);
            for ss = 1:length(ix)
                tmpWf = dat_filt(:,spk_ix_keep(ss)+wave_x);
                this_event_waves(ss,cx,:,:) = tmpWf';
            end
            %waveFormsMean(curUnitInd,:,:) = squeeze(nanmean(waveForms(curUnitInd,:,:,:),2));
            %disp(['Completed ' int2str(curUnitInd) ' units of ' int2str(numUnits) '.'])
        end
    end
    sp.waves    = this_event_waves;
    sp.wv_mn    = squeeze(nanmean(this_event_waves));
    sp.wv_std   = squeeze(nanstd(this_event_waves));
    sp.wave_x   = wave_x;
    sp.mua      = is_mua;
    sp.single   = is_single;
    sp.tt1      = tt;
    sp.ch1      = ch1;
    sp.wave_t_s = wave_x/fs;

    sp.nwaves   = nwaves;
end


