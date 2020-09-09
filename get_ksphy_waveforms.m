function waves = get_ksphy_waveforms(sess_name)
% function get_ksphy_waveforms(sess_name)
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
%%
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
wave_x      = -6:25;
event_clus  = [];
event_ind   = [];
event_waves_uv = [];
nwaves  = 2000;
ch2tt   = @(ch, bb) ceil(ch / 4) + (bb-1)*nchperb/4;
% Loop over bundles to get cluster information and figure out which
% tetrodes we need to load to get cluster waveforms
clear sp S
%S(nbundles) = struct();
clu_inc = 0;
for bb = 1:nbundles;
    % load up spike times and cluster ids from phy
    bundle_dir      = bundle_dirfun(bb);
    clus_info_path  = fullfile(bundle_dir,'cluster_info.tsv');
    sp          = loadKSdir(bundle_dir);
    sp.mua      = sp.cgs == 1;
    sp.single   = sp.cgs == 2;

    % Find which tetrode each cluster is on using cluster info file
    fid = fopen(clus_info_path);
    C = textscan(fid, '%s%s%s%s%s%s%s%s%s%s%s%s');
    assert(strcmp(C{1}(1), 'id')); % these ids match the phy gui
    clu_info_id = cellfun(@str2num,C{1}(2:end));
    assert(strcmp(C{6}(1), 'ch')); % channel w/ strongest template
    clu_info_ch = cellfun(@str2num, C{6}(2:end));
    ncids = length(sp.cids);
    sp.ch1  = nan(size(sp.cids)); % 1 indexed channel id (note: it's 0 indexed in phy)
    sp.tt1  = nan(size(sp.cids)); % 1 indexed tetrode id
    sp.cid1 = nan(size(sp.cgs)); % 1 indexed cluster id. numbers span bundles so we don't have non-unique cluster ids
    sp.nspk = nan(size(sp.cgs)); % how many spikes are in each cluster
    
    
    % loop over good/mua clusters
    for cc = 1:ncids
        clu_inc     = clu_inc + 1;
        clu_ix      = sp.cids(cc); % get phy cluster id
        info_ix     = clu_info_id == clu_ix; % find index for this cluster in info file
        sp.ch1(cc)  = clu_info_ch(info_ix) + 1; % best channel for cluster (indexed from 1)
        sp.tt1(cc)  = ch2tt(sp.ch1(cc),bb); % convert best channel num to best tetrode num
        sp.cid1(cc) = clu_inc;
        sp.nspk(cc) = sum(sp.clu == clu_ix);     
    end
    S(bb) = sp;
end
%%
clear sp

% create a filter for the waveforms
assert(sum(diff([S.sample_rate]))==0) % check that all the files have same sampling rate
fs          = S(1).sample_rate;
[n,f0,a0,w] = firpmord([0 1000 6000 6500]/(fs/2), [0 1 0.1], [0.01 0.06 0.01]);
spk_filt    = firpm(n,f0,a0,w,{20});

cids    = [S.cids];
tt1     = [S.tt1];
ch1     = [S.ch1];
st      = vertcat(S.st);
clu     = vertcat(S.clu);
ncids   = length(cids);

    
    
event_waves = nan(ncids, nwaves, length(wave_x), 4);

% Load each tetrode and compute waveforms
active_tts = unique(tt1);
clu_inc = 0; 
for tt = 1:length(active_tts)
    this_tt     = active_tts(tt);
    active_clu  = cids(tt1 == this_tt)
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
        clu_inc     = clu_inc + 1;
        spk_ix_keep = nan(1,nwaves);
        this_clu    = active_clu(cc);
        cx          = cids == this_clu;
        this_spk_ix = round(st(clu==this_clu)*fs);
        this_nspk   = length(this_spk_ix);
        nspk(cx)    = this_nspk;
        rp_spk_ix   = this_spk_ix(randperm(this_nspk));
        keep        = 1:min([nwaves this_nspk]);
        spk_ix_keep(keep) = sort(rp_spk_ix(keep));
        % wave_inds   = spk_ix_keep(:) + repmat(wave_x,nwaves,1);
        for ss = 1:length(keep)
            tmpWf = dat_filt(:,spk_ix_keep(ss)+wave_x);
            event_waves(cx,ss,:,:) = tmpWf';
        end
        %waveFormsMean(curUnitInd,:,:) = squeeze(nanmean(waveForms(curUnitInd,:,:,:),2));
        %disp(['Completed ' int2str(curUnitInd) ' units of ' int2str(numUnits) '.'])
    end
end

waves.wv_mn    = squeeze(nanmean(event_waves,2));
waves.wv_std   = squeeze(nanstd(event_waves,[],2));
waves.wave_t_s = wave_x/fs;
waves.fs = fs;
waves.wave_x = wave_x;

waves.event_clus = [];
event_ind
event_ts
