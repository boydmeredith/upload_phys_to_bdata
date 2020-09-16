function waveS = get_ksphy_waveforms(sess_name)
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

ratname     = 'H191';
phys_dir    = fullfile(brody_dir,'/RATTER/PhysData');
raw_dir     = fullfile(phys_dir, 'Raw');
sorted_dir  = fullfile(phys_dir, 'Sorted');
sess_dir    = fullfile(sorted_dir, 'Ahmed/SpikeGadgets/', ratname, sess_name);
mda_dir     = fullfile(raw_dir, 'Ahmed/SpikeGadgets/', [sess_name '.mda']);
bndl_dirfun = @(bb) fullfile(sess_dir,sprintf('%s_bundle%i',sess_name,bb));
mda_filefun = @(tt) fullfile(mda_dir,sprintf('%s.nt%i.mda',sess_name,tt) );
clus_notes_path = fullfile(sess_dir,'cluster_notes.txt')
save_name   = fullfile(sess_dir,'waves.mat');
uv_per_bit  = 1;
warning('uv per bit conversion ratio unknown')
nchperb     = 32;
nbundles    = 4;
wave_x      = -6:25;
nwaves      = 2000;
ch2tt       = @(ch, bb) ceil(ch / 4) + (bb-1)*nchperb/4;
% Loop over bundles to get cluster information and figure out which
% tetrodes we need to load to get cluster waveforms
%S(nbundles) = struct();
sess_date   = sess_name((0:7)+regexp(sess_name,'\d\d\d\d\d\d\d'))
clus_notes_hdr = sprintf('%s\n%s\n%s',datestr(sess_date),...
        ratname, experimenter)
    "\nTT%i - "
"%i single \n "
"%i multi "

for bb = 1:nbundles;
    % load up spike times and cluster ids from phy
    bundle_dir  = bndl_dirfun(bb);
    cinf_path  = fullfile(bundle_dir,'cluster_info.tsv');
    sp          = loadKSdir(bundle_dir);
    sp.mua      = sp.cgs == 1;
    sp.single   = sp.cgs == 2;

    % Find which tetrode each cluster is on using cluster info file
    fid = fopen(cinf_path);
    C   = textscan(fid, '%s%s%s%s%s%s%s%s%s%s%s%s');
    assert(strcmp(C{1}(1), 'id')); % these ids match the phy gui
    
    cinfo_id = cellfun(@str2num,C{1}(2:end));
    assert(strcmp(C{6}(1), 'ch')); % channel w/ strongest template
    
    cinf_ch = cellfun(@str2num, C{6}(2:end));
    ncids   = length(sp.cids);
    sp.ch1  = nan(size(sp.cids)); % 1 indexed channel id (note: it's 0 indexed in phy)
    sp.tt1  = nan(size(sp.cids)); % 1 indexed tetrode id
    sp.cid1 = nan(size(sp.cgs)); % 1 indexed cluster id. numbers span bundles so we don't have non-unique cluster ids
    sp.nspk = nan(size(sp.cgs)); % how many spikes are in each cluster

    % loop over good/mua clusters
    for cc = 1:ncids
        clu_ix      = sp.cids(cc); % get phy cluster id
        info_ix     = cinfo_id == clu_ix; % find index for this cluster in info file
        sp.ch1(cc)  = cinf_ch(info_ix) + 1; % best channel for cluster (indexed from 1)
        sp.tt1(cc)  = ch2tt(sp.ch1(cc),bb); % convert best channel num to best tetrode num
        sp.nspk(cc) = sum(sp.clu == clu_ix);   
    end
    S(bb) = sp;
end
clear sp

% create a filter for the waveforms
assert(sum(diff([S.sample_rate]))==0) % check that all the files have same sampling rate
fs          = S(1).sample_rate;
[n,f0,a0,w] = firpmord([0 1000 6000 6500]/(fs/2), [0 1 0.1], [0.01 0.06 0.01]);
spk_filt    = firpm(n,f0,a0,w,{20});


% Compute each cluster's mean waveform by loading relevant tetrode,
% filtering the signal and pulling relevant indices
nactivetts  = length(unique([S.tt1]));
waveS(nactivetts) = struct();

tt_ix = 0;
for bb = 1:nbundles
    active_tts = unique(S(bb).tt1);
    fs = S(bb).sample_rate;
    % loop over trodes in this bundle with clusters
    for tt = 1:length(active_tts)
        % Figure out which clusters are on this tetrode
        tt_ix       = tt_ix + 1;
        this_tt     = active_tts(tt);
        this_mda    = mda_filefun(this_tt);
        this_tt_ind = S(bb).tt1 == this_tt;
        active_clu  = S(bb).cids(this_tt_ind);
        is_mua      = S(bb).mua(this_tt_ind);
        is_single   = S(bb).single(this_tt_ind);
        
        % Load this tetrode and filter it
        fprintf('loading mda file for tetrode %i...',this_tt);
        tic;
        dat    = readmda(this_mda);
        toc;
        fprintf('filtering');
        tic;
        dat_filt    = filtfilt(spk_filt,1,uv_per_bit*dat')';
        toc
        
        % intialize variables of interest for this tetrode
        tt_nspk  = sum(S(bb).nspk(this_tt_ind));
        ev_st    = nan(tt_nspk,1);
        ev_ind   = nan(tt_nspk,1);
        tt_clu   = nan(tt_nspk,1);

        n_clu_on_tt = length(active_clu);
        event_waves = nan(n_clu_on_tt, nwaves, length(wave_x), 4);
        spk_ix_keep = nan(n_clu_on_tt,nwaves); % indices used in mean waveform

        end_ind     = 0; % keep track of last entry into tetrode
        for cc = 1:n_clu_on_tt           
            this_cid    = active_clu(cc);
            this_st     = S(bb).st(S(bb).clu == this_cid);
            this_spk_ix = round(this_st * fs);
            this_nspk   = length(this_st);
            start_ind   = end_ind + 1;
            end_ind     = end_ind + this_nspk;
            ind         = start_ind:end_ind;
            ev_st(ind)  = this_st;
            ev_ind(ind) = this_spk_ix;
            tt_clu(ind) = ones(size(this_st)) * cc;
            keep        = 1:min([nwaves this_nspk]);
            rp_spk_ix   = this_spk_ix(randperm(this_nspk));
            spk_ix_keep(cc,keep) = sort(rp_spk_ix(keep));

            for ss = 1:length(keep)
                tmpWf = dat_filt(:,spk_ix_keep(cc,ss)+wave_x);
                event_waves(cc,ss,:,:) = tmpWf';
            end
            
            
        end

        % put this trodenum and mda filename into wave struct
        waveS(tt_ix).mua        = mua;
        waveS(tt_ix).single     = this_single;
        waveS(tt_ix).mda_fn     = this_mda;
        waveS(tt_ix).trodenum   = this_tt;
        waveS(tt_ix).event_ind  = ev_ind;
        waveS(tt_ix).event_ts   = ev_st;
        waveS(tt_ix).event_clus = tt_clu;
        waveS(tt_ix).phy_cids   = active_clu;
        waveS(tt_ix).fs         = fs;
        waveS(tt_ix).event_wave = event_waves;
        waveS(tt_ix).wave_x     = wave_x;
        waveS(tt_ix).wave_t_s   = wave_x/fs;
        waveS(tt_ix).wv_mn      = -squeeze(nanmean(event_waves,2));
        waveS(tt_ix).wv_std     = -squeeze(nanstd(event_waves,[],2));
        waveS(tt_ix).wave_ind   = spk_ix_keep;
    end
end


save(save_name,'waveS');

