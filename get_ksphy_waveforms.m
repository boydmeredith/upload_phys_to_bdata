brody_dir   = '/Volumes/brody';
if ~exist(brody_dir),
    error(sprintf('can''t find brody directory: %s',brody_dir));
end

ratname     = 'H191'
phys_dir    = fullfile(brody_dir,'/RATTER/PhysData');
raw_dir     = fullfile(phys_dir, 'Raw');
sorted_dir  = fullfile(phys_dir, 'Sorted');
sess_name   = 'data_sdc_20190905_170428_fromSD';
sess_dir    = fullfile(sorted_dir, 'Ahmed/SpikeGadgets/', ratname, sess_name);
mda_dir     = fullfile(raw_dir, 'Ahmed/SpikeGadgets/', [sess_name '.mda']);
%% loop over tetrodes and load up relevant file

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

bundle_dirfun  = @(bb) fullfile(sess_dir,sprintf('%s_bundle%i',sess_name,bb));
mda_filefun  = @(tt) fullfile(mda_dir,sprintf('%s.nt%i.mda',sess_name,tt) );
nchperb = 32;
nbundles = 4;
nbundles = 1; warning('you''ve got nbundles==1')
wave_x      = -6:25;
event_clus  = [];
event_ind   = [];
event_waves_uv = [];
nwaves   = 2000;




for bb = 1:nbundles;
    %% load up spike times and cluster ids from phy
    bundle_dir      = bundle_dirfun(bb);
    clus_info_path  = fullfile(bundle_dir,'cluster_info.tsv');

    sp          = loadKSdir(bundle_dir);
    fs          = sp.sample_rate;
    
    [n,f0,a0,w] = firpmord([0 1000 6000 6500]/(sp.sample_rate/2), [0 1 0.1], [0.01 0.06 0.01]);
    spk_filt    = firpm(n,f0,a0,w,{20});

    uv_per_bit  = 1;
    warning('uv per bit set arbitrarily')
    
    spike_times = sp.st; % what units? how to convert?
    phy_clus    = sp.clu;
    %%
    fid = fopen(clus_info_path);
    C = textscan(fid, '%s%s%s%s%s%s%s%s%s%s%s%s');
    % the first column cluster info file is cell ids - these ids correspond to
    % what you see in the phy gui
    assert(strcmp(C{1}(1), 'id'));
    clu_info_id = cellfun(@str2num,C{1}(2:end));
    % the 6th column tells you what channel has the strongest template
    assert(strcmp(C{6}(1), 'ch'));
    clu_info_ch = cellfun(@str2num, C{6}(2:end));
    
    ncids = length(sp.cids);
    ch1 = nan(size(sp.cids)); % 1 indexed channel id (note: it's 0 indexed in phy)
    tt1 = nan(size(sp.cids)); % 1 indexed tetrode id
    for cc = 1:ncids
        clu_ix      = sp.cids(cc);
        info_ix     = find(clu_info_id == clu_ix);
        % add the 1 indexed channel and tetrode ids to the cluster struct
        ch1(cc)  = clu_info_ch(info_ix) + 1; 
        tt1(cc)  = ceil(ch1(cc) / 4) + (bb-1)*nchperb;
    end
    
    % load the mda file
    active_tts = unique(tt1);
    for tt = 1:length(active_tts)
        this_tt     = active_tts(tt);
        active_clu  = sp.cids(tt1 == this_tt);
        this_mda    = mda_filefun(this_tt);
        tic; dat         = readmda(this_mda); toc % takes 10 minutes on my computer accessing data remotely
        dat_filt    = filtfilt(spk_filt,1,uv_per_bit*dat')';
        %%
        for clu = 1:length(active_clu)
            
            spk_ix_keep = nan(1,nwaves);
            this_clu    = active_clu(clu);
            this_spk_ix = round(sp.st(sp.clu==this_clu)*sp.sample_rate);
            this_nspk   = length(this_spk_ix);
            rp_spk_ix   = this_spk_ix(randperm(this_nspk));
            ix          = 1:min([nwaves this_nspk]);
            spk_ix_keep(ix) = sort(rp_spk_ix(ix));
            
            % wave_inds   = spk_ix_keep(:) + repmat(wave_x,nwaves,1);
            this_event_waves = nan(nwaves, length(wave_x), 4);
            for ss = 1:length(ix)
                tmpWf = dat_filt(:,spk_ix_keep(ss)+wave_x);
                this_event_waves(ss,:,:) = tmpWf';

            end
            wave_mean = squeeze(nanmean(this_event_waves));
            %waveFormsMean(curUnitInd,:,:) = squeeze(nanmean(waveForms(curUnitInd,:,:,:),2));
            %disp(['Completed ' int2str(curUnitInd) ' units of ' int2str(numUnits) '.'])
        end
    end
end
wave_t_s    = wave_x/fs;
%% load the cluster info 




%% plot an example
ii      = 4;
clu_ii  = sp.cids(ii)

temp_ii = unique(sp.spikeTemplates(sp.clu==clu_ii))
ci      = temp_ii + 1
cc      = 4*sp.tt1(ii)+[-3:0];
ch      = sp.ch1(ii);
temps   = (single(sp.spikeTemplates(cc)).*squeeze(sp.temps(ci,:,cc))')';
temp    = single(sp.spikeTemplates(ch)).*squeeze(sp.temps(ci,:,sp.ch1(ii)))';

%temps   = squeeze(sp.temps(ci,:,cc));
%temp    = squeeze(sp.temps(ci,:,sp.ch1(ii)));
figure(1); clf
axes 
set(gca,'ColorOrder',bone(5),'nextplot','replacechildren')

plot(-temps,'linewidth',3)
hold on
plot(-temp,'m','linewidth',3)
%%
 % this is the number shown in phy
%clu_ix = 123

temp_ix = unique(sp.spikeTemplates(sp.clu==clu_ix))
ci      = temp_ix+1;
ch_ix = clu_info_ch(info_ix)
tt_ix = ceil((ch_ix + 1)/4)
chs     = (tt_ix-1)*4+[0:3];
temps   = squeeze(sp.temps(ci,:,chs+1));
figure(1); clf
plot(-temps)
hold on
plot(-squeeze(sp.temps(ci,:,ch_ix+1)),'k')
%%


temp_ix = unique(sp.spikeTemplates);
ci      = temp_ix(3)+1;
temps   = squeeze(sp.temps(ci,:,:));
ch_ix   = clu_info_ch(ci-1)
[a b ] = min(temps,[],2);


figure(1); clf;
plot(squeeze(temps),'color',[1 1 1].*.7)
hold on
temps   = squeeze(sp.temps(ci,:,13:16));
plot(squeeze(temps))



%%
mdapath = fullfile(bundle_dir,[sess_name '_bundle1.mda'])
recfull = readmda(mdapath);
%%
p.dataType  = 'int16'; % ?
p.nCh       = 32;
p.wfWin     = -6:25;
p.nwf       = 2000; % number of waveforms to average
p.spikeTimes    = [];
p.spikeClusters = [];


% sp.st are spike times in seconds
% sp.clu are cluster identities
gwfparams.spikeTimes = ceil(sp.st(sp.clu==155)*30000); % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = sp.clu(sp.clu==155);

if 0
    save_ms_waves_mda(mdapath, clustspath,save_path,ratname,datestr,trodenum)
end
