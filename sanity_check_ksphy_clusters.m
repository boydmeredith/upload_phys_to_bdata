sess_name   = 'data_sdc_20190905_170428_fromSD';
brody_dir   = 'Y:\';
if ~exist(brody_dir),
    error(sprintf('can''t find brody directory: %s',brody_dir));
end

expmtr      = 'Tyler';
ratname     = 'H191';
phys_dir    = fullfile(brody_dir,'/RATTER/PhysData');
raw_dir     = fullfile(phys_dir, 'Raw');
sorted_dir  = fullfile(phys_dir, 'Sorted');
sess_dir    = fullfile(sorted_dir, 'Ahmed/SpikeGadgets/', ratname, sess_name);
waves_path  = fullfile(sess_dir,'waves.mat');
load(waves_path,'waveS');
%%
trodenums   = [waveS.trodenum];
for ttind = 1:length(trodenums)
    unique_clu = unique(waveS(ttind).event_clus);
    for cc = 1:length(unique_clu)
        %%
        phy_ind     = waveS(ttind).phy_cids(cc);
        clus_ind    = waveS(ttind).event_clus==cc;
        ts_fsm      = waveS(ttind).event_ts_fsm(clus_ind);
        beh_dat = load(waveS(1).sess_match.goodpath);
        sides   = beh_dat.saved.SidesSectionDyn_previous_sides == 'r';
        viol    = beh_dat.saved.PBupsDyn_violation_history;
        hits    = beh_dat.saved.PBupsDyn_hit_history;
        peh     =  [beh_dat.saved_history.ProtocolsSection_parsed_events{:}];
        cpk_ts  = cellfun(@(x) x.wait_for_cpoke1(end), {peh.states});
        rh = cellfun(@(x) ~isempty(x.right_reward), {peh.states});
        lh = cellfun(@(x) ~isempty(x.left_reward), {peh.states});
        
        this_waves = squeeze(waveS(ttind).event_wave(cc,:,:,:));
        mean_waves = nanmean(waveS(ttind).wa,4);
        
        pks        = max(this_waves,[],3);
        
        lh_cpk_ts = cpk_ts(rh);
        rh_cpk_ts = cpk_ts(lh);
        bin_size = .001;
        sd       = .15;
        xx       = ceil(5*sd/bin_size);
        krn=normpdf(-xx:xx,0,sd/bin_size);
        krn(1:xx)=0;
        [ylh x] = spike_filter(lh_cpk_ts, ts_fsm, krn, 'pre',2,'post',3);
        [yrh x] = spike_filter(rh_cpk_ts, ts_fsm, krn, 'pre',2,'post',3);
        
        label = sprintf('trode %i cluster %i phy id %i',trodenums(ttind),cc,phy_ind);
        
        figure(1); clf
        set(figure(1),'position',[ 400   525   300   300])
        plot(x,nanmean(ylh),'r')
        hold on
        plot(x,nanmean(yrh),'b')
        title({'choice psth' label })
        
        hf = figure(2); clf
        set(hf,'position',[ 10   525   363   528])
        subplot(321)
        plot(pks(:,1),pks(:,2),'.')
        xlabel('ch 1')
        ylabel('ch 2')
        title({'peak amps' label})
        subplot(322)
        plot(pks(:,1),pks(:,3),'.')
        xlabel('ch 1')
        ylabel('ch 3')
        subplot(323)
        plot(pks(:,1),pks(:,4),'.')
        xlabel('ch 1')
        ylabel('ch 4')
        subplot(324)
        plot(pks(:,2),pks(:,3),'.')
        xlabel('ch 2')
        ylabel('ch 3')
        subplot(325)
        plot(pks(:,2),pks(:,4),'.')
        xlabel('ch 2')
        ylabel('ch 4')
        subplot(326)
        plot(pks(:,3),pks(:,4),'.')
        xlabel('ch 3')
        ylabel('ch 4')
        axis(get(gcf,'children'),'square')

        % look at average wave form
        figure(3); clf
        
        set(figure(3),'position',[ 400   75   200   400])
        
        subplot(212)
        plot(-mean_waves(:,1),'color',[.8 .0 .0],'linewidth',2)
        hold on
        plot(-mean_waves(:,2),'color',[.8 .65 .65],'linewidth',1)
        plot(-mean_waves(:,3),'color',[.65 .65 .8 ],'linewidth',1)
        
        plot(-mean_waves(:,4),'color',[0 0 .8],'linewidth',2)
                subplot(211)
        plot(mean_waves(:,1),'color',[.8 .0 .0],'linewidth',2)
        hold on
        plot(mean_waves(:,2),'color',[.8 .65 .65],'linewidth',1)
        plot(mean_waves(:,3),'color',[.65 .65 .8 ],'linewidth',1)
        
        plot(mean_waves(:,4),'color',[0 0 .8],'linewidth',2)
                title({'mean waveform' label })

        figure(4); clf
        set(figure(4),'position',[ 10   175   300   300])
        dt = .5;
        histogram((diff(ts_fsm)*1e3),0:dt:1000)
        title({'ISI' label})
        xlim([0 25])


        %%
        pause
    end
end

%% look at a psth

%% look at cluster plots


