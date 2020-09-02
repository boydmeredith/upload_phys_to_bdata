function err = upload_spikes(sp)
% function err=process_spikes(sp)

% Takes in a struct sp describing curated spikes from a recording
% it should have fields:
%   eibid       key into eibs table for this recording device
%   sessid      key into sessions table for this recording session
%   ts_fsm_s    nspikes x 1  spike times in seconds synced to the behavioral (FSM time)
%   clust_id    nspikes x 1  cluster id associated with each spike 
%   trodenum    which tetrode was it recorded on (indexed from 1) 
%               trodenum-1 will be inserted into the sc_num field in the
%               cells table
%   rt          1x2     regression parameters for syncing the
%                       phys timestamps into fsm time
%   notes_path  path to a file containing notes on each cluster -
%               keywords 'good'/'single', 'multi', 'bad'
%   waveS       struct containing waveform mean, std, clusterid for each
%               clust_id
%   rec_path    path to file used to create spikes - probably ntt, ncs or mda

% unpack the variables and make sure we have everything we need
eibid       = sp.eibid;
sessid      = sp.sessid;
ts_fsm_s    = sp.ts_fsm_s;
clust_id    = sp.clust_id;
trodenum    = sp.trodenum;
rt          = sp.rt;
notes_path  = sp.notes_path;
waves_mean  = sp.waveS.mn;
waves_std   = sp.waveS.std;
waves_id    = sp.waveS.id;
rec_path    = sp.rec_path;

sc_num  = trodenum - 1;
assert(all(ismember(waves_id,clust_id)));

% get the ratname corresponding to the sessid
ratname = bdata('select ratname from sessions where sessid="{S}"',sessid);
ratname = ratname{1};

[fldr, ntt_fn, ext]   = fileparts(rec_path);
ntt_fn  = [ntt_fn ext];
fid     = fopen(fullfile(rec_path),'r');
hd      = fread(fid,16384);  % The header of all nlx files is 16384 bytes
hd      = char(hd');   % reformat the header to a human readable format.
    
% [ts, cell_n,  waves, param, sc_n ,hd]=nlx_spike2mat(ntt_fn);
% convert the waves into microvolts
[wsc, chans, cheetah_ver] = extract_header(hd);
%wsc = repmat(wsc(:),1,32);
%warning('multiplying waves by wsc')
%waves = waves.*wsc;

%% checkif these have been processed
chans_num = str2num(chans);


%assert(sc_num == (trodenum -1))
%chans  = sc_num*4 + (0:3);

already_done=bdata('select count(*) from cells where sessid="{Si}" and sc_num="{Si}"', sessid, sc_num);
if already_done>0
    warning('already uploaded cells for this session and channel')
    keyboard
end

%% Try to load cutting notes, and syncing parameters in the phys_sess table
in_phys_sess=bdata('select count(sessid) from phys_sess where sessid="{Si}"',sessid);
if in_phys_sess==0
    bdata(['insert into phys_sess (sessid,ratname, cutting_notes, '...
    'sync_fit_m,sync_fit_b) values ("{Si}","{S}", "{F}","{S12}","{S12}")'],...
    sessid,ratname,notes_path,rt(1),rt(2))
end


try
    % update channels table for this tetrode
    in_channels = bdata('select count(sessid) from bdata.channels where sessid="{Si}" and file_name="{S}"',sessid,ntt_fn);
    if in_channels>0
        keyboard
    else
        sqlstr=['insert into bdata.channels (sessid, ad_channels, header, file_name, path_name) ' ...
            ' values ("{Si}","{S}","{S}","{S}","{S}")'];
        bdata(sqlstr, sessid, chans, hd, ntt_fn, fldr);
        channelid=bdata('select last_insert_id()');
    end
    %%  For each cell in the spike file
    
    unique_clusters=unique(clust_id);
    for cc = 1:numel(unique_clusters)
        cl_n = unique_clusters(cc);

        % first the cell table
        nSpikes     = sum(clust_id==cl_n);
        cluster_in_file = cl_n;
        cluster = cc; % need to sc_num+1 because the channels are 0 indexed.
        sqlstr  = ['insert into bdata.cells (ratname, sessid, channelid,' ...
            'sc_num, cluster, nSpikes, filename, cluster_in_file, eibid)' ...
            'values ("{S}","{Si}","{Si}","{Si}","{Si}","{Si}","{S}","{Si}","{Si}")'];
        bdata(sqlstr,ratname,sessid, channelid, sc_num, cluster, nSpikes, ntt_fn, cl_n, eibid);
        cellid = bdata('select last_insert_id()');
        % then the spike
        this_clust = clust_id == cl_n;
        c_ts   = ts_fsm_s(this_clust);
        % squeeze breaks for single electrodes
        % 				w.mn=squeeze(mean(waves(cell_n==cl_n,:,:)));
        % permute like this removes a leading singleton dimension
        w.mn    = squeeze(mean(waves(this_clust,:,:)))';
 
        try
            % squeeze breaks for single electrodes
            %     				w.std=squeeze(std(waves(cell_n==cl_n,:,:)));
            % permute like this removes a leading singleton dimension
            w.std = squeeze(std(waves(this_clust,:,:)))';
        catch
            try
                sprintf(['not enough memory to compute std of waveform ...'
                    'recomputing sample std']);
                w.std = std(waves(this_clust.*(rand(size(this_clust))<.5),:,:))';
            catch
                sprintf('failed')
            end
        end;
        % insert cluster into spiketimes
        sqlstr = ['insert into bdata.spktimes (cellid, sessid, ts, wave) '...
            'values ("{Si}","{Si}","{M}","{M}")'];
        bdata(sqlstr, cellid, sessid, c_ts, w);
        
        
    end
catch
    showerror(lasterror);
end


err = 0;


%% get the parameters to convert the waveforms to uV

function [wscale,c_nums,cheetahv]=extract_header(hd)

n=regexp(hd,'InputRange','end');
eol=find(hd(n:end)==13,1,'first');
wscalen=str2num(hd(n+1:n+eol-1)); % THIS IS A TOTAL HACK
n=regexp(hd,'ADMaxValue','end');
eol=find(hd(n:end)==13,1,'first');
wscaled=str2num(hd(n+1:n+eol-1)); % THIS IS A TOTAL HACK
wscale=wscalen/wscaled;
n=regexp(hd,'-ADChannel','end');
c_nums=strtok(hd(n+1:n+12),'-');


n=regexp(hd,'CheetahRev','end');
eol=find(hd(n:end)==13,1,'first');
crs=hd(n+1:n+eol-1);
lastp=find(crs=='.',1,'last');
cheetahv=str2num(crs(1:lastp-1)); % THIS IS A TOTAL HACK






