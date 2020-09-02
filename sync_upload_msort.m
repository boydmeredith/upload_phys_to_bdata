function sync_upload_msort(ratname,date_str,which_trodes)

if nargin < 1
    ratname         = 'H153';
end
if nargin < 2
    date_str         = '2019-09-19';
end
if nargin < 3
    which_trodes    = 1:16;    
end
experimenter    = 'Tyler';

%% set up directories
brody_dir       = '/Volumes/brody';
phys_dir        = fullfile(brody_dir, 'RATTER/PhysData/Raw');
if ~exist(brody_dir)
    error('ERROR: brody dir not mounted')
end
sorted_data_dir = fullfile(brody_dir,'jtb3/projects/long_pbups/data/phys',ratname,date_str);
ratdir          = fullfile(phys_dir,experimenter,ratname);
sessdirname     = dir(fullfile(ratdir,[date_str '*']));
sessdir         = fullfile(ratdir,sessdirname.name);

clus_notes_path = fullfile(sorted_data_dir,'cluster_notes.txt');
ntt_fn_temp     = @(trodenum) fullfile(sessdir, ['TT' num2str(trodenum) '.ntt']);
ms_fn_temp      = @(trodenum) fullfile(sorted_data_dir, ...
    sprintf('%s_%s_TT%i_waves.mat',ratname,date_str,trodenum));

if ~exist(clus_notes_path,'file')
    in = input(['Warning: couldn''t find cluster notes.'...
        ' Want to create them? (y/n)'],'s');
    if in == 'y';
        annotate_session(ratname, date_str);
    else
        return
    end
end

missing_files=0;
for tt = 1:length(which_trodes)
    trodenum = which_trodes(tt);
    if ~exist(ms_fn_temp(trodenum))
        fprintf('can not find waves file for trode %i\n ', trodenum)
        missing_files = missing_files+1;
        which_trodes(tt) = nan;
    end
end
if missing_files
   x = input('should we continue (y/n)','s')
   if lower(x) ~= 'y'
       error('missing files')
   else
       which_trodes(isnan(which_trodes)) = []
   end
end
%%

[n sessid]  = bdata(['select n_done_trials, sessid from sessions where '...
    'sessiondate="{S}" and ratname="{S}"'],date_str,ratname)
if length(sessid) > 1
    fprintf('multiple behavioral sessions from this date. need to identify recording match')
    in=input('want to choose the session?(y/n)','s');
    if lower(in) == 'y'
        in = input('which session?','s');
        desired_sess = str2num(in);
        if ismember(desired_sess,sessid)
            sessid = desired_sess;
        else
            sessid = [];
        end
    else
        sessid = [];
    end
    if isempty(sessid)
        return
    end
end

[eibid, eib_num]=bdata('select eibid, eib_num from ratinfo.eibs where ratname="{S}"',ratname);
if isempty(eibid)
    error('You must have an eibid for this rat');
elseif length(eibid) > 1
    error('multiple eibdids for this rat');
end

%%
% get syncing parameters
S       = get_sessdata(sessid);
pd      = S.pd{1};
peh     = get_peh(sessid);

ev_fn   = fullfile(sessdir, 'Events.nev'); 

res     = sync_nlx_fsm_tbm(ev_fn, peh);
rt      = res.nlx_to_fsm_betas;

for trodenum = which_trodes
    waveS           = load(ms_fn_temp(trodenum));
    us_per_sample   = 1e6/waveS.fs;
    t0              = res.t0;
    
    warning(['assuming events file has the same first timepoint as spikes...'...
        'can''t think of why this wouldn''t hold, but maybe?'])
    ts_nlx_us   = (waveS.event_ind-1) * us_per_sample + t0;
    ts_fsm_s    = ts_nlx_us*rt(1) + rt(2);
    waves_uv    = waveS.event_waves_uv;

    if ~isempty(waveS.event_clus)
        upload_msort_spikes(ts_fsm_s, waves_uv, ...
            waveS.event_clus, trodenum, ...
            sessid, rt, clus_notes_path, ntt_fn_temp(trodenum), eibid)
    end
end
%%
% incorporate cutting notes from phys sess table into cells tables
update_cutting_notes(sessid);

% calculate bad ISIs and overlap
post_process_spikes(sessid);


%{ 
this is everything that happens after sync_nlx_fsm_tbm in sync_nlx_fsm
[err] =	process_spikes(sessid,pwd, r1,force,eibid);
fprintf('Finished processing spikes\n');



% calculate bad ISIs and overlap
post_process_spikes(sessid,force);

%% Process video tracker coordinates
fprintf('Processing Video\n')
[verr] = process_video(sessid, pwd, r1, forcevideo);
fprintf('Finished processing Video\n')

%% Process video tracker targets
fprintf('Processing Video\n')
%[verr] = process_raw_video(sessid, pwd, r1, forcevideo);
fprintf('Finished processing Video\n')

%% Traverse up
cd(olddir)
%}