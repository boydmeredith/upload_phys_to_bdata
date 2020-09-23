function res = find_wireless_sess(sess, varargin)
% function res = find_ttl_match(sess, ratlist, behav_dir, mda_dir)  
% given a set of ttls and a date, search through a list of rats to find the best physiology session to match these ttls
p = inputParser()
addParameter(p,'overwrite',0)
addParameter(p,'ratlist',{'H191','H176'})
addParameter(p,'expmtr','Ahmed')
addParameter(p,'behav_dir','')
addParameter(p,'mda_dir','')
addParameter(p,'fs',30000)
parse(p,varargin{:})

brody_dir = 'Y:';
overwrite   = p.Results.overwrite;
ratlist     = p.Results.ratlist;
expmtr      = p.Results.expmtr;
behav_dir   = p.Results.behav_dir;
mda_dir     = p.Results.mda_dir;
fs          = p.Results.fs;

if isempty(behav_dir)
    behav_dir = fullfile(brody_dir, 'RATTER/SoloData/Data', expmtr);
end
if  isempty(mda_dir)
    phys_dir    = fullfile(brody_dir, 'RATTER/PhysData/Raw/Ahmed/SpikeGadgets/');
    mda_dir     = [phys_dir sess '.mda/'];
end

save_path   = fullfile(mda_dir, 'ttl_match.mat')
if exist(save_path) & ~overwrite
    dio = load(save_path,'res','sess');
    if strcmp(dio.sess, sess)
        res = dio.res;
        return;
    else 
        clear dio
    end
end

dio_file    = [mda_dir sess '.dio_RFINPUT.dat'];
if ~exist(dio_file)
    dio_file = [strrep(mda_dir,'.mda','.DIO') sess '.dio_RFINPUT.dat'];

end
fi=dir(dio_file);
if(isempty(fi))
  error(sprintf('dio missing %s',dio_file))
end

if strcmp(sess(1:4),'data')
    ratname = '';
    date_str = sess(12:17);
    fprintf('no ratname, will look for best match on date %s',date_str)
    
else
    ratname = sess(1:4);
    ratlist = {ratname};
    date_str = sess([14 15 6 7 9 10]);
    fprintf('ratname %s, date %s',ratname,date_str)
end
sessiondate = ['20' date_str([1 2]) '-' date_str([3 4]) ...
    '-' date_str([5 6])];

% load ttls from wireless 
dio     = readTrodesExtractedDataFile(dio_file);
c       = dio.fields(1).data;
b       = dio.fields(2).data;
ttls    = double(c(b==1))/fs;

inc     = 0;
for rr = 1:length(ratlist)
    fntemp      = fullfile(behav_dir, ratlist{rr}, ['data_*' date_str '*.mat']);
    ratfiles    = dir(fntemp);
    if isempty(ratfiles)
        fprintf(['couldn''t find a match for ' fntemp])
        continue
    end
    for ff=1:length(ratfiles)
        inc             = inc +1;
        rats{inc}       = ratlist{rr};
        this_behav_name = ratfiles(ff).name;
        this_behav_path = fullfile(ratfiles(ff).folder, this_behav_name);
        filename{inc}   = this_behav_name;
        behav_path{inc} = this_behav_path;
        load(this_behav_path, 'saved_history');
        
        parsed_events   = saved_history.ProtocolsSection_parsed_events;
        % load behavior TTLS for alignment with ephys
        ttls_fsm = nan(1,length(parsed_events));
        for i=1:length(parsed_events)
            if(isfield(parsed_events{i}.states,'sending_trialnum') && ...
                    length(parsed_events{i}.states.sending_trialnum)>0 )
                ttls_fsm(i) = parsed_events{i}.states.sending_trialnum(1);
            else
                ttls_fsm(i) = NaN;
            end
        end
        ttls_fsm=ttls_fsm';
        
        % get interpulse interval of ttls
        trode_gaps  = diff(ttls);
        fsm_gaps    = diff(ttls_fsm);
        
        % look only at ttls delivered 5-10 seconds apart
        trode_gaps_valid    = find(trode_gaps > 5 & trode_gaps < 10);
        % look for correspondance between the phys and behav sess ITIs        
        ttls_fsmall=[];
        ttl_trodesall=[];
        for i=1:length(trode_gaps_valid)
            fsm_trodes_diffs=find(abs(fsm_gaps - trode_gaps(trode_gaps_valid(i))) < 0.05);
            for k=1:length(fsm_trodes_diffs)
                indb=fsm_trodes_diffs(k);
                if (indb+5)>length(fsm_gaps)
                    continue
                end
                if (trode_gaps_valid(i)+5)>length(trode_gaps)
                    continue
                end
                vec1    = trode_gaps(trode_gaps_valid(i)+(0:5));
                vec2    = fsm_gaps(indb+(0:5));
                
                %if one of the distances is too long (timeout), abort
                if(max(abs(diff(vec1)))>30)
                    continue
                end
                
                if(norm(vec1-vec2)>.1)
                    continue
                else
                    ttls(trode_gaps_valid(i));
                    ttl_trodesall=[ttl_trodesall; ttls(trode_gaps_valid(i)+(0:5))];
                    ttls_fsmall=[ttls_fsmall; ttls_fsm(indb+(0:5))];
                end
            end
        end
        
        rt{inc}     = polyfit(ttl_trodesall,ttls_fsmall,1);
        
        ttl_trode_fsm = ttl_trodesall*rt{inc}(1)+rt{inc}(2);
        
        if(isempty(ttl_trode_fsm))
            max_residual(inc)   = 100;
            totaldur(inc)       = 0;
        else
            max_residual(inc)   = max(abs(ttl_trode_fsm-ttls_fsmall));
            totaldur(inc)       = (max(ttls_fsmall)-min(ttls_fsmall))/60;
        end
    end
end
    
% figure it out if any of the behav sessions are acceptable matches
match_ind  = find(max_residual<0.02 & totaldur>5);

if(length(match_ind)<1)
    warning('no match!');
    save_path   = fullfile(mda_dir, 'no_ttl_match.mat');
elseif(length(match_ind)>1)
    disp(valore(match_ind))
    disp(filename(match_ind))
    warning('too many matches!')
    save_path   = fullfile(mda_dir, 'multiple_ttl_matches.mat');
else
    res.goodval     = max_residual(match_ind);
    res.gooddur     = totaldur(match_ind);
    res.spk2fsm_rt  = rt{match_ind};
    res.spk2fsm_fn  = @(spk_ts) spk_ts.*res.spk2fsm_rt(1) + res.spk2fsm_rt(2);
    res.goodpath    = behav_path{match_ind};
    res.behfile     = filename{match_ind};
    res.ratname     = rats{match_ind};
    if bdata('connect') > 0
        sessid          =  bdata(['select sessid from sessions where data_file="'  res.behfile(1:end-4) '"']);
        res.sessid      = sessid;
    else
    res.sessiondate = sessiondate;
end

    res.date_str    = date_str;
    res.save_path   = save_path;
    save(save_path, 'res', 'sess')
