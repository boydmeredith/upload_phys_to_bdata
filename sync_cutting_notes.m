function sync_cutting_notes(sessid)
% this function reads the cutting notes already uploaded to the phys_sess table and
% uses it to update the cutting comments and single boolean in the cells
% table
% TO DO: check that cluster and cluster in file are assigned properly

[cn]=bdata('select cutting_notes from phys_sess where sessid="{S}"',sessid);


fprintf(['\nIncorporating the following cutting notes: \n\n' char(cn{1}') '\n'])

S=parse_cutting_notes(cn{1});
if isempty(S)
    fprintf(1,'No cutting notes for session %d\n',sessid);

else
    for sx=1:numel(S)
        cellid=bdata(['select cellid from cells  where sessid="{S}" and '...
            'sc_num="{S}" and cluster="{S}"'],sessid, S(sx).TT-1, S(sx).SC);
        bdata('call update_cutting("{S}","{S}","{S}")',cellid, ...
            S(sx).single, S(sx).cutting_comment);
    end
    fprintf(1,'Done session %d. notes for %d cells addded.\n',sessid,numel(S));
end