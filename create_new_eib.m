function eibid = create_new_eib(ratname, eib_num, made_by, notes, region)
% function eibid = create_new_eib(ratname, eib_num, made_by, notes, region)
%  adds new rat to eibs table so that spikes can be synced in database
%  example usage add_rat_to_eibs('H037','AEH001','Ahmed','32 channels','bl FOF')

% before using, make sure that you understand that structure of the table
% you can examine the structure using bdata('explain ratinfo.eibs')
% you can gather up example entries using:
%   [eibid,ratname,eib_num,made_by,made_on,notes,region] = bdata('select * from ratinfo.eibs');

eibid = max(bdata('select eibid from ratinfo.eibs'))+1;


bdata(['insert into ratinfo.eibs ' ... 
    '(eibid, ratname, made_by, notes, region, eib_num) values ' ...
    '("{Si}", "{S}",   "{S}",   "{S}", "{S}",  "{S}")'], ...
    eibid, ratname, made_by, notes, region, eib_num);