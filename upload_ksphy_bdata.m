function cellids = upload_ksphy_bdata()
%%
overwrite   = 1;
bin_dir     = 'D:\Ahmed\data_sdc_20190905_170428_fromSD.mda';
res         = get_ksphy_results('bin_dir',bin_dir,'overwrite',1);
eibid       = select('')
res(1).clusnotespath = res(1).clusnotepath;
upload_phys_bdata(res);


%% put a breakpoint in upload ksphy bdata and run this
