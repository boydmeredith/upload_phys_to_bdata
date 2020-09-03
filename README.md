# Code for uploading sorted units to bdata
The core functions of this code are based on `sync_nlx_fsm`, which I think was written by Jeff Erlich. BIt would look in a given folder for `cut_*.ntt` files and corresponding `.ncs` files, pull out spike times, cluster ids, waveforms, cutting notes and other good stuff and stick them into the relevant tables in bdata. This code modularizes that process separating the step of getting the data into a form we can upload and then uploading that data. 

The goal is also to provide functionality to convert `.mda` `.ncs` and `.ntt` files sorted either manually, in `kilosort` and curated with `phy`
or with `mountainsort` into a form that the `upload_spikes.m` can work with.