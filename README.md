# Code for uploading sorted units to bdata
The core functions of this code are based on `sync_nlx_fsm`, which I think was written by Jeff Erlich. 
It would look in a given folder for `cut_*.ntt` files and corresponding `.ncs` files, 
pull out spike times, cluster ids, waveforms, cutting notes and other good stuff and stick them into the 
relevant tables in bdata. It assumed you were using a unilateral neuralynx array.

This code modularizes by creating a common format that we put the data in before updating bdata.

## Okay, how do I use it
In order to use this code, you need to create a `spkS` structure and then call 
`upload_phys_bdata(spkS)`. 

`upload_phys_bdata` Expects `spkS` to contain the following fields:
-   `ratname`         str 
-   `sessid`          int     key to get relevant behavior session from sessions table
-   `eibid`           int     key into ratinfo.eibids for this recording device
-   `trodenum`        int     tetrode number on this eib indexed from 1
-   `recpath`         str     path to the recording for this tetrode
-   `fs`              int     sampling rate of recording
-   `hd`              str     header of recording - 
-   `clusnotespath`   str     path to a txt file containing notes on each cluster
-   `sync_fit_m`      1x1     regression slope for converting phys ttls to FSM time  
-   `sync_fit_b`      1x1     regression intercept for converting phys ttls to FSM time  
-   `event_ts_fsm`    nspikes x 1     spike times in FSM time
-   `event_clus`      nspikes x 1     cluster id for each spike 
-   `waves_mn`        nclusts x nchpertrode x ntimepts    mean waveform in uv
-   `waves_std`       nclusts x nchpertrode x ntimepts    std of waveform
-   `waves_ind`       nwaves  x 1     indices of waves used to compute mean
-   `waves_clus`      nclusts x 1     matches waves_mn to event_clus and to cluster
                                    notes
## How do I get a `spkS` structure
You can make `spkS` based on any kind of data! And, if you want to do it for a spike gadgets
recording with cluster cutting performed in kilosort and curated in phy, we can even help. 
There's also a version for neuralynx recordings and mountainsort spike sorting, but you have to ask 
Tyler about it.

### SpikeGadgets -> Kilosort2 -> Phy
To get a `spkS` struct for this pipeline, run 
`spkS = get_ksphy_results('bin_dir',<path/to/your/binary_parent_directory>)`. 
At the moment, you will have to do some work to sort the paths out, but that should be fixed soon.
