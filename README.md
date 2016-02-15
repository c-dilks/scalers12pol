RELATIVE LUMINOSITY ANALYSIS
----------------------------

Installation
------------
  1. run `install` to make directories and build symlinks; you should
     have this current directory as a subdirectory of root12fms
  2. build scaler bit reader by typing `make clean` and then `make`
  3. download scaler data, and drop it in the directory `./sca2013`
     (see instructions under `detailed version`)
  4. download spin pattern data, and drop it in `./spinpat`



Analysis Procedure
------------------

1. Download scaler files from HPSS
  - 2013 instructions
    - code contained in `hpss` subdirectory locally; on rcas found in `~/scalers2013`
      - `BuildList` uses `goodruns.dat` (see trgmon code) to build a list of
        scaler files to retrieve from HPSS; the variable `size` controls
        how many scaler files will be listed in each HPSS request list
        - note that a board number is needed; see usage note
        - the target is set to `$HOME/scratch/sca$YEAR/...`
      - run `hpss_user.pl -f [request list]` to add the list of files
        and respective targets to the HPSS request queue and wait
      - monitor status with `hpss_user.pl -w`
  - 2012 instructions
    - scaler files on HPSS for run12 were saved including a UNIX timestamp, which
      makes it difficult to generate a list of files
    - first, type `hsi` and cd to `/home/starsink/raw/daq/2012`
    - then type `out > FILE_LIST`; this will create a file called `FILE_LIST` in 
      your current local (RCAS) directory and pipe all output from hsi to this file
    - run12 rellum uses board 12 data, so to list all the board 12 scaler files, type
      `ls */*/*_12_*.sca` and be patient
    - when the command is done, exit HPSS and check `FILE_LIST` for the proper output
    - execute `BuildList_2012` to build lists of files to submit to the data carousel; 
      `files_to_retrieve*.lst` will be created, with 150 requests per file
      - submit using `hpss_user.pl -f [list]`
      - check carousel status using `hpss_user.pl -w`
  - 2011 instructions
    - need runlist; I didn't do run QA for run11, so I wrote a couple short scripts
      to extract the run number and fill numbers present in the available output
      files; run `../../get_run_list.C` followed by `../../append_fill_numbers`;
      the file `runlist_with_fills` can be copied to `scalers11t/goodruns.dat`
      (and to wherever you need it in order to download the scaler files)
    - scaler files on HPSS for run11 were saved including a UNIX timestamp, which
      makes it difficult to generate a list of files
    - first, type `hsi` and cd to `/home/starsink/raw/daq/2012`
    - then type `out > FILE_LIST`; this will create a file called `FILE_LIST` in 
      your current local (RCAS) directory and pipe all output from hsi to this file
    - run11 rellum uses board 3&4 data, so, for example,  to list all the board 3 scaler files,
      type `ls */*/*_3_*.sca` and be patient 
      - board 3 is read out once every run
      - board 4 is read out once every 1000 seconds and once at the end of the run
    - when the command is done, exit `HPSS` and check `FILE_LIST` for the proper output
    - execute `BuildList_2011` to build lists of files to submit to the data carousel; 
      `files_to_retrieve*.lst` will be created, with 150 requests per file
      - submit using `hpss_user.pl -f [list]`
      - check carousel status using `hpss_user.pl -w`

2. Read the scalers scaler bit reader reader from Zilong (32-bit for run13 and 24-bit for run12)
  - 32-bit scaler reader (for 2013)
    (obtained from `~zchang/2013-03-scaler/codes_scaler/MyCodes/scaler2_reader_bit.c`)
    - execute `read_scalers`
      - reads all scaler files in `/GreyDisk1/sca2013`
      - executes `scaler2_reader_bit.exe*` via condor
      - outputs the read files in `datfiles` directory
      - explanation of `scaler2_reader_bit.exe` (see make file for compilation)
        - reads 32-bit scaler files obtained from hpss (downloaded to `/GreyDisk1/sca2013`)
        - outputs corresponding file in `datfiles` directory with the following columns:
          - bunch crossing (bXing) number
          - BBC(0-7) (see zilong's scaler bit definitions below)
          - ZDC(0-7) 
          - VPD(0-3)
          - total bXings
  - 24-bit scaler reader (for 2012)
    - `read_scalers` in `scalers12` directory executes `sca_read_bin.o` via condor
      - reads all scaler files in `/GreyDisk1/sca2012`
      - outputs the read files in `datfiles` directory
      - `datfile` columns
        - bunch crossing (bXing) number
        - BBC(0-7) (see zilong's scaler bit definitions below)
        - ZDC(0-7) 
        - VPD(0-7)
        - total bXings
  - 24-bit scaler reader (for 2011)
    - since we have the choice of using board 3 or 4, we select one
      by passing a board number to `read_scalers`
    - before running this, make two directories: `datfiles_bd3` and `datfiles_bd4`;
      make `datfiles` a symlink to the directory that corresponds to the board
      you're analyzing; this symlink should remain unchanged throughout the rest of the
      analysis (in other words, to change board number, you must redo the analysis
      starting from the datfiles)
      - board 3 reads out once per run
      - board 4 reads out once every 1000 seconds and at the end of the run
    - `read_scalers` in `scalers11t` directory executes `sca_read_bin.o` via condor
      - reads all scaler files in `/GreyDisk1/sca2011t`
      - outputs the read files in `datfiles` directory
      - datfile columns
        - bunch crossing (bXing) number
        - BBC[0-7] (see zilong's scaler bit definitions below)
        - ZDC[0-7] 
        - VPD[0-7]
        - total bXings
    - see `bit_doc.txt` for further information
    - if there are multiple scaler files for a run (board 4 seems to read out 
      every 1,000 seconds), you can run `datadd` to add the columns of each run;
      the original, un-added datfiles are backed up into `datfiles/orig`

3. Obtain spin patterns
  - Verify we have the fill numbers (`fill.txt`) and the corresponding
    spin patterns for fills listed in `goodruns.dat`; if you are unsure: 
    - run `getspinpat` to recreate `fill.txt` from `goodruns.dat`
      and download spin patterns from CDEV to `./spinpat`
    - append `$FMSTXT/fill.txt` with `./fill.txt` (then cat through uniq to 
      remove duplicated lines.. not necessary, but it's nice to do)
      - `cat fill.txt >> $FMSTXT/fill.txt`
      - `cat $FMSTXT/fill.txt | uniq > $FMSTXT/fill.tmp`
      - `mv $FMSTXT/fill.{tmp,txt}`
    - copy downloaded spinpats to `$FMSTXT/spinpat/`
  - SEE SPECIAL NOTE BELOW REGARDING F17600
  - run `spin_cogging` to create `spinpat/[fill].spin` files from the downloaded CDEV files
    - since STAR spin is opposite source spin and source spin is same as CDEV spin,
      this script implements a sign flip
  - SPECIAL NOTE F17600: I ran a modifier script called `mod_17600`
    to fix spin pattern for this fill.. see the comments in the
    script for further details; my cdev file is the modified 
    spin pattern for the last 14 runs (`*.bad` marks the unmodified 
    pattern from CDEV, where the runs should be omitted from analyses)

4. Accumulate the scaler data into one table: `run accumulate`
  - bunch kicking: you must first generate a list of kicked bunches using `bunch_kicker`,
    if one does not exist; the list is in the text file `kicked`, with columns
    fill, bx, spinbit
    - use the variable `run_randomizer`; if it's zero, spinbit equalisation does not run
      (see spinbit equalization section below)
    - bunches which are manually removed are listed in the beginning of the script
    - explanation of algorithm is in the comments
  - execute `accumulate`
    - collects all the datfiles into `datfiles/acc.dat`, with columns:
      (** and filters out bad part of 17600)
      - run index
      - runnumber
      - fill
      - run time (seconds)
      - bunch crossing (bx)
      - BBC[1-8]
      - ZDC[1-8]
      - VPD[1-4]
      - total bXings
      - blue spin
      - yellow spin
  - `accumulate` then creates `counts.root,` which contains useful trees, using `mk_tree.C`
  - `mk_tree.C` reads `datfiles/acc.dat` file
    - BXING SHIFT CORRECTIONS ARE IMPLEMENTED HERE (FOR RUN 12 ONLY!!!)
    - acc tree: simply the `acc.dat` table converted into a tree
    - sca tree: restructured tree containing containing branches like 
      - bbc east, bbc west, bbc coincidence 
      - similar entries for zdc and vpd
      - run number, fill number, bunch crossing, spin bit
      - num_runs = number of runs in a fill
      - kicked bunches: 
        certain bunches which are empty according to scalers but filled according to
        cdev are labelled as 'kicked' in the output tree; bXings which are kicked
        are not projected to any distributions in the rellum4 output

5. Compute the relative luminosity
  - `rellum4.C` is the analysis script that reads `counts.root`
  - objects in rdat root file
    - `c_spin_pat` -- `R_spin = same spin / diff spin`
    - `c_raw_{bbc,zdc,vpd}` = raw scaler counts vs var
    - `c_acc_{bbc,zdc,vpd}` = accidentals corrected scaler counts vs var
    - `c_mul_{bbc,zdc,vpd}` = multiples corrected scaler counts vs var
    - `c_fac_{bbc,zdc,vpd}` = correction factor (mult/raw) vs var
    - `c_R#_{bbc,zdc,vpd}` = relative luminosity vs. var
    - `c_mean_R#` = mean rellum over EWX (bbc,zdc,vpd on same canvas)
    - `c_R#_zdc_minus_vpd` = difference between zdc and vpd
    - `c_deviation_R#_{bbc,zdc,vpd}` = rellum minus mean rellum
    - `rate_dep_R#_{bbc*,zdc*,vpd*}` = rellum vs. multiples corrected rate
    - `c_rate_fac_{bbc*,zdc*,vpd*}` = correction factor vs. multiples corrected rate
                                     (tprofile from `rate_fac` for each spinbit)
  - rellum looping scripts for relative luminosity analysis
    - `rellum_all`
      - basically used to run `rellum4.C` for various independent
        variables etc. 
      - outputs pngs in `png_rellum`, ready to be copied to protected area to link
        to scalers web page
    - `rellum_fills`
      - runs `rellum4.C` for all fills separately and output pdfs in subdirectories of
        `pdf_bXings_fills`; this is for looking at fill dependence of bXing 
        distributions
      - execute `ghost_script_fills` afterward to combine all the pdfs into 
        `pdf_bXings_fills/*.pdf`
      - this was created to search for the origin of pathologies which cause disagreement
        between R3 and mean R3 (and disagreement between ZDC & VPD?)
      - also produces `matrix` trees (see matrix section below)
    - `rellum_runs`
      - analagous to `rellum_fills` but for run-by-run bXing distributions
      - use `ghost_script_runs` for pdf concatenation
      - also produces `matrix` trees (see matrix section below)


6. Combine all the data into a tree to pass to asymmetry analysis
  - run `sumTree.C`, which builds `sums.root` from `counts.root`, which sums the counts
    for each run
    - determines spin pattern types that were collided (see spin pattern
      recognition section below)
    - can now run `nbx_check.C` to test whether the variable `tot_bx` actually
      makes sense with respect to the run time, by plotting `tot_bx/(bXing rate)` vs.
      `run time`; the slope of a linear fit to this should equal unity
  - run `combineAll.C`, which combines `sums.root` and `rdat_i.root` into a final tree, which
    can then be passed to asymmetry analysis code
  - run `make_run_list.C` to make list of run numbers & run indices



Spinbit Equalizing Running 
--------------------------
 - empty bXings are omitted manually in `bunch_kicker` 
 - number of bXings per spinbit is usually unequal (usu. 24,24,26,24); spinbit equalization 
   randomly removes the minimum number of bXings in order to equalize the number of 
   bXings per spinbit
 - to turn on spinbit equalization, in `bunch_kicker`, set `run_randomizer=1`
   - if `run_randomizer!=1`, then no bXings other than empty ones will be kicked
 - `kicked` will now be populated with more bXings to remove; proceed with normal 
   rellum analysis
   - NOTE FOR WEB PAGE: be sure to upload the pngs to protected are in proper directories!
     - `scalers2013/png_rellum` is not spinbit equalized
     - `scalers2013/png_rellum_se` is the spinbit_equalized



Other Useful Scripts
--------------------
 - `draw_spin_patterns.C` -- draws the spin patterns to `pattern.pdf`
 - `draw_fill_vs_run.C` -- draws fill index vs. run index
   - useful for looking for fill structure in plots where 
     run index is the independent variable
 - `nbx_check.C` plots total bXings / bXing rate vs. run time from `sums.root`
   - tau := total bXings / bXing rate (tau should equal run time t)
   - useful to make sure `tot_bx` variable makes sense
   - also looks for runs / fills which have tau != t
 - `nbx_check_2.C` plots the total # bXings vs. bXing no. for each run 
   into a pdf, called `nbx_vs_bxing.pdf` --> odd structure?
   - no satisfactory explanation has been found for this structure
   - Zilong confirmed it in his analysis
   - It's symmetric around bXing 60 for almost all the runs, except for a few (~5)
     which have suddent jumps
   - it's a very small effect and unlikely impacts asymmetry analysis, but its cause
     is not yet understood


Making bXing Distributions
--------------------------

1. Run `accumulate` with no bXings removed
  - echo 0 0 0 > kicked
  - accumulate

2. Draw linear bXing distributions
  - `rellum_fills with drawLog=0 and zoomIn=0`
  - `ghost_script_fills`
  - `mv pdf_bXings_fills/*.pdf htmlfiles/pdf_bXings_fills_lin`

3. Draw linear bXing distributions zooming in on abort gaps
  - `rellum_fills with drawLog=0 and zoomIn=1`
  - `ghost_script_fills`
  - `mv pdf_bXings_fills/*.pdf htmlfiles/pdf_bXings_fills_lin_zoom`

4. Draw logarithmic bXing distributions
  - `rellum_fills with drawLog=1 and zoomIn=0`
  - `ghost_script_fills`
  - `mv pdf_bXings_fills/*.pdf htmlfiles/pdf_bXings_fills_log`



Matrix Subdirectory
-------------------
- Running `rellum4.C` with `var="bx"` and with `specificFill>0` XOR `specificRun>0` will
  produce `matx` tree files, found in `matrix/rootfiles/*.root`
  - this can be done for each fill or run using `rellum_fills` or `rellum_runs`
  - the `matx` tree contains scales, corrected scales, and correction factors for each 
    cbit, tbit, and bXing
- execute `hadd matrix/rootfiles/all.root matrix/rootfiles/matx*.root` to merge the matx trees
- `DrawMatrix.C` draws the desired matrix and produces `matrix.root` and `for_mathematica`
  - `matrix.root` contains the matx tree and the matrix `mat`
  - `for_mathematica` contains the matrix `mat` in text form for reading with mathematica
  - singular value decomposition (SVD) is then performed using `SVD.nb`
- bunch fitting
  - the main code for bunch fitting is `BF.C`, which requires `rootfiles/all.root`, 
    `../counts.root`, and `../sums.root`
  - you need to specify the ratio of scalers to bunch fit over
    - numerator tbit and cbit (see `rellum4.C` for definitions of tbit and cbit)
    - denominator tbit and cbit
    - evaluateChiSquare = true will try to draw Chi2 profiles (not working well yet...)
    - specificPattern != 0 will only consider the specified spin pattern
  - it's best to just use `RunPatterns`, which runs `BF.C` under various interesting
    conditions, for all spin patterns; the following files are produced:
    - `fit_result.[num].[den].root`: bunch fit results for all spin patterns, where the fit
      is done to `r^i=num/den`
    - `pats/fit_result.[num].[den].pat[pat].root`: bunch fit results for spin pattern `pat`
    - `colour.[num].[den].root`: bunch fit results, with colour code according to spin 
      patterns (see the TCanvas `legend` in the ROOT file)



Zilong's Scaler Bit Definitions for Run 13
------------------------------------------

- BBC & ZDC -- 3 bits -- `[coincidence][west][east]`
  - `0` - `000` -- no triggers
  - `1` - `001` -- east
  - `2` - `010` -- west
  - `3` - `011` -- west + east
  - `4` - `100` -- coin
  - `5` - `101` -- coin + east
  - `6` - `110` -- coin + west
  - `7` - `111` -- coin + west + east

- VPD -- 2 bits -- `[west][east]` (no coincidence)
  - `0` - `00` -- no triggers
  - `1` - `01` -- east
  - `2` - `10` -- west
  - `3` - `11` -- west + east



Spin Pattern Recognition
------------------------

- There are 4 types of spin patterns for Run13:
  - pattern 1: + + - - + + - - 
  - pattern 2: - - + + - - + + 
  - pattern 3: + + - - - - + +
  - pattern 4: - - + + + + - -

- each fill is given an "overall" spin pattern no.
  - N := overall spin pattern no.
  - Nb := blue spin pattern no.
  - Ny := yellow spin pattern no.
  - N = 10 * Nb + Ny

- For run13, the following patterns were collided:
  - [13, 14, 23, 24, 31, 32, 41, 42]
