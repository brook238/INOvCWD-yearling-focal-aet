extensions [ csv gis ]          ;MODIFIED FOR SENSITIVITY ANALYSIS 2 MAY 2020 RECORD EACH TRANSMISSION  ;@#@#@#IMP ADDITION IN DEER MATING MODEL; line 305,7,9 cleanup due
                            ;IMP CHANGE IN MATING SUBMODEL - UPDATE CURRENT MODEL;;@#@#@#PostharvestTargetedculling:
globals
[
  dxp                                            ;x coordinate of the patch where a yearling deer belongs before dispersal (natal range)
  dxn                                            ;x coordinate of the patch where a yearling deer reaches after dispersing from its natal range
  dyp                                            ;y coordinate of the patch where a yearling deer belongs before dispersal (natal range)
  dyn                                            ;y coordinate of the patch where a yearling deer reaches after dispersing from its natal range
  dd                                             ;predicted dispersal distance from log-normal pdf
  ndd                                            ;counter-number of deer dispersing out of the model landscape
  tfamh                                          ;counter-fawn male deer hunted
  tymh                                           ;counter-yearling male deer hunted
  tamh                                           ;counter-adult male deer hunted
  tfafh                                          ;counter-fawn female deer hunted
  tyfh                                           ;counter-yearling female deer hunted
  tafh                                           ;counter-adult female deer hunted
  tfmt                                           ;counter-fawn male deer tested
  tymt                                           ;counter-yearling male deer tested
  tamt                                           ;counter-adult male deer tested
  tfft                                           ;counter-fawn female deer tested
  tyft                                           ;counter-yearling female deer tested
  taft                                           ;counter-adult female deer tested
  tgroid                                         ;stores groid (group-id) of an individual during implementation of certain submodels
  ttgroid                                        ;stores groid (group-id) of an individual during implementation of certain submodels
  n_leaders_lost                                 ;counter-doe social group leaders losing their leadership status
  twho                                           ;stores 'who' of an individual during implementation of certain submodels
  counter1                                       ;counter for fawns less than 2 month old while implementing hunting or non-hunting mortality for a female deer
  tgr                                            ;counter for group members during group formation as well as fission
  oldm                                           ;proportion of old males (above 229 when d = 1) in the adult male deer population
  oldf                                           ;proportion of old females (above 229 when d = 1) in the adult female deer population
  tmgroid                                        ;stores mgroid (male group id) during execution of certain submodels
  dcwdm                                          ;counter- CWD infected male adults in the harvest
  dcwdmy                                         ;counter- CWD infected male yearlings in the harvest
  dcwdmf                                         ;counter- CWD infected male fawns in the harvest
  dcwdf                                          ;counter- CWD infected female adults in the harvest
  dcwdfy                                         ;counter- CWD infected female yearlings in the harvest
  dcwdff                                         ;counter- CWD infected female fawns in the harvest
  tcwdm                                          ;counter- CWD infected male adults in the sample
  tcwdmy                                         ;counter- CWD infected male yearlings in the sample
  tcwdmf                                         ;counter- CWD infected male fawns in the sample
  tcwdf                                          ;counter- CWD infected female adults in the sample
  tcwdfy                                         ;counter- CWD infected female yearlings in the sample
  tcwdff                                         ;counter- CWD infected female fawns in the sample
  pdcwd                                          ;probability of detecting CWD
  mom                                            ;mom's who
  ttmomid                                        ;for identifying adult siblings
  tmfh                                           ;counter-male fawn deer harvested
  tffh                                           ;counter-female fawn deer harvested
  tmyh                                           ;counter-male yearling deer harvested
  tfyh                                           ;counter-female yearling deer harvested
  d
  year
  ;MOOvPOP( Agent-based model of deer population dynamics) generated deer population is used to initialize this model ('import-world'),
  ;hence global variables from MOOvPOP are also included.
  vals                                          ;list to store reporters for output file
  vals1                                         ;list to store reporters for output file
  vals2                                         ;list to store reporters for output file
  vals3
  output2
  region
  recommended_parameter_values
  post_harvest_density
  sexratio
  adultprop
  yearlingprop
  nt
  ntplus1
  lambda
  max-dist-spark                                ;farthest inf patch from the initial focus
  contactmatrix
  ;first-detection
  cum-trans                                  ;2May2020 for sa record each transmission event
  iteration
  ;@#@#@#
  tc-cwdmf                                   ;@#@#@# target culled -cwd positive male fawns
  tc-tmf                                     ;@#@#@# target culled and tested male fawns (all target culled animals are tested)
  tc-cwdff
  tc-tff
  tc-cwdmy
  tc-tmy
  tc-cwdfy
  tc-tfy
  tc-cwdm
  tc-tm
  tc-cwdf
  tc-tf
  ttc
  tc-cwd
  forest-cover-raster ;20220322. Stores forest cover raster used to set projection
  output-raster ;20220322. Stores raster of infected patches
  ]
patches-own
[
  forest-percent
  do                                             ;deer occupancy 0/1male/2female)
  border                                         ;identifying border patches
  dfp
  dh
  cwd-d                                       ;cwd detected @#@#@#
  infected? ;20220322. Stores value output to CWD distribution rasters. 0 = No CWD infected deer, 1 = At least one infected deer ever occurred in patch.
  ]
breed [ deers deer ]
deers-own
[
  sex
  aim                                            ;age in months
  cwd                                            ;infection status, 1 if infected, 0 if uninfected
  cwdm                                           ;duration of CWD - time from exposure to death
  cwdi                                           ;duration of pre-infectious phase - time from exposure to the onset of infectious phase
  cwdc                                           ;duration of pre-clinical phase - time from exposure to onset of clinical signs
  cwdpr                                          ;cwd progression (in months)
  momid                                          ;mother's who number
  gl                                             ;1 group leader; 0 follower
  groid
  gr                                             ;basically a leader attribute-group size excluding the leader (0 to n); if -1 the deer (gl = 0) is in a group; -2 if the deer is not in a group
  ml                                             ;male bachelor group leader ml = 1 otherwise 0
  mgroid                                         ;male bachelor group -2 at birth, -1 after dispersal and = who of leader when join group as adults
  nm                                             ;max number of matings per doe during rut period (between 1 to 3)
  anm                                            ;counter for actual number of matings
  ]

to setup
  clear-all
  setup-landscape
  reset-ticks
  ask deers with [ groid > 0 ] [
    if distance deer groid > 3 [
      set groid -1
      set gr -2
      ]
    ]
  ask deers [
    ;set cwdm 20 + random 6                      ;total duration of the disease in a deer - from infection until mortality
    set cwdi 6 + random 5                     ;7 + random5;6 + random5 to 5 to 5 sa3  ;3 + random 2 ;duration after which a deer becomes infectious - from exposure to becoming infectious
    set cwdc 21 + random 5                    ;22 + random 5; 21 + random 5 to 20 + random 5  ; 15 + random 4;exposure to overt clinical signs in GG- derived from Johnson et al., 2011 duration after which a deer starts exhibiting clinical signs of CWD
    set cwdm cwdc + random 2
    ]
  set contactmatrix csv:from-file "../data/contactstructure.csv"
  if (file-exists? (word "../results/cwdinfdy" cwd_region "_seedinf_" seed-infection ".csv") = FALSE) [
    file-open (word "../results/cwdinfdy" cwd_region "_seedinf_" seed-infection ".csv")
    file-print "Year,AdultMale,YearlingMale,FawnMale,AdultFemale,YearlingFemale,FawnFemale,TotalPreHarvestPop,AdultMaleCWD,YearlingMaleCWD,FawnMaleCWD,AdultFemaleCWD,YearlingFemaleCWD,FawnFemaleCWD,TotalCWD,AdultMaleHarvest,YearlingMaleHarvest,FawnMaleHarvest,AdultFemaleHarvest,YearlingFemaleHarvest,FawnFemaleHarvest,AdultMaleTested,YearlingMaleTested,FawnMaleTested,AdultFemaleTested,YearlingFemaleTested,FawnFemaleTested,CWDhuntedDeer,CWDtestedDeer,DetProb,ObsPrev,CWDArea,MaxDistSpark,cum-trans,iteration,TotalTargetedCull,TCCWD+";#2May2020 added cum-trans#7Mar18added CWDArea
    file-close
    ]
  ;----------------------------------------------------------------------
;Uncomment the following lines to write contactstructure.csv file. Update values and execute setup.
;  csv:to-file "contactstructure.csv" [["age" "mom-min" "mom-max" "mf-min" "mf-max" "ff-min" "ff-max" "mfullsib-min" "mfullsib-max" "ffullsib-min" "ffullsib-max" "mnonsib-min" "mnonsib-max" "fnonsib-min" "fnonsib-max" "grdoe-min" "grdoe-max" "ndoe-min" "ndoe-max" "nbuck-min" "nbuck-max"]
;                                 ["3 mo" 60 90 0 0 0 0 60 90 60 90 0 0 0 0 5 10 0 0 0 0]
;    ["4-6 mo male" 10 20 0 0 0 0 30 50 30 50 10 20 5 10 5 10 0 0 0	0]
;    ["4-6 mo female" 20	30 0 0 0 0 30 50 30	50 5 10	20 30	10 20	0	0	0	0]
;    ["7-8 mo male" 0 5 0 0 0 0 10 20 0 5 10 20 0 5 0 5 0 0 0 0]
;    ["7-8 mo female" 5 10	0	0	0	0	0	5	10 20	0	5	10 20	0	5	0	0	0	0]
;    ["9-12 mo male" 5	10 0 0 0 0 10 20 5 10 10 20 5 10 0 5 0 0 0 0]
;    ["9-12 mo female" 10 20	0	0	0	0	5	10 10	20 5 10	10 20	5	10 0 5 0 0]
;    ["13-14 mo male" 0 0 0 0 0 0 10	20 5 10	10 20	5	10 0 0 5 10	0	0]
;    ["13-14 mo female nf" 0	0	0	0	0	0	5	10 10	20 5 10	10 20	0	0	5	10 0 0]
;    ["13-14 mo female wf" 0	0	60 90	60 90	0	0	0	0	0	0	0	0	0	0	0	0	0	0]
;    ["15 mo female" 10 20	60 90	60 90	0	0	10 20	0	0	0	0	5	10 0 5 0 0]
;    ["16-18 mo female" 10	20 10	20 20	30 0 0 10	20 0 0 0 0 5 10	0	5	0	0]
;    ["19-20 mo female" 0 5 0 5 5 10	0	0	0	5	0	0	0	0	5	10 5 8 0 0]
;    ["21-24 mo female" 10	20 5 10	10 20	0	0	10 20	0	0	0	0	5	20 0 5 0 0]
;    ["25-26 mo female nf" 0	0	0	0	0	0	0	0	5	10 0 0 0 0 5 10	0	0	5	10]
;    ["25-26 mo female wf" 0	0	60 90	60 90	0	0	0	0	0	0	0	0	0	0	0	0	0	0]
;    [">26 mo female (gest.period)" 0 10	5	10 10	20 0 0 0 0 0 0 0 0 5 10	0	5	0	0]
;    [">26 mo female (weaning period)" 0	10 60	90 60	90 0 0 0 0 0 0 0 0 5 10	0	5	0	0]
;    [">26 mo female (prerut period)" 0 10	10 20	20 30	0	0	0	0	0	0	0	0	5	10 0 5 0 0]
;    [">26 mo female (ryt period)" 0 5 0 5 5 10	0	0	0	0	0	0	0	0	0	5	5	8	0	0]
;    [">14 mo male b-gr" 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	5	15]
;    [">14 mo male solitary" 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	5]
;  ]
  ;------------------------------------------------------------------------------
  ;set myhm 0.52   ;APR2 (increase juv harvest)
  ;set myhm 0.33   ;APR1 (decrease juv male harvest)0.32
  set forest-cover-raster gis:load-dataset "../data/kankakee_20220327.asc" ;20220322. Load raster used to set projection.
  gis:set-world-envelope gis:envelope-of forest-cover-raster ;20220322. Set projection.
end
;------------------------------------------------------------
to setup-landscape
  if (cwd_region = "LinnMaconMO") [ import-world "data/PostHarvestPopulationLinnMacon_v2.2.0.csv" ]
  if (cwd_region = "Kankakee") [ import-world "../data/PostHarvestPopulationKankakee_v2.2.0.csv" ]
  ask deers [ set hidden? TRUE ]
  random-seed new-seed
  ask patches [ set infected? 0 ] ;20220322. Sets infected? to zero.
end
;------------------------------------------------------------
to go
  if ticks = 0 [
    set iteration (behaviorspace-run-number)                                  ;5Oct21 added for HPCC runs
    ]
  if (ticks = 300) [ stop ]       ;model is simulated for 10 years set at = 120. To run for 25 years, set at =300
  set d remainder ticks 12 + 1    ;month
  set year floor (ticks / 12) + 1 ;year
  if (d = 1 or d = 7) [
    ask deers [ ht ]
    ;export-view (word "CWD landscape month " ticks ".png")  ;can also export interface
    ]
;------------------------------------------------------------
;uncomment the following lines if you want to stop the model when there is no CWD+ deer in the model landscape
;  let cwdpd count deers with [ cwd = 1 ]
;  if (cwdpd = 0) [
;    user-message (word "No CWD+ deer in the model landscape!")
;    stop
;    ]
;------------------------------------------------------------
  ask deers [
    individual-growth             ;turtle procedure deer age
    deer-die                      ;turtle procedure deer natural mortality
    if cwd = 1 [
      cwd-progression             ;course of CWD in infected deer
      ]
    ]
  ask deers with [ sex = 2 and gl = 1 and gr < 2 ] [ ;doe social group dynamics
    set tgroid groid
    set gl 0
    set groid -1
    ask deers with [ groid = tgroid ] [
      set gr -2
      set groid -1
      set n_leaders_lost n_leaders_lost + 1
      ]
    ]
  if (d = 1) [
    form-bachelor-groups
    set oldm precision (count deers with [ sex = 1 and aim >= 229 ] / ma) 3
    set oldf precision (count deers with [ sex = 2 and aim >= 229 ] / fa) 3
    ]
  if (d > 1 and d < 10) [       ;turtle procedure: bachelor group leaders assess group membership
    ask deers with [ sex = 1 and ml = 1 and gr <= 1 ] [
      set ml 0
      ]
    ]
  if (d = 10) [                                                     ;turtle procedure: bachelor groups break down before the rutting season
    let male-leaders deers with [ sex = 1 and ml = 1 ]
    ask male-leaders [ set gr 0 ]
    ]

  if (d = 5) [
    let male-yearlings deers with [ aim = 13 and sex = 1 ]
    ask male-yearlings [                                             ;turtle procedure: dispersal of male yearlings
      set tgroid groid
      set counter1 0
      if (gr = -1) [ review-group-dynamics ]
      set gr -2
      set groid -1
      if (random-float 1 < yearling-male-dispersal-rate) [
        set mgroid -1
        deer-mdisperse
        ]
      ]
    let female-yearlings deers with [ aim = 13 and sex = 2 ]
    ask female-yearlings [                                            ;turtle procedure: dispersal of female yearlings
      if (random-float 1 < yearling-female-dispersal-rate) [
        set tgroid groid
        set counter1 0
        if (gr = -1) [ review-group-dynamics ]
        set gr -2
        set groid -1
        deer-fdisperse
        ]
      ]
    let breeding-females deers with [ sex = 2 and aim > 12 and cwdpr < cwdc ]
    ask breeding-females [
      fawning
      ]
    ]
;-----------------CWD introduction occurs in the 1st year of model run/CWD spread calculated
  if d = 6 [
    if year = 1 [
      let inf-focus nobody
      set inf-focus patches with [ do = 1 ]
      ask n-of 1 inf-focus [
        set pcolor orange
        let cand-deers []
        if CWD_introduced_by = "adult-deer" [
          set cand-deers deers in-radius 2.5 with [ aim > 24 ]                  ;adult deer in the vicinity (2.5 mile radius) of the selected patch
          ]
        if CWD_introduced_by = "dispersing-male-yearling" [
          set cand-deers deers with [ color = blue and aim < 15 and sex = 1 ]   ;dispersing male yearling
          ]
        if CWD_introduced_by = "dispersing-female-yearling" [
          set cand-deers deers with [ color = blue and aim < 15 and sex = 2 ]   ;dispersing female yearling
          ]
        if CWD_introduced_by = "doe-groupmember" [
          set cand-deers deers with [ sex = 2 and aim > 24 and groid >= 0 ]     ;adult doe belonging to doe social group
          ]
        if CWD_introduced_by = "doe-solitary" [
          set cand-deers deers with [ sex = 2 and aim > 24 and groid < 0 ]      ;adult solitary does
          ]
        if CWD_introduced_by = "buck-groupmember" [
          set cand-deers deers with [ sex = 1 and aim > 24 and mgroid >= 0 ]    ;adult male member of a bachelor group
          ]
        if CWD_introduced_by = "buck-solitary" [
          set cand-deers deers with [ sex = 1 and aim > 24 and mgroid < 0 ]     ;adult male solitary
          ]
        ask min-n-of seed-infection cand-deers [ distance-nowrap myself ] [
          set cwd 1
          set cwdpr (5 + random 6)
          ]
        ]
      ]
    distance-to-farthest-spark
    ]
  if (d = 11) [                                                 ;turtle procedure: male yearling dispersal before the rutting season
    let solitary-male-yearlings deers with [ sex = 1 and aim = 19 and gr = -2 and mgroid = -2 ]
    ask solitary-male-yearlings [
      if (random-float 1 <= yearling-male-dispersal-rate) [
        set mgroid -1
        deer-mdisperse
        ]
      set mgroid -1
      ]
    if year >= 2 [                                                ;to calculate lambda
      set nt ntplus1
      ]
;----------------------------deer mating submodel----------------
    let female-deers deers with [ sex = 2 ]
    ask female-deers [
      set nm (1 + random 3)
      set anm 0
      ]
    let breeding-males deers with [ sex = 1 and aim > 18 and cwdpr < cwdc]
    ask breeding-males [
      deer-mating
      ]
    ]
  if (d = 12) [
    let vals4 0
    let vals5 0
    let hcwd 0
    ;let hcwd-sr 0                              ;@#@#@#5Jan20 uncommented - delete after cleanup
    let tcwd 0
    ;let tcwd-sr 0                              ;@#@#@#5Jan20 uncommented - delete after cleanup
    let op 0
    ;let op-sr 0                                ;@#@#@#5Jan20 uncommented - delete after cleanup
    let phn (mf + my + ma + ff + fy + fa)
    set ntplus1 phn                                                   ;to calculate lambda
    ifelse year >= 2
    [ set lambda precision (ntplus1 / nt) 3 ]
    [ set lambda "NA" ]
    set vals1 (list (year) (ma) (my) (mf) (fa) (fy) (ff) (phn) (mcwd) (mycwd) (mfcwd) (fcwd) (fycwd) (ffcwd) (totcwdd))
    ask deers [
      hunting-mortality
      ]
    post-harvest-targeted-culling                           ;Postharvest targeted culling@#@#@#
    set ttc (tc-tmf + tc-tff + tc-tmy + tc-tfy + tc-tm + tc-tf)  ;@#@#@#
    set tc-cwd (tc-cwdmf + tc-cwdff + tc-cwdmy + tc-cwdfy + tc-cwdm + tc-cwdf)  ;@#@#@#

    set hcwd (dcwdm + dcwdmy + dcwdmf + dcwdf + dcwdfy + dcwdff)
    set tcwd (tcwdm + tcwdmy + tcwdmf + tcwdf + tcwdfy + tcwdff)
    if (tcwd > 0) [
      set pdcwd 1
      ]
    set op precision (tcwd / (tamt + tymt + tfmt + taft + tyft + tfft)) 3
    set vals2 (list (tamh) (tymh) (tfamh) (tafh) (tyfh) (tfafh) (tamt) (tymt) (tfmt) (taft) (tyft) (tfft) (hcwd) (tcwd) (pdcwd) (op) (cwd_area) (max-dist-spark) (cum-trans) (iteration) (ttc) (tc-cwd));#2May2020 added cum-trans#7Mar18 added cwd-area
    set vals (sentence vals1 vals2)
    file-open (word "../results/cwdinfdy" cwd_region "_seedinf_" seed-infection ".csv")
    file-type first vals
    foreach but-first vals [ [ ?1 ] ->
      file-type "," file-type ?1
      ]
    file-print""
    file-close
    ]
  set ndd 0
  if (ticks < 299) [                               ;to run the model for 10 years, set at < 119. To run the model for 25 years aet to <299
    set tamh 0
    set tamt 0
    set tymh 0
    set tymt 0
    set tfamh 0
    set tfmt 0
    set tafh 0
    set taft 0
    set tyfh 0
    set tyft 0
    set tfafh 0
    set tfft 0
    set dcwdm 0
    set dcwdmy 0
    set dcwdmf 0
    set dcwdf 0
    set dcwdfy 0
    set dcwdff 0
    set tcwdm 0
    set tcwdmy 0
    set tcwdmf 0
    set tcwdf 0
    set tcwdfy 0
    set tcwdff 0
    ;@#@#@#
    set tc-cwdmf 0
    set tc-tmf 0
    set tc-cwdff 0
    set tc-tff 0
    set tc-cwdmy 0
    set tc-tmy 0
    set tc-cwdfy 0
    set tc-tfy 0
    set tc-cwdm 0
    set tc-tm 0
    set tc-cwdf 0
    set tc-tf 0
    ]
  set-current-plot "deer population"
  if (ticks = 0) [ clear-plot set-plot-pen-color blue set-plot-x-range 0 130 ]
  plotxy ticks count deers
  ;uncomment the following lines if only patches where CWD+ deer occur at present are to be highlighted
;  ask patches with [ pcolor = yellow ] [
;    if not any? deers-here with [ cwd = 1 ] [
;      set pcolor scale-color green forest-percent 1 0
;      set plabel 0
;      ]
;    ]
  if d = 12 [
    set-current-plot "doe matings"
    let mated-does deers with [ sex = 2 and anm > 0 ]
    histogram [ anm ] of mated-does
    ;print (word "Proportion one mating " precision (count mated-does with [ anm = 1 ] / count mated-does) 2 )
    ;print (word "Proportion two matings " precision (count mated-does with [ anm = 2 ] / count mated-does) 2 )
    ;print (word "Proportion three matings " precision (count mated-does with [ anm = 3 ] / count mated-does) 2 )
  ]

  ask patches with [ dh != -1 ] [
    let pop count deers-here
    let prev 0
    if pop > 0 [ set prev (count deers-here with [ cwd = 1 ] ) / pop ]
    ifelse prev > 0
    [ ifelse pcolor != orange
      [ set pcolor yellow
        set plabel-color 15
        ]
      [ set plabel-color black
        ]
      set plabel count deers-here with [ cwd = 1 ]
      ]
    [ set plabel ""
      ]
    ]

  ;20220322. Updates infected? if patch is newly infected
  ask patches with [(pcolor = yellow or pcolor = orange) and infected? = 0] [
    set infected? 1
  ]

  ;20220322. If month = January, export raster with infected patches
  if export-rasters = true [
    if d = 1 [
      ask patches [
      set output-raster gis:patch-dataset infected?
      ]
    gis:store-dataset output-raster (word "../results/cwd_dist_" year ".asc")
    ]
  ]

  tick
end

to individual-growth                                         ;turtle procedure: age (age in months) increases by 1 time step
  set aim (aim + 1)
end

to cwd-progression
  set cwdpr cwdpr + 1
  if (cwdpr >= cwdm) [ deer-die-cwd ]
  if (cwdpr > cwdi) and (cwdpr < cwdc) [                           ;infectious phase, once an infected deer starts to show clinical signs, we assume no direct contact with susceptible individuals
    let local-deers deers in-radius-nowrap 3                       ;immediate limitation to local deer
    let close-deers local-deers in-radius-nowrap 1.5               ;even closer
    let lgroid groid
    let lmom momid
    if aim = 3 [                                                   ;fawn 3 months old
      ask close-deers with [ who = lmom and cwd = 0 ] [            ;uninfected mom 60-90 contacts ;@#@#@#a3mom
        if random-float 1 < a3mom [
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
      ask close-deers with [ momid = lmom and aim = 3 and cwd = 0 ] [     ;full sibs in cohort > 70 contacts
        if random-float 1 < a3fs [
          set cwd 1                                                       ;@#@#@#a3fs
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
      if (lgroid >= 0) [                                                  ;other female group members
        ask local-deers with [ groid = lgroid and sex = 2 and cwd = 0 and who != lmom and aim > 3 ] [  ;uninfected yearling and adult females, but not mom or full sibs 5-10 contacts
          if random-float 1 < a3grf [                                     ;@#@#@#a3grf
            set cwd 1
            set cum-trans cum-trans + 1                         ;#2May2020 for sa
            ]
          ]
        ]
      ]
    if aim > 3 and aim < 7 [                                      ;fawn 4, 5 or 6 months old
      ifelse sex = 1                                              ;CWD+ fawn - male/female
      [ ask close-deers with [ who = lmom and cwd = 0 ] [         ;susceptible mom
        if (random-float 1 < a456mmom) [                          ;10-20 contacts ;@#@#@#a456mmom
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
        ask close-deers with [ momid = lmom and aim < 7 and cwd = 0 ] [ ;susceptible full sibs
          if (random-float 1 < a456mfs) [                               ;30-50 contacts ;@#@#@#a456mfs
            set cwd 1
            set cum-trans cum-trans + 1                         ;#2May2020 for sa
            ]
          ]
        if lgroid >= 0 [                                           ;only for fawns in a doe social group
          ask local-deers with [ groid = lgroid and who != lmom and momid != lmom and cwd = 0 ] [  ;susceptible group members, excluding mom and full sibs
            if aim < 7 [
              ifelse sex = 1
              [ if (random-float 1 < a456mnsmf) [                  ;10-20 contacts    male non sib gr memb;@#@#@#a456mnsmf
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
                ]
              [ if random-float 1 < a456mnsff [                    ;5-10 contacts   female non sib gr memb;@#@#@#a456mnsff
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
                ]
              ]
            if aim > 15 [
              if random-float 1 < a456mgrf [                       ;5-10 contacts   yearling and adult female gr memb;@#@#@#a456mgrf
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
              ]
            ]
          ]
        ]
      [ ask close-deers with [ who = lmom and cwd = 0 ] [          ;CWD+ female fawn 4,5,6; infecting susceptible mom
        if (random-float 1 < a456fmom) [                           ;20-30 contacts;@#@#@#a456fmom
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
        ask close-deers with [ momid = lmom and aim < 7 and cwd = 0 ] [ ;susceptible full sibs
          if (random-float 1 < a456ffs) [                               ;30-50 contacts;@#@#@#a456ffs
            set cwd 1
            set cum-trans cum-trans + 1                         ;#2May2020 for sa
            ]
          ]
        if lgroid >= 0 [
          ask local-deers with [ groid = lgroid and who != lmom and momid != lmom and cwd = 0 ] [  ;susceptible group members, excluding mom and full sibs
            if aim < 7 [
              ifelse sex = 1
              [ if (random-float 1 < a456fnsmf) [                   ;5-10 contacts    male non sib gr memb;@#@#@#a456fnsmf
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
                ]
              [ if random-float 1 < a456fnsff [                     ;20-30 contacts   female non sib gr memb;@#@#@#a456fnsff
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
                ]
              ]
            if aim > 15 [
              if random-float 1 < a456fgrf [                         ;10-20 contacts   yearling and adult female gr memb;@#@#@#a456fgrf
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
              ]
            ]
          ]
        ]
      ]
    if aim = 7 or aim = 8 [
      ifelse sex = 1                                             ;male CWD+ age 7 or 8 mo
      [ ask close-deers with [ who = lmom and cwd = 0 ] [        ;CWD+ male fawn 7,8 mo infecting susceptible mom
        if random-float 1 < a78mmom [                            ;0-5 contacts;@#@#@#a78mmom
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
        ask close-deers with [ momid = lmom and aim >= 7 and aim < 9 and cwd = 0 ] [ ;susceptible full sibs
          ifelse sex = 1
          [ if random-float 1 < a78mfsm [                        ;10-20 contacts;@#@#@#a78mfsm
            set cwd 1
            set cum-trans cum-trans + 1                         ;#2May2020 for sa
            ]
          ]
          [ if random-float 1 < a78mfsf [                        ;0-5 contacts;@#@#@#a78mfsf
            set cwd 1
            set cum-trans cum-trans + 1                         ;#2May2020 for sa
            ]
            ]
        ]
        if lgroid >= 0 [
          ask local-deers with [ groid = lgroid and who != lmom and momid != lmom and cwd = 0 ] [  ;susceptible group members, excluding mom and full sibs
            if aim = 7 or aim = 8 [
              ifelse sex = 1
              [ if (random-float 1 < a78mnsmf) [                  ;10-20 contacts    male non sib gr memb;@#@#@#a78mnsmf
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
                ]
              [ if random-float 1 < a78mnsff [                    ;0-5 contacts   female non sib gr memb;@#@#@#a78mnsff
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
                ]
              ]
            if aim > 18 [
              if random-float 1 < a78mgrf [                       ;0-5 contacts   yearling and adult female gr memb;@#@#@#a78mgrf
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
              ]
            ]
          ]
        ]
      [ ask close-deers with [ who = lmom and cwd = 0 ] [         ;CWD+ female fawn 7,8 mo infecting susceptible mom
        if random-float 1 < a78fmom [                             ;5-10 contacts ;@#@#@#a78fmom
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
        ask close-deers with [ momid = lmom and aim > 6 and aim < 9 and cwd = 0 ] [ ;susceptible full sibs
          ifelse sex = 1
          [ if random-float 1 < a78ffsm [                         ;0-5 contacts  ;@#@#@#a78ffsm
            set cwd 1
            set cum-trans cum-trans + 1                         ;#2May2020 for sa
            ]
          ]
          [ if random-float 1 < a78ffsf [                          ;10-20 contacts ;@#@#@#a78ffsf
            set cwd 1
            set cum-trans cum-trans + 1                         ;#2May2020 for sa
            ]
            ]
          ]
        if lgroid >= 0 [
          ask local-deers with [ groid = lgroid and who != lmom and momid != lmom and cwd = 0 ] [  ;susceptible group members, excluding mom and full sibs
            if aim = 7 or aim = 8 [
              ifelse sex = 1
              [ if random-float 1 < a78fnsmf [                     ;0-5 contacts    male non sib gr memb ;@#@#@#a78fnsmf
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
                ]
              [ if random-float 1 < a78fnsff [                     ;10-20 contacts   female non sib gr memb;@#@#@#a78fnsff
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
                ]
              ]
            if aim > 18 [
              if random-float 1 < a78fgrf [                        ;0-5 contacts   yearling and adult female gr memb ;@#@#@#a78fgrf
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
              ]
            ]
          ]
        ]
      ]
    if aim > 8 and aim < 13 [                                       ;CWD+ fawns age 9, 10, 11, 12
      ifelse sex = 1
      [ ask close-deers with [ who = lmom and cwd = 0 ] [           ;CWD+ male fawn 9-12 mo infecting susceptible mom
        if random-float 1 < a912mmom [                              ;5-10 contacts;@@##@@##a912mmom
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
        ask close-deers with [ momid = lmom and aim >= 9 and aim < 13 and cwd = 0 ] [ ;susceptible full sibs
          ifelse sex = 1
          [ if random-float 1 < a912mfsm [                           ;10-20 contacts;@#@#@#a912mfsm
            set cwd 1
            set cum-trans cum-trans + 1                         ;#2May2020 for sa
            ]
          ]
          [ if random-float 1 < a912mfsf [                           ;5-10 contacts;@#@#@#a912mfsf
            set cwd 1
            set cum-trans cum-trans + 1                         ;#2May2020 for sa
            ]
            ]
          ]
        if lgroid >= 0 [
          ask local-deers with [ groid = lgroid and who != lmom and momid != lmom and cwd = 0 ] [  ;susceptible group members, excluding mom and full sibs
            if aim > 8 and aim < 13 [
              ifelse sex = 1
              [ if (random-float 1 < a912mnsmf) [                      ;10-20 contacts    male non sib gr memb;@#@#@#a912mnsmf
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
                ]
              [ if random-float 1 < a912mnsff [                        ;5-10 contacts   female non sib gr memb;@#@#@#a912mnsff
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
                ]
              ]
            if aim > 20 [
              if random-float 1 < a912mgrf [                           ;0-5 contacts   yearling and adult female gr memb;@#@#@#a912mgrf
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
              ]
            ]
          ]
        ]
      [ ask close-deers with [ who = lmom and cwd = 0 ] [               ;CWD+ female fawn 9-12 mo infecting susceptible mom
        if random-float 1 < a912fmom [                                  ;10-20 contacts;@#@#@#a912fmom
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
        ask close-deers with [ momid = lmom and aim >= 9 and aim < 13 and cwd = 0 ] [ ;susceptible full sibs
          ifelse sex = 1
          [ if random-float 1 < a912ffsm [                               ;5-10 contacts;@#@#@#a912ffsm
            set cwd 1
            set cum-trans cum-trans + 1                         ;#2May2020 for sa
            ]
          ]
          [ if random-float 1 < a912ffsf [                               ;10-20 contacts;@#@#@#a912ffsf
            set cwd 1
            set cum-trans cum-trans + 1                         ;#2May2020 for sa
            ]
            ]
        ]
        if lgroid >= 0 [
          ask local-deers with [ groid = lgroid and who != lmom and momid != lmom and cwd = 0 ] [  ;susceptible group members, excluding mom and full sibs
            if aim > 8 and aim < 13 [
              ifelse sex = 1
              [ if (random-float 1 < a912fnsmf) [                          ;5-10 contacts    male non sib gr memb;@#@#@#a912fnsmf
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
                ]
              [ if random-float 1 < a912fnsff [                            ;10-20 contacts   female non sib gr memb;@#@#@#a912fnsff
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
                ]
              ]
            if aim > 20 [
              if random-float 1 < a912fgrf [                               ;5-10 contacts   yearling and adult female gr memb;@#@#@#a912fgrf
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
              ]
            ]
        ]
        if lgroid >= 0 [                               ;if the infected female fawn is a group member, it can also transmit the infection to one of the non-group female yearling/adult in the neighborhood
          if any? local-deers with [ sex = 2 and aim > 20 and groid != lgroid and who != lmom and cwd = 0 ] [
            ask n-of 1 local-deers with [ sex = 2 and aim > 20 and groid != lgroid and who != lmom and cwd = 0 ] [
              if random-float 1 < a912fngrf [                                 ;0-5 contacts ;@#@#@#a912fngrf
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
              ]
            ]
          ]
        if lgroid = -1 [                              ;if the infected female fawn is solitary, it can also infect one of the yearling/adult in the neighborhood
          if any? local-deers with [ sex = 2 and aim > 20 and who != lmom and momid != lmom and cwd = 0 ] [
            ask n-of 1 local-deers with [ sex = 2 and aim > 20 and who != lmom and momid != lmom and cwd = 0 ] [
              if random-float 1 < a912fngrf [                                 ;0-5 contacts ;@#@#@#a912fngrf
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
              ]
            ]
          ]
        ]
      ]
    if aim > 12 and aim < 15 [                                               ;yearling age 13, 14
      ifelse sex = 1                                                         ;male
      [ ask close-deers with [ momid = lmom and aim > 12 and aim < 15 and cwd = 0 ] [ ;susceptible full sibs;
        ifelse sex = 1                                                       ;male/female
        [ if random-float 1 < a1314mfsm [                                    ;10-20 contacts ;fullsib male;@#@#@#a1314mfsm
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
          ]
        [ let myid who
          if not any? close-deers with [ momid = myid ] [                    ;female full sib without fawns
            if random-float 1 < a1314mfsf [                                  ;5-10 contacts;@#@#@#a1314mfsf
              set cwd 1
              set cum-trans cum-trans + 1                         ;#2May2020 for sa
              ]
            ]
          ]
        ]
        if lgroid >= 0 [
          let male-year-here local-deers with [ cwd = 0 and aim > 12 and aim < 15 and momid != lmom and sex = 1 ]   ;susceptible non-sib male yearlings in the vicinity
          let female-year-here local-deers with [ cwd = 0 and aim > 12 and aim < 15 and momid != lmom and sex = 2 ] ;susceptible non-sib female yearlings in the vicinity
          let female-ad-here local-deers with [ cwd = 0 and aim > 24 and aim < 27 and sex = 2 ] ;susceptible young adult females in the vicinity
          let num-year-mcontact 0
          let num-year-fcontact 0
          let num-ad-fcontact 0
          let counter11 0
          while [ count male-year-here > 0 and num-year-mcontact < 4 ] ;contacts upto 3 other male yearlings in the vicinity
          [ ask n-of 1 male-year-here [
            if random-float 1 < a1314mnsmy [                           ;10-20 contacts;@#@#@#a1314mnsmy
              set cwd 1
              set cum-trans cum-trans + 1                         ;#2May2020 for sa
            ]
            set num-year-mcontact num-year-mcontact + 1
            ]
          ]
          while [ count female-year-here > 0 and num-year-fcontact < 4 and counter11 < 6 ]        ;contacts upto 3 other female yearlings without fawns in the vicinity
          [ ask n-of 1 female-year-here [
            let myid who
            if not any? close-deers with [ momid = myid ] [            ;nonsib yearling female without fawns
              if random-float 1 < a1314mnsfy [                         ;5-10 contacts;@#@#@#a1314mnsfy
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
              set num-year-fcontact num-year-fcontact + 1
              ]
            set counter11 counter11 + 1
            ]
          ]
          set counter11 0
          while [ count female-ad-here > 0 and num-ad-fcontact < 4 and counter11 < 6 ]             ;contacts upto 3 other young ad females without fawns in the vicinity
          [ ask n-of 1 female-ad-here [
            let myid who
            if not any? close-deers with [ momid = myid ] [            ;nonmom young adult female without fawns
              if random-float 1 < a1314mngrf [                         ;5-10 contacts;@#@#@#a1314mngrf
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
              set num-ad-fcontact num-ad-fcontact + 1
              ]
            set counter11 counter11 + 1
            ]
            ]
          set counter11 0
          ]
        ]
      [ let myid who
        ifelse not any? close-deers with [ momid = myid ]                ;female yearling without/with fawns
        [ ask close-deers with [ momid = lmom and aim > 12 and aim < 15 and cwd = 0 ] [ ;susceptible full sibs
          ifelse sex = 1
          [ if random-float 1 < a1314ffsm [                              ;5-10 contacts;@#@#@# a1314ffsm
            set cwd 1
            set cum-trans cum-trans + 1                         ;#2May2020 for sa
            ]
            ]
          [ let myid1 who
            if not any? close-deers with [ momid = myid1 ] [             ;female full sib without fawns
              if random-float 1 < a1314ffsf [                            ;10-20 contacts;@#@#@#a1314ffsf
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
              ]
            ]
          ]
          if lgroid >= 0 [
            let male-year-here local-deers with [ cwd = 0 and aim > 12 and aim < 15 and momid != lmom and sex = 1 ]   ;susceptible non-sib male yearlings in the vicinity
            let num-year-mcontact 0
            ;let counter11 0
            while [ count male-year-here > 0 and num-year-mcontact < 4 ]   ;contacts upto 3 other male yearlings in the vicinity
            [ ask n-of 1 male-year-here [
              if random-float 1 < a1314fnsmy [                             ;5-10 contacts;@#@#@#a1314fnsmy
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
              set num-year-mcontact num-year-mcontact + 1
              ]
              ]
            ask local-deers with [ groid = lgroid and who != lmom and momid != lmom and sex = 2 and cwd = 0 ] [  ;susceptible female group members, excluding mom and full sibs
              if aim > 12 and aim < 15 [
                let myid1 who
                if not any? close-deers with [ momid = myid1 ] [          ;non sib female yearling in group without fawns
                  if random-float 1 < a1314fnsfy [                        ;10-20 contacts;@#@#@#a1314fnsfy
                    set cwd 1
                    set cum-trans cum-trans + 1                         ;#2May2020 for sa
                    ]
                  ]
                ]
              if aim > 24 and aim < 27 [                                  ;young females in group age 25 and 26
                let myid1 who
                if not any? close-deers with [ momid = myid1 ] [          ;without fawns
                  if random-float 1 < a1314fgrf [                         ;5-10 contacts   yearling and adult female gr memb;@#@#@#a1314fgrf
                    set cwd 1
                    set cum-trans cum-trans + 1                         ;#2May2020 for sa
                    ]
                  ]
                ]
              ]
            ]
          ]
        [ ask close-deers with [ momid = myid and cwd = 0 ] [
          if random-float 1 < a1314fof [                                   ;60-90 contacts;@#@#@#a1314fof   ;yearling female to own fawns
            set cwd 1
            set cum-trans cum-trans + 1                         ;#2May2020 for sa
            ]
          ]
          ]
        ]
      ]
    if aim = 15 and sex = 2 [                                              ;female yearling 15 mo age
      let myid who
      ask close-deers with [ who = lmom and cwd = 0 ] [                    ;mom
        ifelse (groid >= 0 and groid = lgroid)                             ;both in the same group
        [ if random-float 1 < a15fmom [                                    ;10-20 contacts if in the same group;@#@#@#a15fmom
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
          ]
        [ if random-float 1 < (a15fmom / 2) [                              ;0-10 contacts;@#@#@#a15fmom/2
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
      ]
      ask close-deers with [ momid = myid and cwd = 0 ] [
        if random-float 1 < a15fof [                                       ;60-90 contacts;@#@#@#a15fof
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
      ask close-deers with [ momid = lmom and aim = 15 and sex = 2 and cwd = 0 ] [ ;susceptible full sib females
        ifelse (groid >= 0 and groid = lgroid)                                     ;both in the same group
        [ if random-float 1 < a15ffsf [                                     ;10-20 contacts if in the same group;@#@#@#a15ffsf
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
          ]
        [ if random-float 1 < (a15ffsf / 2) [                               ;0-10 contacts;@#@#@#<a15ffsf/2
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
      ]
      if lgroid >= 0 [
        ask local-deers with [ groid = lgroid and who != lmom and momid != lmom and cwd = 0 ] [  ;susceptible group members, excluding mom and full sibs
          if aim = 15 or aim > 26 [
            if random-float 1 < a15fgrf [                                   ;5-10 contacts;@#@#@#a15fgrf
              set cwd 1
              set cum-trans cum-trans + 1                         ;#2May2020 for sa
              ]
            ]
          ]
        ]
      if lgroid >= 0 [                               ;if the infected female is a group member, it can also transmit the infection to one of the non-group female adult in the neighborhood
        if any? local-deers with [ sex = 2 and aim > 26 and groid != lgroid and who != lmom and cwd = 0 ] [
          ask n-of 1 local-deers with [ sex = 2 and aim > 26 and groid != lgroid and who != lmom and cwd = 0 ] [
            if random-float 1 < a15fngrf [                                 ;@#@#@#a15fngrf
              set cwd 1
              set cum-trans cum-trans + 1                         ;#2May2020 for sa
              ]
            ]
          ]
        ]
      if lgroid = -1 [                              ;if the infected female is solitary, it can also infect one of the adult in the neighborhood
        if any? local-deers with [ sex = 2 and aim > 26 and who != lmom and momid != lmom and cwd = 0 ] [
          ask n-of 1 local-deers with [ sex = 2 and aim > 26 and who != lmom and momid != lmom and cwd = 0 ] [
            if random-float 1 < a15fngrf [                                 ;@#@#@#a15fngrf
              set cwd 1
              set cum-trans cum-trans + 1                         ;#2May2020 for sa
              ]
            ]
          ]
        ]
      ]
    if aim > 15 and aim < 19 and sex = 2 [                                 ;female yearling 16-18 mo age
      let myid who
      ask close-deers with [ who = lmom and cwd = 0 ] [                    ;mom
        ifelse (groid >= 0 and groid = lgroid)                             ;both in the same group
        [ if random-float 1 < a1618fmom [                                  ;10-20 contacts if in the same group;@#@#@#a1618fmom
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
          ]
        [ if random-float 1 < (a1618fmom / 2) [                            ;0-10 contacts;@#@#@#a1618fmom/2
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
      ]
      ask close-deers with [ momid = myid and cwd = 0 ] [
        ifelse sex = 1
        [ if random-float 1 < a1618fofm [                                  ;10-20 contacts;@#@#@#a1618fofm
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
        [ if random-float 1 < a1618foff [                                  ;20-30 contacts;@#@#@#a1618foff
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
      ]
      ask close-deers with [ momid = lmom and aim > 15 and aim < 19 and sex = 2 and cwd = 0 ] [ ;susceptible full sib females
        ifelse (groid >= 0 and groid = lgroid)                             ;both in the same group
        [ if random-float 1 < a1618ffsf [                                  ;10-20 contacts if in the same group;@#@#@#a1618ffsf
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
          ]
        [ if random-float 1 < (a1618ffsf / 2) [                           ;0-10 contacts;@#@#@#a1618ffsf/2
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
      ]
      if lgroid >= 0 [
        ask local-deers with [ groid = lgroid and who != lmom and momid != lmom and cwd = 0 ] [  ;susceptible group members, excluding mom and full sibs
          if aim > 15 [
            if random-float 1 < a1618fgrf [                              ;5-10 contacts;@#@#@#a1618fgrf
              set cwd 1
              set cum-trans cum-trans + 1                         ;#2May2020 for sa
              ]
            ]
          ]
        ]
      if lgroid >= 0 [                               ;if the infected female is a group member, it can also transmit the infection to one of the non-group female adult in the neighborhood
        if any? local-deers with [ sex = 2 and aim > 27 and groid != lgroid and who != lmom and cwd = 0 ] [
          ask n-of 1 local-deers with [ sex = 2 and aim > 27 and groid != lgroid and who != lmom and cwd = 0 ] [
            if random-float 1 < a1618fngrf [                              ;@#@#@#a1618fngrf
              set cwd 1
              set cum-trans cum-trans + 1                         ;#2May2020 for sa
              ]
            ]
          ]
        ]
        if lgroid = -1 [                              ;if the infected female is solitary, it can also infect one of the adult in the neighborhood
        if any? local-deers with [ sex = 2 and aim > 27 and who != lmom and momid != lmom and cwd = 0 ] [
          ask n-of 1 local-deers with [ sex = 2 and aim > 27 and who != lmom and momid != lmom and cwd = 0 ] [
            if random-float 1 < a1618fngrf [                              ;@#@#@#a1618fngrf
              set cwd 1
              set cum-trans cum-trans + 1                         ;#2May2020 for sa
              ]
            ]
          ]
        ]
      ]
    if aim > 18 and aim < 21 and sex = 2 [                                 ;female yearling 19, 20 mo age
      let myid who
      ask close-deers with [ who = lmom and cwd = 0 ] [                    ;mom
        ifelse (groid >= 0 and groid = lgroid)                             ;both in the same group
        [ if random-float 1 < a1920fmom [                                  ;0-5 contacts if in the same group;@#@#@#a1920fmom
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
          ]
        [ if random-float 1 < a1920fmom [                                  ;0-5 contacts;#@#@#@@#@#@#a1920fmom
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
      ]
      ask close-deers with [ momid = myid and cwd = 0 ] [                  ;own fawns
        ifelse sex = 1
        [ if random-float 1 < a1920fofm [                                  ;0-5 contacts male fawn;@#@#@#a1920fofm
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
        [ if random-float 1 < a1920foff [                                  ;5-10 contacts;@#@#@#a1920foff
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
      ]
      ask close-deers with [ momid = lmom and aim > 18 and aim < 21 and sex = 2 and cwd = 0 ] [ ;susceptible full sib females
        ifelse (groid >= 0 and groid = lgroid)                             ;both in the same group
        [ if random-float 1 < a1920ffsf [                                  ;0-5 contacts if in the same group;@#@#@#a1920ffsf
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
          ]
        [ if random-float 1 < a1920ffsf [                                  ;0-5 contacts;@#@#@#a1920ffsf
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
      ]
      if lgroid >= 0 [
        ask local-deers with [ groid = lgroid and who != lmom and momid != lmom and cwd = 0 ] [  ;susceptible group members, excluding mom and full sibs
          if aim > 15 [
            if random-float 1 < a1920fgrf [                                ;5-10 contacts;@#@#@#a1920fgrf
              set cwd 1
              set cum-trans cum-trans + 1                         ;#2May2020 for sa
              ]
            ]
          ]
        ]
      if lgroid >= 0 [                               ;if the infected female is a group member, it can also transmit the infection to one of the non-group female adult in the neighborhood
        if any? local-deers with [ sex = 2 and aim > 27 and groid != lgroid and who != lmom and cwd = 0 ] [
          ask n-of 1 local-deers with [ sex = 2 and aim > 27 and groid != lgroid and who != lmom and cwd = 0 ] [
            if random-float 1 < a1920fngrf [                              ;5-8 contacts  ;@#@#@#a1920fngrf
              set cwd 1
              set cum-trans cum-trans + 1                         ;#2May2020 for sa
              ]
            ]
          ]
        ]
        if lgroid = -1 [                              ;if the infected female is solitary, it can also infect one of the adult in the neighborhood
        if any? local-deers with [ sex = 2 and aim > 27 and who != lmom and momid != lmom and cwd = 0 ] [
          ask n-of 1 local-deers with [ sex = 2 and aim > 27 and who != lmom and momid != lmom and cwd = 0 ] [
            if random-float 1 < a1920fngrf [                               ;5-8 contacts;@#@#@#a1920fngrf
              set cwd 1
              set cum-trans cum-trans + 1                         ;#2May2020 for sa
              ]
            ]
          ]
        ]
      ]
    if aim > 20 and aim < 25 and sex = 2 [                                 ;female yearling 21-24 mo age
      let myid who
      ask close-deers with [ who = lmom and cwd = 0 ] [                    ;mom
        ifelse (groid >= 0 and groid = lgroid)                             ;both in the same group
        [ if random-float 1 < a2124fmom [                                  ;10-20 contacts if in the same group;#$#$#$a2124fmom
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
          ]
        [ if random-float 1 < (a2124fmom / 2) [                            ;0-10 contacts;@#@#@#a2124fmom/2
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
      ]
      ask close-deers with [ momid = myid and cwd = 0 ] [
        ifelse sex = 1
        [ if random-float 1 < a2124fofm [                                 ;5-10 contacts;@#@#@#a2124fofm
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
        [ if random-float 1 < a2124foff [                                ;10-20 contacts;@#@#@#a2124foff
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
      ]
      ask close-deers with [ momid = lmom and aim > 20 and aim < 25 and sex = 2 and cwd = 0 ] [ ;susceptible full sib females
        ifelse (groid >= 0 and groid = lgroid)                           ;both in the same group
        [ if random-float 1 < a2124ffsf [                                ;10-20 contacts if in the same group;@#@#@#a2124ffsf
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
          ]
        [ if random-float 1 < (a2124ffsf / 2) [                          ;0-10 contacts;@#@#@#a2124ffsf/2
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
      ]
      if lgroid >= 0 [
        ask local-deers with [ groid = lgroid and who != lmom and momid != lmom and cwd = 0 ] [  ;susceptible group members, excluding mom and full sibs
          if aim > 20 [
            if random-float 1 < a2124fgrf [                              ;5-20 contacts;@#@#@#a2124fgrf
              set cwd 1
              set cum-trans cum-trans + 1                         ;#2May2020 for sa
              ]
            ]
          ]
        ]
      if lgroid >= 0 [                               ;if the infected female is a group member, it can also transmit the infection to one of the non-group female adult in the neighborhood
        if any? local-deers with [ sex = 2 and aim > 32 and groid != lgroid and who != lmom and cwd = 0 ] [
          ask n-of 1 local-deers with [ sex = 2 and aim > 32 and groid != lgroid and who != lmom and cwd = 0 ] [
            if random-float 1 < a2124fngrf [                             ;@#@#@#a2124fngrf
              set cwd 1
              set cum-trans cum-trans + 1                         ;#2May2020 for sa
              ]
            ]
          ]
      ]
      if lgroid = -1 [                              ;if the infected female is solitary, it can also infect one of the adult in the neighborhood
        if any? local-deers with [ sex = 2 and aim > 32 and who != lmom and momid != lmom and cwd = 0 ] [
          ask n-of 1 local-deers with [ sex = 2 and aim > 32 and who != lmom and momid != lmom and cwd = 0 ] [
            if random-float 1 < a2124fngrf [                         ;0-5 contacts;@#@#@#a2124fngrf
              set cwd 1
              set cum-trans cum-trans + 1                         ;#2May2020 for sa
              ]
            ]
          ]
        ]
      ]
    if (aim > 24 and aim < 27 and sex = 2) [                        ;female, young adult 25,26 mo
      let myid who
      ifelse not any? close-deers with [ momid = myid ]             ;female young adult without/with fawns
      [ ask close-deers with [ momid = lmom and sex = 2 and aim > 24 and aim < 27 and cwd = 0 ] [ ;susceptible full sibs females
        let myid1 who
        if not any? close-deers with [ momid = myid1 ] [            ;female full sib without fawns
          if random-float 1 < a2526ffsf [                           ;5-10 contacts;@#@#@#a2526ffsf
            set cwd 1
            set cum-trans cum-trans + 1                         ;#2May2020 for sa
            ]
          ]
        ]

        if lgroid >= 0 [
          ask local-deers with [ groid = lgroid and who != lmom and momid != lmom and cwd = 0 ] [  ;susceptible group members, excluding mom and full sibs
            if aim > 12 and aim < 15 [                                 ;female yearlings in group
              let myid2 who
              if not any? close-deers with [ momid = myid2 ] [         ;without fawns
                if random-float 1 < a2526fgrf [                        ;5-10 contacts;@#@#@#a2526fgrf
                  set cwd 1
                  set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
              ]
            ]
            if aim > 24 and aim < 27 [                                  ;young females in group age 25 and 26
              let myid1 who
              if not any? close-deers with [ momid = myid1 ] [          ;without fawns
                if random-float 1 < a2526fgrf [                         ;5-10 contacts   yearling and adult female gr memb;@#@#@#a2526fgrf
                  set cwd 1
                  set cum-trans cum-trans + 1                         ;#2May2020 for sa
                  ]
                ]
              ]
            ]
          let local-myearlings local-deers with [ sex = 1 and aim > 12 and aim < 15 ]    ;@#@#@#$$$local yearlings upto 5 to simulate temporary mixed group
          if any? local-myearlings [
            let n-my 0
            ifelse count local-myearlings > 5
            [ set n-my 5 ]
            [ set n-my count local-myearlings ]
            ask n-of random n-my local-myearlings [
              if random-float 1 < a2526fmy [                            ;5-10 contacts   local male yearling;@#@#@#a2526fmy
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
            ]
          ]
          ]
        ]
      [ ask close-deers with [ momid = myid and cwd = 0 ] [
        if random-float 1 < a2526fof [                                 ;@#@#@#a2526fof
          set cwd 1
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
        ]
      ]

    if aim > 26 and sex = 2 [                                          ;adult female > 26 mo age
      let myid who
      ask close-deers with [ who = lmom and cwd = 0 ] [                ;mom
        ifelse (groid >= 0 and groid = lgroid)                         ;both in the same group
        [ if d < 5 [
          if random-float 1 < a27fd4mom [                              ;0-10 contacts if in the same group;@#@#@#a27fd4mom
            set cwd 1
            set cum-trans cum-trans + 1                         ;#2May2020 for sa
            ]
          ]
          if d > 6 and d < 11 [
            if random-float 1 < a27fd10mom [                             ;0-10 contacts if in the same group;@#@#@#a27fd10mom
              set cwd 1
              set cum-trans cum-trans + 1                         ;#2May2020 for sa
            ]
          ]
          if d > 10 [
            if random-float 1 < a27fd12mom [                          ;@#@#@#a27fd12mom
              set cwd 1
              set cum-trans cum-trans + 1                         ;#2May2020 for sa
            ]
          ]
        ]
        [ if d < 5 or d > 6 [
          if random-float 1 < (a27fd4mom / 2) [                       ;0-5 contacts ;@#@#@#a27fd4mom/2
            set cwd 1
            set cum-trans cum-trans + 1                         ;#2May2020 for sa
            ]
          ]
          ]
      ]
      ask close-deers with [ momid = myid and cwd = 0 ] [
        ifelse sex = 1
        [ if d < 5 [
          if random-float 1 < a27fd4ofm [                            ;5-10 contacts;@#@#@#a27fd4ofm
            set cwd 1
            set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
          ]
          if d >= 5 and d <= 7 [                                    ;@#@#@#11Oct19=7changed t0 >=5 and <= 7
            if random-float 1 < a27fd7ofm [                         ;60-90 contacts;@#@#@#a27fd7ofm
              set cwd 1
              set cum-trans cum-trans + 1                         ;#2May2020 for sa
              ]
            ]
          if d > 7 and d < 11 [
            if random-float 1 < a27fd10ofm [                         ;10-20 contacts;@#@#@#a27fd10ofm
              set cwd 1
              set cum-trans cum-trans + 1                         ;#2May2020 for sa
            ]
          ]
          if d > 10 [
            if random-float 1 < a27fd12ofm [                        ;@#@#@#a27fd12ofm
              set cwd 1
              set cum-trans cum-trans + 1                         ;#2May2020 for sa
            ]
          ]
          ]
        [ if d < 5 [
          if random-float 1 < a27fd4off [                            ;10-20 contacts;@#@#@#a27fd4off
            set cwd 1
            set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
          ]
          if d >= 5 and d <= 7 [                                               ;@#@#@#11Oct19=7changed t0 >=5 and <= 7
            if random-float 1 < a27fd7off [                        ;60-90 contacts;@#@#@#a27fd7off
              set cwd 1
              set cum-trans cum-trans + 1                         ;#2May2020 for sa
              ]
            ]
          if d > 7 and d < 11 [
            if random-float 1 < a27fd10off [                       ;20-30 contacts;@#@#@#a27fd10off
              set cwd 1
              set cum-trans cum-trans + 1                         ;#2May2020 for sa
            ]
          ]
          if d > 10 [
            if random-float 1 < a27fd12off [                       ;@#@#@#a27fd12off
              set cwd 1
              set cum-trans cum-trans + 1                         ;#2May2020 for sa
              ]
            ]
          ]
        ]

      if lgroid >= 0 [
        ask local-deers with [ groid = lgroid and who != lmom and cwd = 0 ] [  ;susceptible group members, excluding mom
          if aim > 24 [
            if (d < 5) or (d > 6 and d < 11) [
              if random-float 1 < a27fd10grf [                      ;5-10 contacts;@#@#@#a27fd10grf    ;gestation and prerut
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
              ]
            if d > 10 [
              if random-float 1 < a27fd12grf [                      ;@#@#@#a27fd12grf
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
              ]
            ]
          ]
        ]
      if lgroid >= 0 [                               ;if the infected female is a group member, it can also transmit the infection to one of the non-group female adult in the neighborhood
        if any? local-deers with [ sex = 2 and aim > 32 and groid != lgroid and who != lmom and cwd = 0 ] [
          ask n-of 1 local-deers with [ sex = 2 and aim > 32 and groid != lgroid and who != lmom and cwd = 0 ] [
            if (d < 5) or (d > 6 and d < 11) [
              if random-float 1 < a27fd10ngrf [                     ;0-5 contacts;@#@#@#a27fd10ngrf
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
              ]
            if d > 10 [
              if random-float 1 < a27fd12ngrf [                    ;5-8 contacts;@#@#@#a27fd12ngrf
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
              ]
            ]
          ]
        ]
      if lgroid = -1 [                              ;if the infected female is solitary, it can also infect one of the adult in the neighborhood
        if any? local-deers with [ sex = 2 and aim > 32 and who != lmom and cwd = 0 ] [
          ask n-of 1 local-deers with [ sex = 2 and aim > 32 and who != lmom and cwd = 0 ] [
            if (d < 5) or (d > 6 and d < 11) [
              if random-float 1 < a27fd10ngrf [                  ;0-5 contacts@#@#@#a27fd10ngrf
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
              ]
            if d > 10 [
              if random-float 1 < a27fd12ngrf [                   ;5-8 contacts;@#@#@#a27fd12ngrf
                set cwd 1
                set cum-trans cum-trans + 1                         ;#2May2020 for sa
                ]
              ]
            ]
          ]
        ]
      ]
    if aim > 14 and sex = 1 and d < 11 [                             ;cwd+ is yearling or adult MALE and non rut season.     http://community.deergear.com/the-hunt/bachelor-group-behavior/
      ifelse (mgroid < 0)                                            ;not a bachelor group member
      [ let potgrm close-deers with [ mgroid = -1 and cwd = 0 ]      ;uninfected solitary males in the vicinity
        if any? potgrm [
          ask n-of 1 potgrm [                                        ;susceptible, solitary males in the neighborhood
            if (random-float 1 < a15msb) [                           ;1 to less than 5 encounters per month ;17Nov17 0.19-0.38 changed to 0-0.19; probability of infection 0.19-0.38
              set cwd 1                                              ;@@##@@##a15msb
              set cum-trans cum-trans + 1                         ;#2May2020 for sa
              ]
            ]
          ]
        ]
      [ ask local-deers with [ mgroid = [ mgroid ] of myself and cwd = 0 ] [ ;susceptible male deer in bachelor group ;@@##22Sep19 changed who= to mgroid =
        if (random-float 1 < a15mgrb) [                               ;5-15 contacts
          set cwd 1                                                   ;@#@#@#a15mgrb
          set cum-trans cum-trans + 1                         ;#2May2020 for sa
          ]
        ]
        ]
      ]
    ]
end
to deer-die-CWD
  let lgroid groid
  let lwho who
 ;----------------------------------------------------------------- fawns upto 6 months
  ifelse (aim < 6.5)
  [ ifelse (sex = 1)
    [ set counter1 0
      if (gr = -1) [
        review-group-dynamics
        ]
      die
      ]
    [ set counter1 0
      if (gr = -1) [
        review-group-dynamics
        ]
      die
      ]
    ]
 ;------------------------------------------------------------------ fawns 7 to 12 months
  [ ifelse (aim < 12.5)
    [ ifelse (sex = 1)
      [ set counter1 0
        if (gr = -1) [
          review-group-dynamics
          ]
        die
        ]
      [ set counter1 0
        if (gr = -1) [
          review-group-dynamics
          ]
        die
        ]
      ]
  ;------------------------------------------------------------------ yearlings 13 to 24 months
    [ ifelse (aim < 24.5)
      [ ifelse (sex = 1)
        [ set counter1 0
          if (gr = -1) [
            review-group-dynamics
            ]
          die
          ]
        [ ask deers in-radius-nowrap 1.5 with [ momid = lwho and aim < 2.5 ] [ ;my fawns
          set counter1 counter1 + 1
          die
          ]
          ifelse (gl > 0)
          [ new-group-leader ]
          [ if (gr = -1) [
            review-group-dynamics
            ]
            ]
          set counter1 0
          die
          ]
        ]
 ;----------------------------------------------------------------------- male 25 to 240 and more than 240
      [ ifelse (sex = 1)
        [ if (ml = 1) [
          let lmgroid mgroid
          ask deers in-radius-nowrap 3 with [ mgroid = lmgroid ] [ set mgroid -1 ]
          ]
          die
        ]
 ;--------------------------------------------------------------------- female 25 to 240
        [ ask deers in-radius-nowrap 1.5 with [ momid = lwho and aim < 2.5 ] [
          set counter1 counter1 + 1
          die
          ]
          ifelse (gl > 0)
          [ new-group-leader ]
          [ if (gr = -1) [ review-group-dynamics ]
            ]
          set counter1 0
          die
          ]
        ]
      ]
    ]
end
to deer-die                                               ;turtle procedure: non-hunting mortality
  let rn precision (random-float 1) 3
 ;------------------------------------------------------------- fawns upto 6 months
  ifelse (aim < 6.5) [
    let lmort mf6nhm
    if(sex = 2) [ set lmort ff6nhm ]
    if rn < lmort [
      if (gr = -1) [
        review-group-dynamics
        ]
      die
      ]
    ]
 ;---------------------------------------------------------------  7 to 12 months
  [ ifelse (aim < 12.5) [
    let lmort mf12nhm
    if (sex = 2) [ set lmort ff12nhm ]
    if rn < lmort [
      if (gr = -1) [ review-group-dynamics ]
      if (sex = 1 and mgroid > 0) [ review-bachelor-group ]
      die
      ]
    ]
 ;------------------------------------------------------------- 13 to 24
    [ ifelse (aim < 24.5) [
      ifelse (sex = 1) [
        if rn < mynhm [
          if (gr = -1) [ review-group-dynamics ]
          if (mgroid > 0) [ review-bachelor-group ]
          die
          ]
        ]
      [ if rn < fynhm [
        ask deers in-radius-nowrap 2 with [ momid = [ who ] of myself and aim < 2.5 ] [
          review-group-dynamics
          die
          ]
        ifelse (gl > 0) [ new-group-leader ] [ if (gr = -1) [ review-group-dynamics ] ]
        die
        ]
        ]
      ]
;--------------------------------------------------------------- male 25 to 240 and more than 240
      [ ifelse (sex = 1)
        [ let lmort 0.8
          if aim < 240 [ set lmort precision (manhm - oldm) 3 ]
          if rn < lmort [
            ifelse ml = 1
            [ask deers with [ mgroid = [ groid ] of myself ] [ set mgroid -1 ] ]
            [ review-bachelor-group ]
            die
            ]
          ]
;------------------------------------------------------------------- female 25 to 240
        [ let lmort 0.8
          if aim < 240 [ set lmort precision (fanhm - oldf) 3 ]
          if rn < lmort [
            ask deers in-radius-nowrap 2 with [ momid = [ who ] of myself and aim < 2.5 ] [
              review-group-dynamics
              die
              ]
            ifelse (gl > 0) [ new-group-leader ] [ if (gr = -1) [ review-group-dynamics ] ]
            die
          ]
        ]
      ]
    ]
  ]
end
to new-group-leader
  let lgroid groid
  let ngroid -1
  let my-group deers in-radius-nowrap 3 with [ groid = lgroid and who != lgroid ]
  let pot-groupleaders my-group with [ aim > 18 and sex = 2 ]
  ifelse (any? pot-groupleaders)
  [ ifelse (any? pot-groupleaders with [ aim > 29 ] )
    [ set ngroid ( [ who ] of one-of pot-groupleaders with [ aim > 29 ] ) ]
    [ set ngroid ( [ who ] of one-of pot-groupleaders) ] ]
  [ let other-groupleaders-here deers-here with [ gl = 1 and gr < 3 and who != lgroid ]
    if (any? other-groupleaders-here)
    [ set ngroid ( [ who ] of one-of other-groupleaders-here) ]
  ]
  ifelse ngroid = -1
  [ ask my-group [
    set gr -2
    set groid -1
    move-to min-one-of patches with [ do = 1 and count deers-here < 40 ] [ distance-nowrap myself ]          ;to prevent solitary females clusteriing on patches
    ]
    ]
  [ ask deer ngroid [
      set gl 1
      if gr < 0 [ set gr 0 ]
      set gr (gr + count my-group)
    ]
    ask my-group [ set groid ngroid ]
  ]
  set n_leaders_lost (n_leaders_lost + 1)
end

to review-bachelor-group
  if (d = 1 or d > 10) [ stop ]
  if is-turtle? deer mgroid [
    ask deer mgroid [
      set gr (gr - 1)
      if gr <= 1 [ set ml 0 ]
    ]
  ]
end

to form-bachelor-groups
  let tmbg ((my + ma) / 4)                                      ;appropriate # of group leaders, mean bchelor group size is 4
  let male-leaders deers with [ sex = 1 and ml = 1 ]
  let diff (tmbg - count male-leaders)
  if (diff > 0) [                                               ;designate new group leaders to fill up slots
    let new-leaders n-of diff deers with [ sex = 1 and aim > 32 and ml = 0 ]
    ask new-leaders [
      set ml 1
      set mgroid who
      ]
    set male-leaders (turtle-set male-leaders new-leaders)    ;add to leaders group
    ]
  ask male-leaders [                                         ;Fill out groups
    let lgroid mgroid
    let pot-groupmembers-here deers in-radius-nowrap 1.5 with [ sex = 1 and mgroid = -1 and ml = 0 ]
    let group-members pot-groupmembers-here with [ mgroid = lgroid ]
    let to-fill max (list 0 (round mean-bachelor-group-size - count group-members))   ;needed members
    let possible-fills count pot-groupmembers-here
    if to-fill > 0 [                                          ;fill to target if possible
      let new-members n-of (min list to-fill possible-fills) pot-groupmembers-here
      ask new-members [
        set mgroid lgroid
        ]
      set group-members (turtle-set group-members new-members)
      ]
    set gr count group-members
    if gr <= 1 [ set ml 0 ]
    ]
end

to deer-mdisperse                                                          ;turtle procedure: male yearling dispersal
  ask patch-here [
    if (dfp < 0.72) [
      let mdd (35.07 - (48.14 * (dfp)))                                    ;relationship between mean dispersal distance and proportion of forest cover to estimate the mean dispersal distance Long et al., 2005; Difenbach et al., 2008
      let sddd sqrt(e ^ (3.51 + (0.077 * mdd)))
      let lv ln (1 + (sddd ^ 2) / (mdd ^ 2))
      let lm ln mdd - (lv / 2)
      let ls sqrt lv
      set dd (exp (random-normal lm ls) * 0.6214)                           ;this is in miles, 1 patch fd is 1 mile, so we convert the dispersal distance to miles
      ]
    ]
  let counter-md 0
  rt random 360
  while [ counter-md < round dd ]
  [ ask patch-here [
      set dxp (pxcor)
      set dyp (pycor)
    ]
    fd 1
    set counter-md counter-md + 1
    ask patch-here [
      set dxn (pxcor)
      set dyn (pycor)
      ]
    if (abs(dxn - dxp) > 1 or abs(dyn - dyp) > 1)
    [ set ndd ndd + 1
      set cwd 0                                                             ;deer dispersing INTO the landscape are CWD-free.
      ;set cwdm 20 + random 6
      ;set cwdi 3 + random 2
      ;set cwdc 15 + random 4
      set cwdi 6 + random 5                       ;3 + random 2 ;duration after which a deer becomes infectious - from exposure to becoming infectious
      set cwdc 21 + random 5                      ; 15 + random 4;exposure to overt clinical signs in GG- derived from Johnson et al., 2011 duration after which a deer starts exhibiting clinical signs of CWD
      set cwdm cwdc + random 2
      set momid 0 ]
    ]
  finalize-home-patch
end

to deer-fdisperse
  let counter-fd 0
  rt random 360
  set dd round (random-normal mean-female-dispersal-distance stddev-dispersal-distance)
  while [ counter-fd < dd ]
  [ ask patch-here [
      set dxp (pxcor)
      set dyp (pycor) ]
    fd 1
    set counter-fd counter-fd + 1
    ask patch-here [
      set dxn (pxcor)
      set dyn (pycor)
      ]
    if (abs(dxn - dxp) > 1 or abs(dyn - dyp) > 1) [
      set ndd ndd + 1
      set cwd 0
      ;set cwdm 20 + random 6
      ;set cwdi 3 + random 2
      ;set cwdc 15 + random 4
      set cwdi 6 + random 5                       ;3 + random 2 ;duration after which a deer becomes infectious - from exposure to becoming infectious
      set cwdc 21 + random 5                      ; 15 + random 4;exposure to overt clinical signs in GG- derived from Johnson et al., 2011 duration after which a deer starts exhibiting clinical signs of CWD
      set cwdm cwdc + random 2
      set momid 0
      ]
    ]
  finalize-home-patch
end

to finalize-home-patch
  let fhp 0
  ask patch-here [
    if (do != 1)[
      set fhp 1
    ]
    ]
  if (fhp > 0)[
    set dxp (pxcor)
    set dyp (pycor)
    move-to min-one-of patches with [ do = 1 and count deers-here < 40 ] [ distance myself ]
    set color blue
    ask patch-here [
      set dxn (pxcor)
      set dyn (pycor)
      ]
    if (abs(dxn - dxp) > 1 or abs(dyn - dyp) > 1)                                                  ; dispersing deer goes out of the model landscape
    [ set ndd ndd + 1
      set cwd 0
      ;set cwdm 20 + random 6
      ;set cwdi 3 + random 2
      ;set cwdc 15 + random 4
      set cwdi 6 + random 5                       ;3 + random 2 ;duration after which a deer becomes infectious - from exposure to becoming infectious
      set cwdc 21 + random 5                      ; 15 + random 4;exposure to overt clinical signs in GG- derived from Johnson et al., 2011 duration after which a deer starts exhibiting clinical signs of CWD
      set cwdm cwdc + random 2
      set momid 0
      ]
    ]
end

to fawning
  let grfis 0
  if (aim = 13) [                                     ;turtle procedure: fawn female breeding
    if (random 100 <= juvenile-pregnancy-rate) [
      set grfis 0
      set tgroid groid
      ifelse (tgroid >= 0 and is-turtle? deer tgroid)               ;doe social group size adjustment in response to fawning
      [ ask deer tgroid [
        ifelse (gr > (doe-group-size-regulator - 2))
        [ set gr (gr - 1)
          set tgroid -1
          set tgr -2
          set grfis 1
        ]
        [ set gr (gr + 1)
          set tgroid groid
          set tgr -1
          ]
        ]
        if (grfis = 1) [
          set groid -1
          set gr -2
          move-to min-one-of patches with [ do = 1 and count deers-here < 40 ] [ distance-nowrap myself ] ;to prevent clustering of solitary females on a patch
          ]
        ]
      [ set tgroid -1
        set tgr -2
        ]
      deer-reproduce
      ]
    ]
  if (aim > 24) [                                         ;turtle procedure: adult doe breeding
    if (random 100 <= adult-pregnancy-rate) [
      set grfis 0
      set tgroid groid
      ifelse (gl > 0)                                     ;turtle procedure: new group formation after breeding
      [ set gr (gr + 2)
        if (gr > doe-group-size-regulator) [              ;turtle procedure: new group formation as a response to doe breeding
          let xgr gr - doe-group-size-regulator
          let group-members deers in-radius-nowrap 3 with [ groid = tgroid and sex = 2 and aim > 13 and gl = 0 ]      ;@@@test how many group members are outside this range
          let ngr 0
          ifelse (count group-members >= xgr)
          [ set ngr xgr ]
          [ set ngr count group-members ]
          ask n-of ngr group-members [
            new-group-formation
            ]
          ]
        set tgr -1
        ]
      [ ifelse (groid < 0 or not is-turtle? deer groid)
        [ ifelse (aim > 36 and n_leaders_lost > 0)         ;to prevent declining numbers of females with (gl = 1 )
          [ set gl 1
            set gr 2
            set groid who
            set tgroid who
            set tgr -1
            set n_leaders_lost (n_leaders_lost - 1)
            move-to min-one-of patches with [ do = 1 and count deers-here < 40 ] [ distance-nowrap myself ]          ;to prevent solitary females clusteriing on patches
            ]
          [ set tgroid -1
            set tgr -2
            ]
          ]
        [ ask deer tgroid [
          ifelse (gr > 4)
          [ set gr (gr - 1)
            set tgroid -1
            set tgr -2
            set grfis 1
            ]
          [ set gr (gr + 2)
            set tgroid groid
            set tgr -1
            ]
          ]
          if (grfis = 1) [
            set groid -1
            set gr -2
            move-to min-one-of patches with [ do = 1 and count deers-here < 40 ] [ distance-nowrap myself ]          ;to prevent solitary females clusteriing on patches
            ]
          ]
        ]
      deer-reproduce
      ]
    ]
  if (gl = 1 and gr < 4) [                                    ;turtle procedure: adjustment of doe social group size after the breeding season
    set tgroid groid
    let solitary-adult-females-here deers in-radius-nowrap 1.5 with [ sex = 2 and gr = -2 and aim >= 13 and cwdpr < cwdc ]
    let sd count solitary-adult-females-here
    let sd1 0
    set tgr 0
    while [ sd > 0 and sd1 < 3 ]
    [ ask n-of 1 solitary-adult-females-here [
      let tmomid who
      set groid tgroid
      set gr -1
      set tgr (tgr + 1)
      let my-fawns deers in-radius-nowrap 1.5 with [ momid = tmomid and aim = 1 ]
      if any? my-fawns [
        ask my-fawns [
          set groid tgroid
          set gr -1
          set tgr (tgr + 1)
          ]
        ]
      ]
      set sd1 sd1 + 1
      ]
    set gr (gr + tgr)
    if (gr = 0) [
      set gl 0
      set groid -1
      set gr -2
      ]
    ]
end

to deer-reproduce
  set mom who
  ifelse (aim < 13.5)
  [ hatch-deers 1 [
    set aim 1
    set gl 0
    set cwd 0
    ;set cwdm 20 + random 6
    ;set cwdi 3 + random 2
    ;set cwdc 15 + random 4
    set cwdi 6 + random 5                       ;3 + random 2 ;duration after which a deer becomes infectious - from exposure to becoming infectious
    set cwdc 21 + random 5                      ; 15 + random 4;exposure to overt clinical signs in GG- derived from Johnson et al., 2011 duration after which a deer starts exhibiting clinical signs of CWD
    set cwdm cwdc + random 2
    set momid mom
    set groid tgroid
    set gr tgr
    ifelse random 100 < 51
    [ set sex 1
      set mgroid -2 ]
    [ set sex 2 ]
    ]
    ]
  [ hatch-deers 2 [
    set aim 1
    set gl 0
    set cwd 0
    ;set cwdm 20 + random 6
    ;set cwdi 3 + random 2
    ;set cwdc 15 + random 4
    set cwdi 6 + random 5                       ;3 + random 2 ;duration after which a deer becomes infectious - from exposure to becoming infectious
    set cwdc 21 + random 5                      ; 15 + random 4;exposure to overt clinical signs in GG- derived from Johnson et al., 2011 duration after which a deer starts exhibiting clinical signs of CWD
    set cwdm cwdc + random 2
    set momid mom
    set groid tgroid
    set gr tgr
    ifelse random 100 < 51
    [ set sex 1
      set mgroid -2 ]
    [ set sex 2 ]
    ]
    ]
end

to new-group-formation
  move-to min-one-of patches with [ do = 1 and count deers-here < 40 ] [ distance-nowrap myself ]   ;to prevent solitary females clusteriing on patches
  let new-location patch-here
  let lgroid groid
  set groid -1
  set gr -2
  let lgr 0
  let my-fawns deers in-radius-nowrap 2 with [ momid = [ who ] of myself and aim = 1 ]
  if any? my-fawns [
    ask my-fawns [
      set gr -2
      set groid -1
      set lgr (lgr + 1)
      move-to new-location
      ]
  ]
  ask deer lgroid [
    set gr (gr - (lgr + 1))
    ]
end

to distance-to-farthest-spark
  let ini-cwd-focus patches with [ pcolor = orange ]             ;calculate distance to the farthest spark
  if any? ini-cwd-focus [
    ask n-of 1 patches with [ pcolor = orange ] [
      let inf-patches patches with [ pcolor = yellow ]
      let farthest-spark max-one-of inf-patches [ distance-nowrap myself ]
      ifelse is-patch? farthest-spark
      [ set max-dist-spark precision distance-nowrap farthest-spark 1 ]
      [ set max-dist-spark 0 ]
      ]
    ]
end

to hunting-mortality
  if (aim < 10) [
    ifelse (sex = 1)
    [ if (random-float 1 < mf12hm) [
      set tgroid groid
      hunting-mortality-mf12
      ]
      ]
    [ if (random-float 1 < ff12hm) [
      set tgroid groid
      set twho who
      hunting-mortality-ff12
      ]
      ]
    ]
  if (aim = 20) [
    ifelse (sex = 1)
    [ if (random-float 1 < myhm) [
      set tgroid groid
      hunting-mortality-my
      ]
      ]
    [ if (random-float 1 < fyhm) [
      set tgroid groid
      set twho who
      hunting-mortality-fy
      ]
      ]
    ]
  if (aim > 30) [
    ifelse (sex = 1)
    [ if (random-float 1 < mahm) [
      set tgroid groid
      hunting-mortality-ma
      ]
      ]
    [ if (random-float 1 < fahm) [
      set tgroid groid
      set twho who
      hunting-mortality-fa
      ]
      ]
    ]
end

to hunting-mortality-mf12
  if (gr = -1) [
    set counter1 0
    review-group-dynamics
    ]
  ifelse (random-float 1 <= %fawn-male-harvest-tested)
  [ if (cwd = 1) [
    set dcwdmf dcwdmf + 1
    set tcwdmf tcwdmf + 1
    ask patch-here [                         ;@#@#@#
      set cwd-d 1
      ]
    ]
    set tfamh tfamh + 1
    set tfmt tfmt + 1
    die
  ]
  [ if (cwd = 1) [
    set dcwdmf (dcwdmf + 1)
    ]
    set tfamh (tfamh + 1)
    die
    ]
end

to hunting-mortality-ff12
  if (gr = -1) [
    set counter1 0
    review-group-dynamics
    ]
  ifelse (random-float 1 <= %fawn-female-harvest-tested)
  [ if (cwd = 1) [
    set dcwdff dcwdff + 1
    set tcwdff tcwdff + 1
    ask patch-here [                         ;@#@#@#
      set cwd-d 1
      ]
    ]
    set tfafh tfafh + 1
    set tfft tfft + 1
    die
  ]
  [if (cwd = 1) [
    set dcwdff (dcwdff + 1)
    ]
    set tfafh (tfafh + 1)
    die
  ]
end

to hunting-mortality-my
  ifelse (random-float 1 <= %yearling-male-harvest-tested)
  [ if (cwd = 1) [
    set dcwdmy dcwdmy + 1
    set tcwdmy tcwdmy + 1
    ask patch-here [                         ;@#@#@#
      set cwd-d 1
      ]
    ]
    set tymh tymh + 1
    set tymt tymt + 1
    die
    ]
  [ if (cwd = 1) [
    set dcwdmy (dcwdmy + 1)
    ]
    set tymh (tymh + 1)
    die
    ]
end

to hunting-mortality-fy
  ifelse (gl > 0)
  [ new-group-leader ]
  [ if (gr = -1) [
    set counter1 0
    review-group-dynamics
    ]
  ]
  ifelse (random-float 1 <= %yearling-female-harvest-tested)
  [ if (cwd = 1) [
    set dcwdfy (dcwdfy + 1)
    set tcwdfy (tcwdfy + 1)
    ask patch-here [                         ;@#@#@#
      set cwd-d 1
      ]
    ]
    set tyfh tyfh + 1
    set tyft tyft + 1
    die
    ]
  [ if (cwd = 1) [
    set dcwdfy (dcwdfy + 1)
    ]
    set tyfh (tyfh + 1)
    die
    ]
end

to hunting-mortality-ma
  if (ml = 1) [
    let lmgroid mgroid
    ask deers in-radius-nowrap 3 with [ mgroid = lmgroid and ml = 0 ] [
      set mgroid -1
      ]
    ]
  ifelse (random-float 1 <= %adult-male-harvest-tested)
  [ if (cwd = 1) [
    set dcwdm (dcwdm + 1)
    set tcwdm (tcwdm + 1)
    ask patch-here [                         ;@#@#@#
      set cwd-d 1
      ]
    ]
    set tamh tamh + 1
    set tamt tamt + 1
    die
    ]
  [ if (cwd = 1) [
    set dcwdm (dcwdm + 1)
    ]
    set tamh tamh + 1
    die
    ]
end

to hunting-mortality-fa
  if (gl = 1) [ new-group-leader ]
  if (gr = -1) [
    set counter1 0
    review-group-dynamics
    ]
  ifelse (random-float 1 <= %adult-female-harvest-tested)
  [ if (cwd = 1) [
    set dcwdf (dcwdf + 1)
    set tcwdf (tcwdf + 1)
    ask patch-here [                         ;@#@#@#
      set cwd-d 1
      ]
    ]
    set tafh tafh + 1
    set taft taft + 1
    die
    ]
  [ if (cwd = 1) [
    set dcwdf (dcwdf + 1)
    ]
    set tafh tafh + 1
    die
    ]
end

;;;
to post-harvest-targeted-culling                                   ;@#@#@#
  ask patches with [ cwd-d = 1 ] [
    let test-cand deers in-radius-nowrap 1
    ask test-cand [
      if random-float 1 <= targeted-culling-percentage [
        set tgroid groid
        if (aim < 10) [
          ifelse (sex = 1)
          [ targeted-culling-mf12 ]
          [ targeted-culling-ff12 ]
          ]
        if (aim = 20) [
          ifelse (sex = 1)
          [ targeted-culling-my ]
          [ targeted-culling-fy ]
          ]
        if (aim > 30) [
          ifelse (sex = 1)
          [ targeted-culling-ma ]
          [ targeted-culling-fa ]
        ]
        ]
      ]
    ]
end

to targeted-culling-mf12
  if (gr = -1) [
    set counter1 0
    review-group-dynamics
    ]
  if (cwd = 1) [
    set tc-cwdmf tc-cwdmf + 1
    ;set tcwdmf tcwdmf + 1
    ]
  ;set tfamh tfamh + 1
  ;set tfmt tfmt + 1
  set tc-tmf tc-tmf + 1
  die
end

to targeted-culling-ff12
  if (gr = -1) [
    set counter1 0
    review-group-dynamics
    ]
  if (cwd = 1) [
    set tc-cwdff tc-cwdff + 1
    ;set tcwdff tcwdff + 1
    ]
  ;set tfafh tfafh + 1
  ;set tfft tfft + 1
  set tc-tff tc-tff + 1
  die
end

to targeted-culling-my
  if (cwd = 1) [
    set tc-cwdmy tc-cwdmy + 1
    ;set tcwdmy tcwdmy + 1
    ]
  ;set tymh tymh + 1
  ;set tymt tymt + 1
  set tc-tmy tc-tmy + 1
  die
end

to targeted-culling-fy
  ifelse (gl > 0)
  [ new-group-leader ]
  [ if (gr = -1) [
    set counter1 0
    review-group-dynamics
    ]
    ]
  if (cwd = 1) [
    set tc-cwdfy tc-cwdfy + 1
    ;set tcwdfy (tcwdfy + 1)
    ]
  ;set tyfh tyfh + 1
  ;set tyft tyft + 1
  set tc-tfy tc-tfy + 1
  die
end

to targeted-culling-ma
  if (ml = 1) [
    let lmgroid mgroid
    ask deers in-radius-nowrap 3 with [ mgroid = lmgroid and ml = 0 ] [
      set mgroid -1
      ]
    ]
  if (cwd = 1) [
    set tc-cwdm tc-cwdm + 1
    ;set tcwdm (tcwdm + 1)
    ]
  ;set tamh tamh + 1
  ;set tamt tamt + 1
  set tc-tm tc-tm + 1
  die
end

to targeted-culling-fa
  if (gl = 1) [ new-group-leader ]
  if (gr = -1) [
    set counter1 0
    review-group-dynamics
    ]
  if (cwd = 1) [
    set tc-cwdf tc-cwdf + 1
    ;set tcwdf (tcwdf + 1)
    ]
  ;set tafh tafh + 1
  ;set taft taft + 1
  set tc-tf tc-tf + 1
  die
end

;;;
to deer-mating
  let ter 0
  ifelse (aim > 30)
  [ set ter 1.5                      ;males > 30 months of age
    set nm (1 + random 6)
    set anm 0
    ]
  [ set ter 2.5                      ;males 18-30 months of age
    set nm (1 + random 3)
    set anm 0
    ]
  let female-deer-near-me deers in-radius-nowrap ter with [ sex = 2 and anm < nm and cwdpr < cwdc ]
  let pmates (count female-deer-near-me)
  if (pmates > 0) [
    let lcwd cwd
    let avmates min list nm pmates     ;available mates; minimum value from a list of 1) num of matings and 2) potential mates in the vicinity
    set anm avmates
    let exposures 0
    ask n-of avmates female-deer-near-me [
      if (cwd = 1 and lcwd = 0) [
        set exposures (exposures + 1)     ;male is exposed once
        ]
      if (cwd = 0 and lcwd = 1) [        ;chance for female to be infected
        let tpm (5 * transmission-prob)                ;IMP CHANGE TO BE UPDATED IN CURRENT VERSION; ASSUMING 5 CONTACTS PER MATING
        if (random-float 1 < tpm) [                   ;0.06 CHANGED TO TPM
         set cwd 1
         set cum-trans cum-trans + 1                         ;#2May2020 for sa]    ;4-5 contacts
        ]
        ]
      set anm anm + 1                                       ;@#@#@#4Dec19 IMP ADDITION***********************

      ]
    if cwd = 0 and exposures > 0 [
      while [ cwd = 0 and exposures > 0 ]
      [ let tpm-min (5 * transmission-prob)            ;IMP CHANGE TO BE UPDATED IN CURRENT VERSION; ASSUMING MIN 5 AND MAX 10 CONTACTS PER MATING
        let tpm-max (10 * transmission-prob)
        if random-float 1 < (tpm-min + random-float (tpm-max - tpm-min)) [     ;5-10 contacts; CHANGED FROM (0.06 + random-float 0.08) to tpm-min + random-float (tpm-max - tpm-min)
        set cwd 1
        set cum-trans cum-trans + 1                         ;#2May2020 for sa
        ]
        set exposures exposures - 1
        ]
      ]
    ]
end
;turtle procedure: doe social group leader loses leadership status if no group members left
to review-group-dynamics
  if (groid = -1 or gr = -2) [ stop ]
  if not is-turtle? deer groid [
    set groid -1
    set gr -2
    stop
  ]
  ask deer groid [
    set gr (gr - 1)
    if (gr <= 0) [
      set gl 0
      set groid -1
      set gr -2
      set n_leaders_lost (n_leaders_lost + 1)
      ]
    ]
end

;-------------REPORTERS---------------------------------
to-report doe-group-size-regulator
  report 6
end
to-report juvenile-pregnancy-rate
  report 20
end
to-report adult-pregnancy-rate
  report 80
end
to-report yearling-male-dispersal-rate           ;dispersal rates for yearling males range between 46% to 80% (Long et al., 2005)
  report 0.46
end
to-report yearling-female-dispersal-rate
  report 0.22
end
to-report mean-female-dispersal-distance
 report 11
end
to-report stddev-dispersal-distance
 report 4
end
to-report mean-bachelor-group-size
 report (1 + random-normal 4 1)
end
to-report cwd_area
  report count patches with [ pcolor = yellow and count deers-here with [ cwd = 1 ] > 0 ]
end
to-report mf
  report count deers with [ sex = 1 and aim < 12.5 ]
end
to-report my
  report count deers with [ sex = 1 and aim > 12.5 and aim < 25 ]
end
to-report ma
  report count deers with [ sex = 1 and aim > 24 ]
end
to-report ff
  report count deers with [ sex = 2 and aim < 12.5 ]
end
to-report fy
  report count deers with [ sex = 2 and aim > 12.5 and aim < 25 ]
end
to-report fa
  report count deers with [ sex = 2 and aim > 24 ]
end
to-report transmission-prob
  report 0.0128                                       ;SA1 0.0128 CHANGED TO 0.01216/ changed to 0.01344
end
to-report mcwd
  let cwd-deers deers with [ cwd = 1 ]
  report count cwd-deers with [ sex = 1 and aim > 24 ]
end
to-report mycwd
  let cwd-deers deers with [ cwd = 1 ]
  report count cwd-deers with [ sex = 1 and aim <= 24 and aim > 12 ]
end
to-report mfcwd
  let cwd-deers deers with [ cwd = 1 ]
  report count cwd-deers with [ sex = 1 and aim <= 12 ]
end
to-report fcwd
  let cwd-deers deers with [ cwd = 1 ]
  report count cwd-deers with [ sex = 2 and aim > 24 ]
end
to-report fycwd
  let cwd-deers deers with [ cwd = 1 ]
  report count cwd-deers with [ sex = 2 and aim <= 24 and aim > 12 ]
end
to-report ffcwd
  let cwd-deers deers with [ cwd = 1 ]
  report count cwd-deers with [ sex = 2 and aim <= 12 ]
end
to-report totcwdd
  let cwd-deers deers with [ cwd = 1 ]
  report count cwd-deers
end
;------------------------------------------------------------
to-report a3mom                               ;transmission prob from fawn age 3 mo to mom
  let a3 item 1 contactmatrix ;0
  let mommin item 1 a3 ;0
  let mommax item 2 a3 ;1
  let diffmaxmin ((mommax - mommin) * transmission-prob)
  set mommin mommin * transmission-prob
  report (mommin + random-float diffmaxmin)
end
to-report a3fs                               ;transmission prob from fawn age 3 mo to fullsibs
  let a3 item 1 contactmatrix ;0
  let fsmin item 7 a3 ;6
  let fsmax item 8 a3 ;7
  let diffmaxmin ((fsmax - fsmin) * transmission-prob)
  set fsmin fsmin * transmission-prob
  report (fsmin + random-float diffmaxmin)
end
to-report a3grf                               ;transmission prob from fawn age 3 mo to group females
  let a3 item 1 contactmatrix ;0
  let gfmin item 15 a3 ;14
  let gfmax item 16 a3 ;15
  let diffmaxmin ((gfmax - gfmin) * transmission-prob)
  set gfmin gfmin * transmission-prob
  report (gfmin + random-float diffmaxmin)
end
to-report a456mmom                            ;transmission prob from male fawn age 456 mo to mom
  let a456 item 2 contactmatrix ;1
  let mommin item 1 a456 ;0
  let mommax item 2 a456 ;1
  let diffmaxmin ((mommax - mommin) * transmission-prob)
  set mommin mommin * transmission-prob
  report (mommin + random-float diffmaxmin)
end
to-report a456mfs                            ;transmission prob from male fawn age 456 mo to fullsibs
  let a456 item 2 contactmatrix ;1
  let fsmin item 7 a456  ;6
  let fsmax item 8 a456  ;7
  let diffmaxmin ((fsmax - fsmin) * transmission-prob)
  set fsmin fsmin * transmission-prob
  report (fsmin + random-float diffmaxmin)
end
to-report a456mnsmf                            ;transmission prob from male fawn age 456 mo to nonsib male fawns
  let a456 item 2 contactmatrix ;1
  let nsmfmin item 11 a456   ;10
  let nsmfmax item 12 a456   ;11
  let diffmaxmin ((nsmfmax - nsmfmin) * transmission-prob)
  set nsmfmin nsmfmin * transmission-prob
  report (nsmfmin + random-float diffmaxmin)
end
to-report a456mnsff                            ;transmission prob from male fawn age 456 mo to nonsib female fawns
  let a456 item 2 contactmatrix ;1
  let nsffmin item 13 a456  ;12
  let nsffmax item 14 a456  ;13
  let diffmaxmin ((nsffmax - nsffmin) * transmission-prob)
  set nsffmin nsffmin * transmission-prob
  report (nsffmin + random-float diffmaxmin)
end
to-report a456mgrf                            ;transmission prob from male fawn age 456 mo to other females in group
  let a456 item 2 contactmatrix ;1
  let grfmin item 15 a456  ;14
  let grfmax item 16 a456  ;15
  let diffmaxmin ((grfmax - grfmin) * transmission-prob)
  set grfmin grfmin * transmission-prob
  report (grfmin + random-float diffmaxmin)
end
to-report a456fmom                            ;transmission prob from female fawn age 456 mo to mom
  let a456 item 3 contactmatrix ;2
  let mommin item 1 a456 ;0
  let mommax item 2 a456 ;1
  let diffmaxmin ((mommax - mommin) * transmission-prob)
  set mommin mommin * transmission-prob
  report (mommin + random-float diffmaxmin)
end
to-report a456ffs                            ;transmission prob from female fawn age 456 mo to fullsibs
  let a456 item 3 contactmatrix ;2
  let fsmin item 7 a456 ;6
  let fsmax item 8 a456 ;7
  let diffmaxmin ((fsmax - fsmin) * transmission-prob)
  set fsmin fsmin * transmission-prob
  report (fsmin + random-float diffmaxmin)
end
to-report a456fnsmf                            ;transmission prob from female fawn age 456 mo to nonsib male fawns
  let a456 item 3 contactmatrix ;2
  let nsmfmin item 11 a456 ;10
  let nsmfmax item 12 a456 ;11
  let diffmaxmin ((nsmfmax - nsmfmin) * transmission-prob)
  set nsmfmin nsmfmin * transmission-prob
  report (nsmfmin + random-float diffmaxmin)
end
to-report a456fnsff                            ;transmission prob from female fawn age 456 mo to nonsib female fawns
  let a456 item 3 contactmatrix ;2
  let nsffmin item 13 a456
  let nsffmax item 14 a456
  let diffmaxmin ((nsffmax - nsffmin) * transmission-prob)
  set nsffmin nsffmin * transmission-prob
  report (nsffmin + random-float diffmaxmin)
end
to-report a456fgrf                            ;transmission prob from female fawn age 456 mo to other females in group
  let a456 item 3 contactmatrix ;2
  let grfmin item 15 a456
  let grfmax item 16 a456
  let diffmaxmin ((grfmax - grfmin) * transmission-prob)
  set grfmin grfmin * transmission-prob
  report (grfmin + random-float diffmaxmin)
end
to-report a78mmom                            ;transmission prob from male fawn age 7/8 mo to mom
  let a78 item 4 contactmatrix ;3
  let mommin item 1 a78
  let mommax item 2 a78
  let diffmaxmin ((mommax - mommin) * transmission-prob)
  set mommin mommin * transmission-prob
  report (mommin + random-float diffmaxmin)
end
to-report a78mfsm                            ;transmission prob from male fawn age 7/8 mo to male fullsibs
  let a78 item 4 contactmatrix ;3
  let fsmmin item 7 a78
  let fsmmax item 8 a78
  let diffmaxmin ((fsmmax - fsmmin) * transmission-prob)
  set fsmmin fsmmin * transmission-prob
  report (fsmmin + random-float diffmaxmin)
end
to-report a78mfsf                            ;transmission prob from male fawn age 7/8 mo to female fullsibs
  let a78 item 4 contactmatrix ;3
  let fsfmin item 9 a78
  let fsfmax item 10 a78
  let diffmaxmin ((fsfmax - fsfmin) * transmission-prob)
  set fsfmin fsfmin * transmission-prob
  report (fsfmin + random-float diffmaxmin)
end
to-report a78mnsmf                            ;transmission prob from male fawn age 7/8 mo to nonsib male fawns
  let a78 item 4 contactmatrix ;3
  let nsmfmin item 11 a78
  let nsmfmax item 12 a78
  let diffmaxmin ((nsmfmax - nsmfmin) * transmission-prob)
  set nsmfmin nsmfmin * transmission-prob
  report (nsmfmin + random-float diffmaxmin)
end
to-report a78mnsff                            ;transmission prob from male fawn age 7/8 mo to nonsib female fawns
  let a78 item 4 contactmatrix ;3
  let nsffmin item 13 a78
  let nsffmax item 14 a78
  let diffmaxmin ((nsffmax - nsffmin) * transmission-prob)
  set nsffmin nsffmin * transmission-prob
  report (nsffmin + random-float diffmaxmin)
end
to-report a78mgrf                            ;transmission prob from male fawn age 7/8 mo to other females in group
  let a78 item 4 contactmatrix ;3
  let grfmin item 15 a78
  let grfmax item 16 a78
  let diffmaxmin ((grfmax - grfmin) * transmission-prob)
  set grfmin grfmin * transmission-prob
  report (grfmin + random-float diffmaxmin)
end
to-report a78fmom                            ;transmission prob from female fawn age 7/8 mo to mom
  let a78f item 5 contactmatrix ;4
  let mommin item 1 a78f
  let mommax item 2 a78f
  let diffmaxmin ((mommax - mommin) * transmission-prob)
  set mommin mommin * transmission-prob
  report (mommin + random-float diffmaxmin)
end
to-report a78ffsm                            ;transmission prob from female fawn age 7/8 mo to male fullsibs
  let a78f item 5 contactmatrix ;4
  let fsmmin item 7 a78f
  let fsmmax item 8 a78f
  let diffmaxmin ((fsmmax - fsmmin) * transmission-prob)
  set fsmmin fsmmin * transmission-prob
  report (fsmmin + random-float diffmaxmin)
end
to-report a78ffsf                            ;transmission prob from female fawn age 7/8 mo to female fullsibs
  let a78f item 5 contactmatrix ;4
  let fsfmin item 9 a78f
  let fsfmax item 10 a78f
  let diffmaxmin ((fsfmax - fsfmin) * transmission-prob)
  set fsfmin fsfmin * transmission-prob
  report (fsfmin + random-float diffmaxmin)
end
to-report a78fnsmf                            ;transmission prob from female fawn age 7/8 mo to nonsib male fawns
  let a78f item 5 contactmatrix ;4
  let nsmfmin item 11 a78f
  let nsmfmax item 12 a78f
  let diffmaxmin ((nsmfmax - nsmfmin) * transmission-prob)
  set nsmfmin nsmfmin * transmission-prob
  report (nsmfmin + random-float diffmaxmin)
end
to-report a78fnsff                            ;transmission prob from female fawn age 7/8 mo to nonsib female fawns
  let a78f item 5 contactmatrix ;4
  let nsffmin item 13 a78f
  let nsffmax item 14 a78f
  let diffmaxmin ((nsffmax - nsffmin) * transmission-prob)
  set nsffmin nsffmin * transmission-prob
  report (nsffmin + random-float diffmaxmin)
end
to-report a78fgrf                            ;transmission prob from female fawn age 7/8 mo to other females in group
  let a78f item 5 contactmatrix ;4
  let grfmin item 15 a78f
  let grfmax item 16 a78f
  let diffmaxmin ((grfmax - grfmin) * transmission-prob)
  set grfmin grfmin * transmission-prob
  report (grfmin + random-float diffmaxmin)
end
to-report a912mmom                            ;transmission prob from male fawn age 9/12 mo to mom
  let a912m item 6 contactmatrix ;5
  let mommin item 1 a912m
  let mommax item 2 a912m
  let diffmaxmin ((mommax - mommin) * transmission-prob)
  set mommin mommin * transmission-prob
  report (mommin + random-float diffmaxmin)
end
to-report a912mfsm                            ;transmission prob from male fawn age 9/12 mo to male fullsibs
  let a912m item 6 contactmatrix ;5
  let fsmmin item 7 a912m
  let fsmmax item 8 a912m
  let diffmaxmin ((fsmmax - fsmmin) * transmission-prob)
  set fsmmin fsmmin * transmission-prob
  report (fsmmin + random-float diffmaxmin)
end
to-report a912mfsf                            ;transmission prob from male fawn age 9/12 mo to female fullsibs
  let a912m item 6 contactmatrix ;5
  let fsfmin item 9 a912m
  let fsfmax item 10 a912m
  let diffmaxmin ((fsfmax - fsfmin) * transmission-prob)
  set fsfmin fsfmin * transmission-prob
  report (fsfmin + random-float diffmaxmin)
end
to-report a912mnsmf                            ;transmission prob from male fawn age 9/12 mo to nonsib male fawns
  let a912m item 6 contactmatrix  ;5
  let nsmfmin item 11 a912m
  let nsmfmax item 12 a912m
  let diffmaxmin ((nsmfmax - nsmfmin) * transmission-prob)
  set nsmfmin nsmfmin * transmission-prob
  report (nsmfmin + random-float diffmaxmin)
end
to-report a912mnsff                            ;transmission prob from male fawn age 9/12 mo to nonsib female fawns
  let a912m item 6 contactmatrix  ;5
  let nsffmin item 13 a912m
  let nsffmax item 14 a912m
  let diffmaxmin ((nsffmax - nsffmin) * transmission-prob)
  set nsffmin nsffmin * transmission-prob
  report (nsffmin + random-float diffmaxmin)
end
to-report a912mgrf                            ;transmission prob from male fawn age 9/12 mo to other females in group
  let a912m item 6 contactmatrix  ;5
  let grfmin item 15 a912m
  let grfmax item 16 a912m
  let diffmaxmin ((grfmax - grfmin) * transmission-prob)
  set grfmin grfmin * transmission-prob
  report (grfmin + random-float diffmaxmin)
end
to-report a912fmom                            ;transmission prob from female fawn age 9/12 mo to mom
  let a912f item 7 contactmatrix ;6
  let mommin item 1 a912f
  let mommax item 2 a912f
  let diffmaxmin ((mommax - mommin) * transmission-prob)
  set mommin mommin * transmission-prob
  report (mommin + random-float diffmaxmin)
end
to-report a912ffsm                            ;transmission prob from female fawn age 9/12 mo to male fullsibs
  let a912f item 7 contactmatrix ;6
  let fsmmin item 7 a912f
  let fsmmax item 8 a912f
  let diffmaxmin ((fsmmax - fsmmin) * transmission-prob)
  set fsmmin fsmmin * transmission-prob
  report (fsmmin + random-float diffmaxmin)
end
to-report a912ffsf                            ;transmission prob from female fawn age 9/12 mo to female fullsibs
  let a912f item 7 contactmatrix  ;6
  let fsfmin item 9 a912f
  let fsfmax item 10 a912f
  let diffmaxmin ((fsfmax - fsfmin) * transmission-prob)
  set fsfmin fsfmin * transmission-prob
  report (fsfmin + random-float diffmaxmin)
end
to-report a912fnsmf                            ;transmission prob from female fawn age 9/12 mo to nonsib male fawns
  let a912f item 7 contactmatrix ;6
  let nsmfmin item 11 a912f
  let nsmfmax item 12 a912f
  let diffmaxmin ((nsmfmax - nsmfmin) * transmission-prob)
  set nsmfmin nsmfmin * transmission-prob
  report (nsmfmin + random-float diffmaxmin)
end
to-report a912fnsff                            ;transmission prob from female fawn age 9/12 mo to nonsib female fawns
  let a912f item 7 contactmatrix ;6
  let nsffmin item 13 a912f
  let nsffmax item 14 a912f
  let diffmaxmin ((nsffmax - nsffmin) * transmission-prob)
  set nsffmin nsffmin * transmission-prob
  report (nsffmin + random-float diffmaxmin)
end
to-report a912fgrf                            ;transmission prob from female fawn age 9/12 mo to other females in group
  let a912f item 7 contactmatrix ;6
  let grfmin item 15 a912f
  let grfmax item 16 a912f
  let diffmaxmin ((grfmax - grfmin) * transmission-prob)
  set grfmin grfmin * transmission-prob
  report (grfmin + random-float diffmaxmin)
end
to-report a912fngrf                            ;transmission prob from female fawn age 9/12 mo to nongr females in nhood
  let a912f item 7 contactmatrix  ;6
  let ngrfmin item 17 a912f
  let ngrfmax item 18 a912f
  let diffmaxmin ((ngrfmax - ngrfmin) * transmission-prob)
  set ngrfmin ngrfmin * transmission-prob
  report (ngrfmin + random-float diffmaxmin)
end
to-report a1314mfsm                            ;transmission prob from male yearling age 13/14 mo to male fullsibs
  let a1314m item 8 contactmatrix  ;7
  let fsmmin item 7 a1314m
  let fsmmax item 8 a1314m
  let diffmaxmin ((fsmmax - fsmmin) * transmission-prob)
  set fsmmin fsmmin * transmission-prob
  report (fsmmin + random-float diffmaxmin)
end
to-report a1314mfsf                            ;transmission prob from male yearling age 13/14 mo to female fullsibs
  let a1314m item 8 contactmatrix   ;7
  let fsfmin item 9 a1314m
  let fsfmax item 10 a1314m
  let diffmaxmin ((fsfmax - fsfmin) * transmission-prob)
  set fsfmin fsfmin * transmission-prob
  report (fsfmin + random-float diffmaxmin)
end
to-report a1314mnsmy                            ;transmission prob from male yearling age 13/14 mo to nonsib male yearlings
  let a1314m item 8 contactmatrix  ;7
  let nsmfmin item 11 a1314m
  let nsmfmax item 12 a1314m
  let diffmaxmin ((nsmfmax - nsmfmin) * transmission-prob)
  set nsmfmin nsmfmin * transmission-prob
  report (nsmfmin + random-float diffmaxmin)
end
to-report a1314mnsfy                            ;transmission prob from male yearling age 13/14 mo to nonsib female yearlings
  let a1314m item 8 contactmatrix ;7
  let nsffmin item 13 a1314m
  let nsffmax item 14 a1314m
  let diffmaxmin ((nsffmax - nsffmin) * transmission-prob)
  set nsffmin nsffmin * transmission-prob
  report (nsffmin + random-float diffmaxmin)
end
to-report a1314mngrf                            ;transmission prob from male yearling age 13/14 mo to nongr females in nhood
  let a1314m item 8 contactmatrix ;7
  let ngrfmin item 17 a1314m
  let ngrfmax item 18 a1314m
  let diffmaxmin ((ngrfmax - ngrfmin) * transmission-prob)
  set ngrfmin ngrfmin * transmission-prob
  report (ngrfmin + random-float diffmaxmin)
end
to-report a1314ffsm                            ;transmission prob from female yearling age 13/14 mo to male fullsibs
  let a1314f item 9 contactmatrix ;8
  let fsmmin item 7 a1314f
  let fsmmax item 8 a1314f
  let diffmaxmin ((fsmmax - fsmmin) * transmission-prob)
  set fsmmin fsmmin * transmission-prob
  report (fsmmin + random-float diffmaxmin)
end
to-report a1314ffsf                            ;transmission prob from female yearling age 13/14 mo to female fullsibs
  let a1314f item 9 contactmatrix ;8
  let fsfmin item 9 a1314f
  let fsfmax item 10 a1314f
  let diffmaxmin ((fsfmax - fsfmin) * transmission-prob)
  set fsfmin fsfmin * transmission-prob
  report (fsfmin + random-float diffmaxmin)
end
to-report a1314fnsmy                            ;transmission prob from female yearling age 13/14 mo to nonsib male yearlings
  let a1314f item 9 contactmatrix ;8
  let nsmfmin item 11 a1314f
  let nsmfmax item 12 a1314f
  let diffmaxmin ((nsmfmax - nsmfmin) * transmission-prob)
  set nsmfmin nsmfmin * transmission-prob
  report (nsmfmin + random-float diffmaxmin)
end
to-report a1314fnsfy                            ;transmission prob from female yearling age 13/14 mo to nonsib female yearlings
  let a1314f item 9 contactmatrix ;8
  let nsffmin item 13 a1314f
  let nsffmax item 14 a1314f
  let diffmaxmin ((nsffmax - nsffmin) * transmission-prob)
  set nsffmin nsffmin * transmission-prob
  report (nsffmin + random-float diffmaxmin)
end
to-report a1314fgrf                            ;transmission prob from female yearling age 13/14 mo to young ad gr females
  let a1314f item 9 contactmatrix  ;8
  let grfmin item 15 a1314f
  let grfmax item 16 a1314f
  let diffmaxmin ((grfmax - grfmin) * transmission-prob)
  set grfmin grfmin * transmission-prob
  report (grfmin + random-float diffmaxmin)
end
to-report a1314fof                            ;transmission prob from female yearling age 13/14 mo to own fawns
  let a1314f item 10 contactmatrix ;9
  let ofmin item 3 a1314f
  let ofmax item 4 a1314f
  let diffmaxmin ((ofmax - ofmin) * transmission-prob)
  set ofmin ofmin * transmission-prob
  report (ofmin + random-float diffmaxmin)
end
to-report a15fmom                            ;transmission prob from female y age 15 mo to mom
  let a15f item 11 contactmatrix  ;10
  let mommin item 1 a15f
  let mommax item 2 a15f
  let diffmaxmin ((mommax - mommin) * transmission-prob)
  set mommin mommin * transmission-prob
  report (mommin + random-float diffmaxmin)
end
to-report a15fof                            ;transmission prob from female yearling age 15 mo to own fawns
  let a15f item 11 contactmatrix  ;10
  let ofmin item 3 a15f
  let ofmax item 4 a15f
  let diffmaxmin ((ofmax - ofmin) * transmission-prob)
  set ofmin ofmin * transmission-prob
  report (ofmin + random-float diffmaxmin)
end
to-report a15ffsf                            ;transmission prob from female yearling age 15 mo to female fullsibs
  let a15f item 11 contactmatrix  ;10
  let fsfmin item 9 a15f
  let fsfmax item 10 a15f
  let diffmaxmin ((fsfmax - fsfmin) * transmission-prob)
  set fsfmin fsfmin * transmission-prob
  report (fsfmin + random-float diffmaxmin)
end
to-report a15fgrf                            ;transmission prob from female yearling age 15 mo to young ad gr females
  let a15f item 11 contactmatrix ;10
  let grfmin item 15 a15f
  let grfmax item 16 a15f
  let diffmaxmin ((grfmax - grfmin) * transmission-prob)
  set grfmin grfmin * transmission-prob
  report (grfmin + random-float diffmaxmin)
end
to-report a15fngrf                            ;transmission prob from female yearling age 15 mo to non gr females
  let a15f item 11 contactmatrix ;10
  let grfmin item 17 a15f
  let grfmax item 18 a15f
  let diffmaxmin ((grfmax - grfmin) * transmission-prob)
  set grfmin grfmin * transmission-prob
  report (grfmin + random-float diffmaxmin)
end
to-report a1618fmom                            ;transmission prob from female y age 16-18 mo to mom
  let a1618f item 12 contactmatrix  ;11
  let mommin item 1 a1618f
  let mommax item 2 a1618f
  let diffmaxmin ((mommax - mommin) * transmission-prob)
  set mommin mommin * transmission-prob
  report (mommin + random-float diffmaxmin)
end
to-report a1618fofm                          ;transmission prob from female yearling age 16-18 mo to own male fawns
  let a1618f item 12 contactmatrix  ;11
  let ofmmin item 3 a1618f
  let ofmmax item 4 a1618f
  let diffmaxmin ((ofmmax - ofmmin) * transmission-prob)
  set ofmmin ofmmin * transmission-prob
  report (ofmmin + random-float diffmaxmin)
end
to-report a1618foff                          ;transmission prob from female yearling age 16-18 mo to own female fawns
  let a1618f item 12 contactmatrix ;11
  let offmin item 5 a1618f
  let offmax item 6 a1618f
  let diffmaxmin ((offmax - offmin) * transmission-prob)
  set offmin offmin * transmission-prob
  report (offmin + random-float diffmaxmin)
end
to-report a1618ffsf                            ;transmission prob from female yearling age 16-18 mo to female fullsibs
  let a1618f item 12 contactmatrix  ;11
  let fsfmin item 9 a1618f
  let fsfmax item 10 a1618f
  let diffmaxmin ((fsfmax - fsfmin) * transmission-prob)
  set fsfmin fsfmin * transmission-prob
  report (fsfmin + random-float diffmaxmin)
end
to-report a1618fgrf                            ;transmission prob from female yearling age 16-18 mo to young ad gr females
  let a1618f item 12 contactmatrix ;11
  let grfmin item 15 a1618f
  let grfmax item 16 a1618f
  let diffmaxmin ((grfmax - grfmin) * transmission-prob)
  set grfmin grfmin * transmission-prob
  report (grfmin + random-float diffmaxmin)
end
to-report a1618fngrf                            ;transmission prob from female yearling age 16-18 mo to non gr females
  let a1618f item 12 contactmatrix  ;11
  let grfmin item 17 a1618f
  let grfmax item 18 a1618f
  let diffmaxmin ((grfmax - grfmin) * transmission-prob)
  set grfmin grfmin * transmission-prob
  report (grfmin + random-float diffmaxmin)
end
to-report a1920fmom                            ;transmission prob from female y age 19/20 mo to mom
  let a1920f item 13 contactmatrix  ;12
  let mommin item 1 a1920f
  let mommax item 2 a1920f
  let diffmaxmin ((mommax - mommin) * transmission-prob)
  set mommin mommin * transmission-prob
  report (mommin + random-float diffmaxmin)
end
to-report a1920fofm                          ;transmission prob from female yearling age 19/20 mo to own male fawns
  let a1920f item 13 contactmatrix  ;12
  let ofmmin item 3 a1920f
  let ofmmax item 4 a1920f
  let diffmaxmin ((ofmmax - ofmmin) * transmission-prob)
  set ofmmin ofmmin * transmission-prob
  report (ofmmin + random-float diffmaxmin)
end
to-report a1920foff                          ;transmission prob from female yearling age 19/20 mo to own female fawns
  let a1920f item 13 contactmatrix ;12
  let offmin item 5 a1920f
  let offmax item 6 a1920f
  let diffmaxmin ((offmax - offmin) * transmission-prob)
  set offmin offmin * transmission-prob
  report (offmin + random-float diffmaxmin)
end
to-report a1920ffsf                            ;transmission prob from female yearling age 19/20 mo to female fullsibs
  let a1920f item 13 contactmatrix  ;12
  let fsfmin item 9 a1920f
  let fsfmax item 10 a1920f
  let diffmaxmin ((fsfmax - fsfmin) * transmission-prob)
  set fsfmin fsfmin * transmission-prob
  report (fsfmin + random-float diffmaxmin)
end
to-report a1920fgrf                            ;transmission prob from female yearling age 19/20 mo to young ad gr females
  let a1920f item 13 contactmatrix ;12
  let grfmin item 15 a1920f
  let grfmax item 16 a1920f
  let diffmaxmin ((grfmax - grfmin) * transmission-prob)
  set grfmin grfmin * transmission-prob
  report (grfmin + random-float diffmaxmin)
end
to-report a1920fngrf                            ;transmission prob from female yearling age 19/20 mo to non gr females
  let a1920f item 13 contactmatrix ;12
  let grfmin item 17 a1920f
  let grfmax item 18 a1920f
  let diffmaxmin ((grfmax - grfmin) * transmission-prob)
  set grfmin grfmin * transmission-prob
  report (grfmin + random-float diffmaxmin)
end
to-report a2124fmom                            ;transmission prob from female y age 21/24 mo to mom
  let a2124f item 14 contactmatrix ;13
  let mommin item 1 a2124f
  let mommax item 2 a2124f
  let diffmaxmin ((mommax - mommin) * transmission-prob)
  set mommin mommin * transmission-prob
  report (mommin + random-float diffmaxmin)
end
to-report a2124fofm                          ;transmission prob from female yearling age 21/24 mo to own male fawns
  let a2124f item 14 contactmatrix ;13
  let ofmmin item 3 a2124f
  let ofmmax item 4 a2124f
  let diffmaxmin ((ofmmax - ofmmin) * transmission-prob)
  set ofmmin ofmmin * transmission-prob
  report (ofmmin + random-float diffmaxmin)
end
to-report a2124foff                          ;transmission prob from female yearling age 21/24 mo to own female fawns
  let a2124f item 14 contactmatrix  ;13
  let offmin item 5 a2124f
  let offmax item 6 a2124f
  let diffmaxmin ((offmax - offmin) * transmission-prob)
  set offmin offmin * transmission-prob
  report (offmin + random-float diffmaxmin)
end
to-report a2124ffsf                            ;transmission prob from female yearling age 21/24 mo to female fullsibs
  let a2124f item 14 contactmatrix  ;13
  let fsfmin item 9 a2124f
  let fsfmax item 10 a2124f
  let diffmaxmin ((fsfmax - fsfmin) * transmission-prob)
  set fsfmin fsfmin * transmission-prob
  report (fsfmin + random-float diffmaxmin)
end
to-report a2124fgrf                            ;transmission prob from female yearling age 21/24 mo to young ad gr females
  let a2124f item 14 contactmatrix  ;13
  let grfmin item 15 a2124f
  let grfmax item 16 a2124f
  let diffmaxmin ((grfmax - grfmin) * transmission-prob)
  set grfmin grfmin * transmission-prob
  report (grfmin + random-float diffmaxmin)
end
to-report a2124fngrf                            ;transmission prob from female yearling age 21/24 mo to non gr females
  let a2124f item 14 contactmatrix  ;13
  let grfmin item 17 a2124f
  let grfmax item 18 a2124f
  let diffmaxmin ((grfmax - grfmin) * transmission-prob)
  set grfmin grfmin * transmission-prob
  report (grfmin + random-float diffmaxmin)
end
to-report a2526ffsf                            ;transmission prob from female  age 25/26 mo to female fullsibs
  let a2526f item 15 contactmatrix  ;14
  let fsfmin item 9 a2526f
  let fsfmax item 10 a2526f
  let diffmaxmin ((fsfmax - fsfmin) * transmission-prob)
  set fsfmin fsfmin * transmission-prob
  report (fsfmin + random-float diffmaxmin)
end
to-report a2526fgrf                            ;transmission prob from female age 25/26 mo to young ad gr females
  let a2526f item 15 contactmatrix  ;14
  let grfmin item 15 a2526f
  let grfmax item 16 a2526f
  let diffmaxmin ((grfmax - grfmin) * transmission-prob)
  set grfmin grfmin * transmission-prob
  report (grfmin + random-float diffmaxmin)
end
to-report a2526fmy                            ;transmission prob from female age 25/26 mo to male yearlings in mixed group
  let a2526f item 15 contactmatrix  ;14
  let mymin item 19 a2526f
  let mymax item 20 a2526f
  let diffmaxmin ((mymax - mymin) * transmission-prob)
  set mymin mymin * transmission-prob
  report (mymin + random-float diffmaxmin)
end
to-report a2526fof                            ;transmission prob from female age 25/26 mo to own fawns
  let a2526f item 16 contactmatrix ;15
  let ofmin item 3 a2526f
  let ofmax item 4 a2526f
  let diffmaxmin ((ofmax - ofmin) * transmission-prob)
  set ofmin ofmin * transmission-prob
  report (ofmin + random-float diffmaxmin)
end
to-report a27fd4mom                            ;transmission prob from ad female >26mo to mom Jan-Apr gest period
  let a27fd4 item 17 contactmatrix  ;16
  let mommin item 1 a27fd4
  let mommax item 2 a27fd4
  let diffmaxmin ((mommax - mommin) * transmission-prob)
  set mommin mommin * transmission-prob
  report (mommin + random-float diffmaxmin)
end
to-report a27fd10mom                            ;transmission prob from ad female >26mo to mom Jul-Oct pre-rut period
  let a27fd10 item 19 contactmatrix ;18
  let mommin item 1 a27fd10
  let mommax item 2 a27fd10
  let diffmaxmin ((mommax - mommin) * transmission-prob)
  set mommin mommin * transmission-prob
  report (mommin + random-float diffmaxmin)
end
to-report a27fd12mom                            ;transmission prob from ad female >26mo to mom Nov-Dec rut period
  let a27fd12 item 20 contactmatrix ;19
  let mommin item 1 a27fd12
  let mommax item 2 a27fd12
  let diffmaxmin ((mommax - mommin) * transmission-prob)
  set mommin mommin * transmission-prob
  report (mommin + random-float diffmaxmin)
end
to-report a27fd4ofm                            ;transmission prob from female age >26 mo to own m fawns; Jan-Apr gest period
  let a27fd4 item 17 contactmatrix ;16
  let ofmin item 3 a27fd4
  let ofmax item 4 a27fd4
  let diffmaxmin ((ofmax - ofmin) * transmission-prob)
  set ofmin ofmin * transmission-prob
  report (ofmin + random-float diffmaxmin)
end
to-report a27fd7ofm                            ;transmission prob from female age >26 mo to own m fawns; July weaning period
  let a27fd7 item 18 contactmatrix ;17
  let ofmin item 3 a27fd7
  let ofmax item 4 a27fd7
  let diffmaxmin ((ofmax - ofmin) * transmission-prob)
  set ofmin ofmin * transmission-prob
  report (ofmin + random-float diffmaxmin)
end
to-report a27fd10ofm                            ;transmission prob from female age >26 mo to own m fawns; aug-oct prerut period
  let a27fd10 item 19 contactmatrix ;18
  let ofmin item 3 a27fd10
  let ofmax item 4 a27fd10
  let diffmaxmin ((ofmax - ofmin) * transmission-prob)
  set ofmin ofmin * transmission-prob
  report (ofmin + random-float diffmaxmin)
end
to-report a27fd12ofm                            ;transmission prob from female age >26 mo to own m fawns; nov dec rut period
  let a27fd12 item 20 contactmatrix ;19
  let ofmin item 3 a27fd12
  let ofmax item 4 a27fd12
  let diffmaxmin ((ofmax - ofmin) * transmission-prob)
  set ofmin ofmin * transmission-prob
  report (ofmin + random-float diffmaxmin)
end
to-report a27fd4off                           ;transmission prob from female age >26 mo to own f fawns; Jan-Apr gest period
  let a27fd4 item 17 contactmatrix  ;16
  let ofmin item 5 a27fd4
  let ofmax item 6 a27fd4
  let diffmaxmin ((ofmax - ofmin) * transmission-prob)
  set ofmin ofmin * transmission-prob
  report (ofmin + random-float diffmaxmin)
end
to-report a27fd7off                            ;transmission prob from female age >26 mo to own f fawns; July weaning period
  let a27fd7 item 18 contactmatrix ;17
  let ofmin item 5 a27fd7
  let ofmax item 6 a27fd7
  let diffmaxmin ((ofmax - ofmin) * transmission-prob)
  set ofmin ofmin * transmission-prob
  report (ofmin + random-float diffmaxmin)
end
to-report a27fd10off                            ;transmission prob from female age >26 mo to own fawns; aug-oct prerut period
  let a27fd10 item 19 contactmatrix ;18
  let ofmin item 5 a27fd10
  let ofmax item 6 a27fd10
  let diffmaxmin ((ofmax - ofmin) * transmission-prob)
  set ofmin ofmin * transmission-prob
  report (ofmin + random-float diffmaxmin)
end
to-report a27fd12off                            ;transmission prob from female age >26 mo to own fawns; nov dec rut period
  let a27fd12 item 20 contactmatrix ;19
  let ofmin item 5 a27fd12
  let ofmax item 6 a27fd12
  let diffmaxmin ((ofmax - ofmin) * transmission-prob)
  set ofmin ofmin * transmission-prob
  report (ofmin + random-float diffmaxmin)
end
to-report a27fd10grf                            ;transmission prob from female age 25/26 mo to young ad gr females- prerut except nursing
  let a27fd10 item 17 contactmatrix  ;16
  let grfmin item 15 a27fd10
  let grfmax item 16 a27fd10
  let diffmaxmin ((grfmax - grfmin) * transmission-prob)
  set grfmin grfmin * transmission-prob
  report (grfmin + random-float diffmaxmin)
end
to-report a27fd12grf                            ;transmission prob from female age 25/26 mo to young ad gr females rut
  let a27fd12 item 20 contactmatrix  ;19
  let grfmin item 15 a27fd12
  let grfmax item 16 a27fd12
  let diffmaxmin ((grfmax - grfmin) * transmission-prob)
  set grfmin grfmin * transmission-prob
  report (grfmin + random-float diffmaxmin)
end
to-report a27fd10ngrf                            ;transmission prob from female age 25/26 mo to non gr females- prerut except nursing
  let a27fd10 item 17 contactmatrix ;16
  let grfmin item 17 a27fd10
  let grfmax item 18 a27fd10
  let diffmaxmin ((grfmax - grfmin) * transmission-prob)
  set grfmin grfmin * transmission-prob
  report (grfmin + random-float diffmaxmin)
end
to-report a27fd12ngrf                            ;transmission prob from female age 25/26 mo to ngr females rut
  let a27fd12 item 20 contactmatrix ;19
  let grfmin item 17 a27fd12
  let grfmax item 18 a27fd12
  let diffmaxmin ((grfmax - grfmin) * transmission-prob)
  set grfmin grfmin * transmission-prob
  report (grfmin + random-float diffmaxmin)
end
to-report a15mgrb                            ;transmission prob from ad male in b group to gr members
  let a15m item 21 contactmatrix  ;20
  let grbmin item 19 a15m
  let grbmax item 20 a15m
  let diffmaxmin ((grbmax - grbmin) * transmission-prob)
  set grbmin grbmin * transmission-prob
  report (grbmin + random-float diffmaxmin)
end
to-report a15msb                            ;transmission prob from solitary ad male to other male in vicinity
  let a15m item 22 contactmatrix ;21
  let sbmin item 19 a15m
  let sbmax item 20 a15m
  let diffmaxmin ((sbmax - sbmin) * transmission-prob)
  set sbmin sbmin * transmission-prob
  report (sbmin + random-float diffmaxmin)
end
@#$#@#$#@
GRAPHICS-WINDOW
788
103
1248
524
-1
-1
4.0
1
10
1
1
1
0
1
1
1
0
112
0
102
1
1
1
ticks
30.0

BUTTON
451
69
515
102
Setup
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
223
670
293
715
total deer
count deers
17
1
11

MONITOR
222
572
294
617
male deer
count deers with [sex = 1]
17
1
11

MONITOR
223
621
293
666
female deer
count deers with [sex = 2]
17
1
11

SLIDER
149
137
282
170
mf6nhm
mf6nhm
0
.1
0.055
.001
1
NIL
HORIZONTAL

SLIDER
149
174
282
207
ff6nhm
ff6nhm
0
.1
0.055
.001
1
NIL
HORIZONTAL

SLIDER
149
212
283
245
mf12nhm
mf12nhm
0
1
0.05
.001
1
NIL
HORIZONTAL

SLIDER
150
249
283
282
ff12nhm
ff12nhm
0
1
0.05
0.001
1
NIL
HORIZONTAL

SLIDER
290
290
427
323
myhm
myhm
0
1
0.25
.01
1
NIL
HORIZONTAL

SLIDER
291
328
427
361
fyhm
fyhm
0
1
0.15
.01
1
NIL
HORIZONTAL

SLIDER
150
288
283
321
mynhm
mynhm
0
1
0.01
.01
1
NIL
HORIZONTAL

SLIDER
151
326
282
359
fynhm
fynhm
0
1
0.0
0.01
1
NIL
HORIZONTAL

SLIDER
291
366
429
399
mahm
mahm
0
1
0.4
.01
1
NIL
HORIZONTAL

SLIDER
291
403
430
436
fahm
fahm
0
1
0.2
.001
1
NIL
HORIZONTAL

SLIDER
152
363
283
396
manhm
manhm
0
1
0.01
.01
1
NIL
HORIZONTAL

SLIDER
150
400
284
433
fanhm
fanhm
0
1
0.02
.01
1
NIL
HORIZONTAL

BUTTON
569
70
632
103
Go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
46
575
206
708
deer population
months
deer
0.0
1000.0
0.0
15000.0
true
false
"plotxy ticks count deers" ""
PENS
"pen-0" 1.0 0 -5298144 true "" ""

SLIDER
293
136
423
169
mf6hm
mf6hm
0
1
0.0
.01
1
NIL
HORIZONTAL

SLIDER
292
176
424
209
ff6hm
ff6hm
0
1
0.0
.01
1
NIL
HORIZONTAL

SLIDER
291
215
425
248
mf12hm
mf12hm
0
1
0.05
.01
1
NIL
HORIZONTAL

SLIDER
292
251
427
284
ff12hm
ff12hm
0
1
0.02
.01
1
NIL
HORIZONTAL

MONITOR
575
793
711
838
CWD detection probability
pdcwd
0
1
11

MONITOR
786
53
843
98
Year
year
17
1
11

MONITOR
853
53
910
98
Month
remainder ticks 12 + 1
17
1
11

PLOT
218
717
378
837
Doe group size
NIL
NIL
2.0
15.0
0.0
10.0
true
false
"set-plot-x-range 2 15" ""
PENS
"default" 1.0 1 -16777216 true "" "histogram [gr + 1] of deers with [gl = 1]"

PLOT
46
715
206
835
Bachelor group size
NIL
NIL
0.0
15.0
0.0
10.0
true
false
"set-plot-x-range 0 15" ""
PENS
"default" 1.0 1 -16777216 true "" "histogram [gr] of deers with [ml > 0]"

MONITOR
573
574
711
619
CWD infected deer
count deers with [cwd = 1]
17
1
11

SLIDER
436
212
656
245
%fawn-male-harvest-tested
%fawn-male-harvest-tested
0
1
0.05
0.01
1
NIL
HORIZONTAL

SLIDER
434
291
655
324
%yearling-male-harvest-tested
%yearling-male-harvest-tested
0
1
0.2
0.1
1
NIL
HORIZONTAL

SLIDER
436
364
657
397
%adult-male-harvest-tested
%adult-male-harvest-tested
0
1
0.1
0.1
1
NIL
HORIZONTAL

SLIDER
436
253
657
286
%fawn-female-harvest-tested
%fawn-female-harvest-tested
0
1
0.05
0.01
1
NIL
HORIZONTAL

SLIDER
435
328
658
361
%yearling-female-harvest-tested
%yearling-female-harvest-tested
0
1
0.2
0.1
1
NIL
HORIZONTAL

SLIDER
436
400
657
433
%adult-female-harvest-tested
%adult-female-harvest-tested
0
1
0.1
0.1
1
NIL
HORIZONTAL

CHOOSER
5
64
147
109
cwd_region
cwd_region
"Kankakee"
0

TEXTBOX
158
94
293
127
Non-hunting mortality\n(monthly rates)
13
15.0
1

TEXTBOX
307
93
406
125
Hunting mortality\n(annual rates)
13
15.0
1

TEXTBOX
20
149
135
174
Young male fawns
13
0.0
0

TEXTBOX
11
184
131
209
Young female fawns
13
0.0
1

TEXTBOX
24
221
133
239
Older male fawns
13
0.0
1

TEXTBOX
15
257
132
280
Older female fawns
13
0.0
1

TEXTBOX
44
289
130
321
Male yearlings
13
0.0
1

TEXTBOX
31
328
131
353
Female yearlings
13
0.0
1

TEXTBOX
63
365
138
383
Male adults
13
0.0
1

TEXTBOX
48
404
134
436
Female adults
13
0.0
1

MONITOR
574
625
711
670
CWD true prevalence
precision (count deers with [cwd = 1] / count deers) 3
17
1
11

MONITOR
574
682
710
727
CWD area (square miles)
cwd_area
17
1
11

PLOT
388
631
565
772
CWD prevalence
Months
True_prevalence
0.0
10.0
0.0
0.01
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot precision (count deers with [cwd = 1] / count deers) 3"

TEXTBOX
50
439
631
467
-------------------------------------------------------------------------------------
20
0.0
1

SLIDER
387
494
559
527
seed-infection
seed-infection
0
100
12.0
1
1
NIL
HORIZONTAL

TEXTBOX
46
535
626
553
------------------------------------------------------------------------------------
20
0.0
1

TEXTBOX
57
37
207
55
NIL
14
0.0
1

TEXTBOX
9
18
681
56
INdiana Odocoileus virginianus CWD dynamics (INOvCWD) version 2.2.0
20
0.0
1

MONITOR
574
739
710
784
NIL
max-dist-spark
17
1
11

CHOOSER
106
486
300
531
CWD_introduced_by
CWD_introduced_by
"adult-deer" "dispersing-male-yearling" "dispersing-female-yearling" "doe-groupmember" "doe-solitary" "buck-groupmember" "buck-solitary"
1

TEXTBOX
201
460
496
492
CWD introduction in the model deer population
13
15.0
1

PLOT
46
852
227
1006
doe matings
Number of matings
Number of does
1.0
5.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -16777216 true "" ""

SLIDER
458
149
662
182
targeted-culling-percentage
targeted-culling-percentage
0
1
0.0
0.05
1
NIL
HORIZONTAL

SWITCH
977
52
1113
85
export-rasters
export-rasters
0
1
-1000

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

deer
false
0
Polygon -7500403 true true 195 210 210 255 195 240 180 195 165 165 135 165 105 165 75 165 72 211 60 210 60 180 45 150 45 120 30 90 45 105 180 105 225 45 225 60 270 90 255 90 225 90 180 150
Polygon -7500403 true true 73 210 86 251 75 240 60 210
Polygon -7500403 true true 45 105 30 75 30 90 45 105 60 120 45 120
Line -7500403 true 210 60 165 15
Line -7500403 true 225 60 255 45
Line -7500403 true 195 45 210 15
Line -7500403 true 255 45 255 30
Line -7500403 true 255 45 270 30
Line -7500403 true 195 15 180 30

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="kank_bl_100" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="manhm">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mf6nhm">
      <value value="0.055"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ff12hm">
      <value value="0.13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adult-male-harvest-tested">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targeted-culling-percentage">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ff6nhm">
      <value value="0.055"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="CWD_introduced_by">
      <value value="&quot;dispersing-male-yearling&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ff12nhm">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mf12nhm">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%yearling-female-harvest-tested">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%adult-female-harvest-tested">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed-infection">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fahm">
      <value value="0.12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%fawn-male-harvest-tested">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fynhm">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mahm">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ff6hm">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="myhm">
      <value value="0.13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%fawn-female-harvest-tested">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fanhm">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mf12hm">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fyhm">
      <value value="0.21"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mynhm">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="%yearling-male-harvest-tested">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mf6hm">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cwd_region">
      <value value="&quot;Kankakee&quot;"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
