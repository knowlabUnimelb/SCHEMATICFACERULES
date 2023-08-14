# SCHEMATICFACERULES
This project contains data for the experiments published in:  Cheng and Little (2023)

Archive created 14-Aug-23

========================================================================
Folders:
========================================================================
Data
- Contains Raw Data from Cheng2023
------------------------------------------------------------------------

SchematicSFTdata
- \rawdata\ contains raw data from the SFT Logical Rules Experiment
-- columns are:
cols = {'sub', 'con', 'rot', 'ses', 'tri', 'itm', 'eyes', 'mouth', 'rsp', 'cat', 'acc', 'rt'};
-- rot = rotation (not used)
-- items are numbered as follows: 1 = HH, 2 = HL, 3 = LH, 4 = LL, 5 = Ex, 6 = Ix, 7 = Ey, 8 = Iy, 9 = R
-- eyes is value of top half
-- mouth is value of bottom half

- \Rdata\ is output from analysis for input to Houpt's SFT [R] package
-- cols are {'sub', 'con', 'rt', 'acc', 'itm'}

- \modeldata\ is output from analysis for input to DEMCMC modeling
-- files are .mat files
-- cols are {'sub', 'con', 'rot', 'ses', 'blk', 'tri', 'itm', 'eyes', 'mouth', 'rsp', 'cat', 'acc', 'rt'};

- Subject numbers correspond to the paper as follows:
Orientation Alignment   Code    Number
Upright     Aligned     UA1     503
Upright     Aligned     UA2     504
Upright     Aligned     UA3     505
Upright     Aligned     UA4     506
Upright     Misaligned  UM1     601
Upright     Misaligned  UM2     602
Upright     Misaligned  UM3     604
Upright     Misaligned  UM4     606
Inverted    Aligned     IA1     702
Inverted    Aligned     IA2     703
Inverted    Aligned     IA3     704
Inverted    Aligned     IA4     705
Inverted    Misaligned  IM1     802
Inverted    Misaligned  IM2     804
Inverted    Misaligned  IM3     806
Inverted    Misaligned  IM4     807


SchematicCompositeFaceData
- \rawdata\ contains raw data files from complete design composite face task
- Columns in data file are:
cols = {'Subject', 'Block', 'Number', 'Set', 'Cued', 'Resp', 'Congruent', 'Direction', 'Alignment', 'Study', 'Test', 'Correct', 'Response', 'RT'};

-- Number is Nominal Block (not used)
-- Set is Face Set (v1 or v2)
-- Cued is instructed attended face half (1 = Top, 0 = Bottom)
-- Resp is the correct reponse (1 = Same, 0 = Different)
-- Congruent is 1 if congruent trial, 0 if incongruent
-- Direction is 1 if upright, 0 if inverted
-- Alignment is 1 if aligned, 0 if split
-- Study is the study item
-- Test is the test item
-- Correct = 1 if correct, 0 if wrong
-- Response is the actual response = 1 if Same, 0 if DIff
-- RT is response time

- \csv_files\ contain the rawdata in csvfile format

FaceCompMDS [Turk]
- Main data used to derive the coordinates used in the modeling in the paper
