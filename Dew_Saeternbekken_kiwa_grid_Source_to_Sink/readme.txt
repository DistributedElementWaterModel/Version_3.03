
#Generate watercourse hierarchy
riverFlow control_flow.txt

#Define model domain for distributed element water balance model
stationMask control_mask.txt

#Preprocessing for distributed element water balance model
predew control_pre.txt

#Run distributed element water balance model
dew control_dew.txt

#Source-to-sink routing model results in file sts_00080005.var

