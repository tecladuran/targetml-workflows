# CAL 1
samples_directory <- "/storage/projects/TargetML/tgn/Tarragona Pool_Calibration/Calib1"
FileName <- sort(list.files(samples_directory, full.names = FALSE))
SampleID <- FileName
cal1 <- data.frame(
  FileName = FileName,
  SampleID = SampleID,
  cal = rep(1, 22),
  class = c(
    "M1", "blank", 
    "ppb0.1", "blank", 
    "ppb10", "blank", 
    "ppb30", "blank", 
    "ppb1", "blank", 
    "ppb0.5", "blank", 
    "ppb7.5", "blank", 
    "ppb0", "blank", 
    "ppb20", "blank", 
    "ppb5", "blank", 
    "ppb15", "blank"
  )
)

# CAL 2
samples_directory <- "/storage/projects/TargetML/tgn/Tarragona Pool_Calibration/Calib2"
FileName <- sort(list.files(samples_directory, full.names = FALSE))
SampleID <- FileName
cal2 <- data.frame(
  FileName = FileName,
  SampleID = SampleID,
  cal = rep(2, 22),
  class = c(
    "M1", "blank", 
    "ppb7.5", "blank", 
    "ppb5", "blank", 
    "ppb0.1", "blank", 
    "ppb20", "blank", 
    "ppb15", "blank", 
    "ppb30", "blank", 
    "ppb0.5", "blank", 
    "ppb1", "blank", 
    "ppb10", "blank", 
    "ppb0", "blank"
  )
)

# CAL 3
samples_directory <- "/storage/projects/TargetML/tgn/Tarragona Pool_Calibration/Calib3"
FileName <- sort(list.files(samples_directory, full.names = FALSE))
SampleID <- FileName
cal3 <- data.frame(
  FileName = FileName,
  SampleID = SampleID,
  cal = rep(3, 22),
  class = c(
    "M1", "blank", 
    "ppb30", "blank", 
    "ppb0", "blank", 
    "ppb15", "blank", 
    "ppb0.5", "blank", 
    "ppb20", "blank", 
    "ppb1", "blank", 
    "ppb10", "blank", 
    "ppb7.5", "blank", 
    "ppb0.1", "blank", 
    "ppb5", "blank"
  )
)

# Merge and write CSV
df <- rbind(cal1, cal2, cal3)
write.csv(df, "annotations_pool_calibration_tgn.csv", row.names = FALSE)


