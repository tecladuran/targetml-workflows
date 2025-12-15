
samples_directory <- "/storage/projects/TargetML/new_samples_22/day_1"
FileName <- sort(list.files(samples_directory, full.names = FALSE))
SampleID <- FileName

day_1 <- data.frame(
  FileName = FileName,
  SampleID = SampleID,
  class = c(
    "M1", "blank",
    "control", "blank",
    "patient", "blank",
    "patient", "blank",
    "patient", "blank",
    "patient", "blank",
    "patient", "blank",
    "control", "blank",
    "patient", "blank",
    "patient", "blank",
    "patient", "blank",
    "patient", "blank",
    "control", "blank",
    "control", "blank"
  ),
  matrix = c(
    "M1", NA,
    "SU", NA,
    NA, NA,      # 391
    NA, NA,      # 440
    NA, NA,      # 398
    NA, NA,      # 485
    NA, NA,      # 369
    "SU", NA,
    NA, NA,      # 450
    NA, NA,      # 430
    NA, NA,      # 476
    NA, NA,      # 312
    "SU", NA,
    "SU", NA
  ),
  patient_id = c(
    NA, NA,
    NA, NA,
    "391", NA,
    "440", NA,
    "398", NA,
    "485", NA,
    "369", NA,
    NA, NA,
    "450", NA,
    "430", NA,
    "476", NA,
    "312", NA,
    NA, NA,
    NA, NA
  ),
  patient_condition = c(
    NA, NA,
    NA, NA,
    "CRC", NA,
    "CTRL", NA,
    "Adenoma", NA,
    "CRC", NA,
    "CTRL", NA,
    NA, NA,
    "Adenoma", NA,
    "CTRL", NA,
    "CTRL", NA,
    "CRC", NA,
    NA, NA,
    NA, NA
  ),
  control_level = c(
    NA, NA,
    "MQC", NA,
    NA, NA,
    NA, NA,
    NA, NA,
    NA, NA,
    NA, NA,
    "LQC", NA,
    NA, NA,
    NA, NA,
    NA, NA,
    NA, NA,
    "HQC", NA,
    "LLQC", NA
  ),
  day = rep(1, 28)
)



samples_directory <- "/storage/projects/TargetML/new_samples_22/day_2"
FileName <- sort(list.files(samples_directory, full.names = FALSE))
SampleID <- FileName

samples_directory <- "/storage/projects/TargetML/new_samples_22/day_2"
FileName <- sort(list.files(samples_directory, full.names = FALSE))
SampleID <- FileName

day_2 <- data.frame(
  FileName = FileName,
  SampleID = SampleID,
  class = c(
    "M1", "blank",
    "control", "blank",
    "patient", "blank",
    "patient", "blank",
    "patient", "blank",
    "patient", "blank",
    "patient", "blank",
    "control", "blank",
    "patient", "blank",
    "patient", "blank",
    "patient", "blank",
    "patient", "blank",
    "control", "blank",
    "control", "blank"
  ),
  matrix = c(
    "M1", NA,
    "SU", NA,
    NA, NA,       # 242
    NA, NA,       # 295
    NA, NA,       # 409
    NA, NA,       # 288
    NA, NA,       # 231
    "SU", NA,
    NA, NA,       # 333
    NA, NA,       # 339
    NA, NA,       # 507
    NA, NA,       # 291
    "SU", NA,
    "SU", NA
  ),
  patient_id = c(
    NA, NA,
    NA, NA,
    "242", NA,
    "295", NA,
    "409", NA,
    "288", NA,
    "231", NA,
    NA, NA,
    "333", NA,
    "339", NA,
    "507", NA,
    "291", NA,
    NA, NA,
    NA, NA
  ),
  patient_condition = c(
    NA, NA,
    NA, NA,
    "Adenoma", NA,
    "Adenoma", NA,
    "CTRL", NA,
    "Adenoma", NA,
    "CRC", NA,
    NA, NA,
    "Adenoma", NA,
    "CTRL", NA,
    "CRC", NA,
    "CRC", NA,
    NA, NA,
    NA, NA
  ),
  control_level = c(
    NA, NA,
    "MQC", NA,
    NA, NA,
    NA, NA,
    NA, NA,
    NA, NA,
    NA, NA,
    "LQC", NA,
    NA, NA,
    NA, NA,
    NA, NA,
    NA, NA,
    "HQC", NA,
    "LLQC", NA
  ),
  day = rep(2, 28)
)

df <- rbind(day_1, day_2)
df$matrix[df$class == "patient"] <- "urine"
write.csv(df, "annotations_patient_measurements_22.csv", row.names = FALSE)

