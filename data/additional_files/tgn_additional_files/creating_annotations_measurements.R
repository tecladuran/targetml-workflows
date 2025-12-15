# DAY 4
samples_directory <- "/storage/projects/TargetML/tgn/measurements/Day_4"
FileName <- sort(list.files(samples_directory, full.names = FALSE))
SampleID <- FileName
day_4<- data.frame(
  FileName = FileName,
  SampleID = SampleID,
  class = c(
    "M1", "blank",
    "control", "blank",
    "control", "blank",
    "patient", "blank",
    "control", "blank",
    "patient", "blank",
    "patient", "blank",
    "patient", "blank",
    "control", "blank",
    "control", "blank",
    "patient", "blank",
    "patient", "blank",
    "patient", "blank",
    "patient", "blank",
    "control", "blank",
    "control", "blank"
  ),
  matrix = c(
    "M1", NA,               # 01: Fresh M1, 02: blank
    "SU", NA,               # 03: MQC in SU, 04: blank
    "pool", NA,   # 05: MQC in pool, 06: blank
    NA, NA,                 # 07: Patient 247, 08: blank
    "pool", NA,   # 09: LLQC in pool, 10: blank
    NA, NA,                 # 11: Patient 325, 12: blank
    NA, NA,                 # 13: Patient 327, 14: blank
    NA, NA,                 # 15: Patient 361, 16: blank
    "SU", NA,               # 17: LQC in SU, 18: blank
    "pool", NA,   # 19: LQC in pool, 20: blank
    NA, NA,                 # 21: Patient 370, 22: blank
    NA, NA,                 # 23: Patient 380, 24: blank
    NA, NA,                 # 25: Patient 401, 26: blank
    NA, NA,                 # 27: Patient 402, 28: blank
    "SU", NA,               # 29: HQC in SU, 30: blank
    "pool", NA    # 31: HQC in pool, 32: blank
  ),
  patient_id = c(
    NA, NA,
    NA, NA,
    NA, NA,
    "247", NA,
    NA, NA,
    "325", NA,
    "327", NA,
    "361", NA,
    NA, NA,
    NA, NA,
    "370", NA,
    "380", NA,
    "401", NA,
    "402", NA,
    NA, NA,
    NA, NA
  ),
  patient_condition = c(
    NA, NA,
    NA, NA,
    NA, NA,
    "CRC", NA,
    NA, NA,
    "CRC", NA,
    "CRC", NA,
    "CTRL", NA,
    NA, NA,
    NA, NA,
    "CTRL", NA,
    "CTRL", NA,
    "CRC", NA,
    "CRC", NA,
    NA, NA,
    NA, NA
  ),
  control_level = c(
    NA, NA,
    "MQC", NA,
    "MQC", NA,
    NA, NA,
    "LLQC", NA,
    NA, NA,
    NA, NA,
    NA, NA,
    "LQC", NA,
    "LQC", NA,
    NA, NA,
    NA, NA,
    NA, NA,
    NA, NA,
    "HQC", NA,
    "HQC", NA
  ),
  day = rep(4, 32)
)


samples_directory <- "/storage/projects/TargetML/tgn/measurements/Day_5"
FileName <- sort(list.files(samples_directory, full.names = FALSE))
SampleID <- FileName

day_5 <- data.frame(
  FileName = FileName,
  SampleID = SampleID,
  class = c(
    "M1", "blank",
    "control", "blank",
    "control", "blank",
    "patient", "blank",
    "control", "blank",
    "patient", "blank",
    "patient", "blank",
    "patient", "blank",
    "control", "blank",
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
    "pool", NA,
    NA, NA,
    "pool", NA,
    NA, NA,
    NA, NA,
    NA, NA,
    "SU", NA,
    "pool", NA,
    NA, NA,
    NA, NA,
    NA, NA,
    NA, NA,
    "SU", NA,
    "pool", NA
  ),
  patient_id = c(
    NA, NA,
    NA, NA,
    NA, NA,
    "407", NA,
    NA, NA,
    "411", NA,
    "426", NA,
    "432", NA,
    NA, NA,
    NA, NA,
    "453", NA,
    "492", NA,
    "493", NA,
    "494", NA,
    NA, NA,
    NA, NA
  ),
  patient_condition = c(
    NA, NA,
    NA, NA,
    NA, NA,
    "CTRL", NA,
    NA, NA,
    "CRC", NA,
    "CTRL", NA,
    "CTRL", NA,
    NA, NA,
    NA, NA,
    "CRC", NA,
    "CTRL", NA,
    "CTRL", NA,
    "CRC", NA,
    NA, NA,
    NA, NA
  ),
  control_level = c(
    NA, NA,
    "MQC", NA,
    "MQC", NA,
    NA, NA,
    "LLQC", NA,
    NA, NA,
    NA, NA,
    NA, NA,
    "LQC", NA,
    "LQC", NA,
    NA, NA,
    NA, NA,
    NA, NA,
    NA, NA,
    "HQC", NA,
    "HQC", NA
  ),
  day = rep(5, 32)
)


samples_directory <- "/storage/projects/TargetML/tgn/measurements/Day_6"
FileName <- sort(list.files(samples_directory, full.names = FALSE))
SampleID <- FileName

day_6 <- data.frame(
  FileName = FileName,
  SampleID = SampleID,
  class = c(
    "M1", "blank",
    "control", "blank",
    "control", "blank",
    "patient", "blank",
    "patient", "blank",
    "patient", "blank",
    "patient", "blank",
    "control", "blank",
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
    "pool", NA,
    NA, NA,
    NA, NA,
    NA, NA,
    NA, NA,
    "SU", NA,
    "pool", NA,
    NA, NA,
    NA, NA,
    NA, NA,
    NA, NA,
    "SU", NA,
    "pool", NA
  ),
  patient_id = c(
    NA, NA,
    NA, NA,
    NA, NA,
    "514", NA,
    "527", NA,
    "531", NA,
    "568", NA,
    NA, NA,
    NA, NA,
    "590", NA,
    "604", NA,
    "618", NA,
    "628", NA,
    NA, NA,
    NA, NA
  ),
  patient_condition = c(
    NA, NA,
    NA, NA,
    NA, NA,
    "CTRL", NA,
    "CRC", NA,
    "CTRL", NA,
    "CTRL", NA,
    NA, NA,
    NA, NA,
    "CRC", NA,
    "CTRL", NA,
    "CRC", NA,
    "CTRL", NA,
    NA, NA,
    NA, NA
  ),
  control_level = c(
    NA, NA,
    "MQC", NA,
    "MQC", NA,
    NA, NA,
    NA, NA,
    NA, NA,
    NA, NA,
    "LQC", NA,
    "LQC", NA,
    NA, NA,
    NA, NA,
    NA, NA,
    NA, NA,
    "HQC", NA,
    "HQC", NA
  ),
  day = rep(6, 30)
)

samples_directory <- "/storage/projects/TargetML/tgn/measurements/Day_7"
FileName <- sort(list.files(samples_directory, full.names = FALSE))
SampleID <- FileName

day_7 <- data.frame(
  FileName = FileName,
  SampleID = SampleID,
  class = c(
    "M1", "blank",
    "control", "blank",
    "control", "blank",
    "patient", "blank",
    "patient", "blank",
    "patient", "blank",
    "patient", "blank",
    "control", "blank",
    "control", "blank",
    "patient", "blank",
    "patient", "blank",
    "control", "blank",
    "control", "blank"
  ),
  matrix = c(
    "M1", NA,
    "SU", NA,
    "pool", NA,
    NA, NA,
    NA, NA,
    NA, NA,
    NA, NA,
    "SU", NA,
    "pool", NA,
    NA, NA,
    NA, NA,
    "SU", NA,
    "pool", NA
  ),
  patient_id = c(
    NA, NA,
    NA, NA,
    NA, NA,
    "643", NA,
    "690", NA,
    "800", NA,
    "932", NA,
    NA, NA,
    NA, NA,
    "1044", NA,
    "1062", NA,
    NA, NA,
    NA, NA
  ),
  patient_condition = c(
    NA, NA,
    NA, NA,
    NA, NA,
    "CRC", NA,
    "CTRL", NA,
    "CRC", NA,
    "CTRL", NA,
    NA, NA,
    NA, NA,
    "CRC", NA,
    "CRC", NA,
    NA, NA,
    NA, NA
  ),
  control_level = c(
    NA, NA,
    "MQC", NA,
    "MQC", NA,
    NA, NA,
    NA, NA,
    NA, NA,
    NA, NA,
    "LQC", NA,
    "LQC", NA,
    NA, NA,
    NA, NA,
    "HQC", NA,
    "HQC", NA
  ),
  day = rep(7, 26)
)

df <- rbind(day_4, day_5, day_6, day_7)
df$matrix[df$class == "patient"] <- "urine"
write.csv(df, "annotations_patient_measurements_tgn.csv", row.names = FALSE)
