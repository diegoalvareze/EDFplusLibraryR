#   Created by Diego Alvarez-Estevez (http://dalvarezestevez.com)
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use these files except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
EDFreadHeader <- function (filename) {

fid <- file(filename, "rb") #0pen to read in binary mode

general_header_size <- 256#bytes
one_signal_header_size <- 256#bytes

# Reading general header
header <- list(version = as.integer(rawToChar(readBin(con = fid, what = raw(), n = 8, size = 1, signed = FALSE, endian = "little")))) #use raw to read bytes
header[["local_patient_identification"]] <- rawToChar(readBin(fid, raw(), 80, size = 1, signed = FALSE, endian = "little"))
header[["local_recording_identification"]] <- rawToChar(readBin(fid, raw(), 80, size = 1, signed = FALSE, endian = "little"))
header[["startdate_recording"]] <- rawToChar(readBin(fid, raw(), 8, size = 1, signed = FALSE, endian = "little"))
header[["starttime_recording"]] <- rawToChar(readBin(fid, raw(), 8, size = 1, signed = FALSE, endian = "little"))
header[["num_bytes_header"]] <- as.integer(rawToChar(readBin(fid, raw(), 8, size = 1, signed = FALSE, endian = "little")))
header[["reserved"]] <- rawToChar(readBin(fid, raw(), 44, size = 1, signed = FALSE, endian = "little"))
header[["num_data_records"]] <- as.integer(rawToChar(readBin(fid, raw(), 8, size = 1, signed = FALSE, endian = "little")))
header[["duration_data_record"]] <- as.integer(rawToChar(readBin(fid, raw(), 8, size = 1, signed = FALSE, endian = "little")))
header[["num_signals"]] <- as.integer(rawToChar(readBin(fid, raw(), 4, size = 1, signed = FALSE, endian = "little")))

# Reading signal dependent header
header[["signals_label"]] <- character(header$num_signals)
for (i in 1:header$num_signals) header$signals_label[i] <- rawToChar(readBin(fid, raw(), 16, size = 1, signed = FALSE, endian = "little"))
header[["signals_transducer_type"]] <- character(header$num_signals)
for (i in 1:header$num_signals) header$signals_transducer_type[i] <- rawToChar(readBin(fid, raw(), 80, size = 1, signed = FALSE, endian = "little"))
header[["signals_physical_dimension"]] <- character(header$num_signals)
for (i in 1:header$num_signals) header$signals_physical_dimension[i] <- rawToChar(readBin(fid, raw(), 8, size = 1, signed = FALSE, endian = "little"))
header[["signals_physical_min"]] <- character(header$num_signals)
for (i in 1:header$num_signals) header$signals_physical_min[i] <- as.integer(rawToChar(readBin(fid, raw(), 8, size = 1, signed = FALSE, endian = "little")))
header[["signals_physical_max"]] <- character(header$num_signals)
for (i in 1:header$num_signals) header$signals_physical_max[i] <- as.integer(rawToChar(readBin(fid, raw(), 8, size = 1, signed = FALSE, endian = "little")))
header[["signals_digital_min"]] <- character(header$num_signals)
for (i in 1:header$num_signals) header$signals_digital_min[i] <- as.integer(rawToChar(readBin(fid, raw(), 8, size = 1, signed = FALSE, endian = "little")))
header[["signals_digital_max"]] <- character(header$num_signals)
for (i in 1:header$num_signals) header$signals_digital_max[i] <- as.integer(rawToChar(readBin(fid, raw(), 8, size = 1, signed = FALSE, endian = "little")))
header[["signals_prefiltering"]] <- character(header$num_signals)
for (i in 1:header$num_signals) header$signals_prefiltering[i] <- rawToChar(readBin(fid, raw(), 80, size = 1, signed = FALSE, endian = "little"))
header[["signals_num_samples_datarecord"]] <- integer(header$num_signals) #values default to zero
header[["signals_signalOffset"]] <- integer(header$num_signals)
header[["signals_sample_rate"]] <- integer(header$num_signals)
for (i in 1:header$num_signals) 
{
  header$signals_num_samples_datarecord[i] <- as.integer(rawToChar(readBin(fid, raw(), 8, size = 1, signed = FALSE, endian = "little")))
  # NOTE: The two following are not specific EDF header fields, but are practical for EDF handling
  if (header$duration_data_record > 0)
    header$signals_sample_rate[i] <- header$signals_num_samples_datarecord[i] / header$duration_data_record
  else
    header$signals_sample_rate[i] <- -1
  
  if (i > 1)
    header$signals_signalOffset[i] <- header$signals_signalOffset[i - 1] + 2 * header$signals_num_samples_datarecord[i - 1]
}
header[["signals_reserved"]] <- character(header$num_signals)
for (i in 1:header$num_signals) header$signals_reserved[i] <- rawToChar(readBin(fid, raw(), 32, size = 1, signed = FALSE, endian = "little"))

close(fid)

return(header)
}
