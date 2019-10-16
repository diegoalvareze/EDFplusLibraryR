# Saves a list of annoations (as returned by EDFreadAnnotations function) into an EDF+ file.
# Additional data needs to be provided via the parameters, such as the
# patient identification (patId), the recording identification (recId),
# start date of the recording (startdate), start time of the recording
# (starttime) and the output file name (filename)
#
#   annotations2EDFplus(annotations, patId, recId, startdate, starttime, filename)
#
# More info, see EDF+ specification at http://www.edfplus.info/
#
# See also: EDFreadAnnotations
#
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
annotations2EDFplus = function (annotations, patId, recId, startdate, starttime, filename)
{

edf_block_size_limit <- 61440# In bytes

# Compute how many characters we will need to store
writenChars <- 0
num_necessary_blocks <- 1 # At least one block is necessary

# Sort annotations before writing to file
sortIndx <- order(annotations$offset)
annotations$offset <- annotations$offset[sortIndx]
annotations$duration <- annotations$duration[sortIndx]
annotations$label <- annotations$label[sortIndx]

# Note: The time-keeping TAL for the first datablock it always assumes a zero onset,
#       thus the corresponding TAL is: +0'20''20''0'
data <- c(charToRaw('+0'), as.raw(c(20,20,0)));

writenChars <- length(data)

# Each annotation will be included in a separate TAL
for (k in 1:length(annotations$offset)){

    dataForThisAnnotation <- NULL

    # Note: Each offset starts with '+' and finishes with unprintable ASCII '21'
    if (annotations$offset[k] >= 0){
      dataForThisAnnotation <- c(dataForThisAnnotation, c(charToRaw('+'), charToRaw(as.character(annotations$offset[k])), as.raw(21)));
    } else {
      dataForThisAnnotation <- c(dataForThisAnnotation, c(charToRaw(as.character(annotations$offset[k])), as.raw(21)));
    }
    
    # Note: Each duration finishes with unprintable ASCII '20'
    dataForThisAnnotation <- c(dataForThisAnnotation, charToRaw(as.character(annotations$duration[k])), as.raw(20));
    # Note: Each annotation finishes with unprintable ASCII '20' and must not contain
    #       any '20' within. Also because we assume one TAL per annotation
    #       then the unprintable ASCII '0' is added to close the TAL
    dataForThisAnnotation <- c(dataForThisAnnotation, charToRaw(annotations$label[k]), as.raw(c(20,0)));

    # Notice each char "weights" 1-byte
    if (writenChars + length(dataForThisAnnotation) >= edf_block_size_limit){
        # This annotation needs to be set on a new block
        num_necessary_blocks <- num_necessary_blocks + 1;

        # Current data is filled with zeros to the end of the datablock
        data <- c(data, as.raw(integer(edf_block_size_limit - writenChars)));

        # Reset the "writtenChars" counter
        writenChars <- 0;

        # A a new time-keeping TAL is built, taking as reference the offset time of the curent annotation
        if (annotations$offset[k] >= 0){
          time_keep_tal <- c(charToRaw('+'), charToRaw(as.character(annotations$offset[k])), as.raw(c(20,20,0)));
        } else {
          time_keep_tal <- c(charToRaw(as.character(annotations$offset[k])), as.raw(c(20,20,0)));
        }
        # New time-keeping TAL is inserted preceding the annotation (first annotation of the new block)
        dataForThisAnnotation <- c(time_keep_tal, dataForThisAnnotation);
    }

    data <- c(data, dataForThisAnnotation);
    writenChars <- writenChars + length(dataForThisAnnotation)
}

#Note: Unused bytes of the 'EDF Annotations' signal in the remainder of the data record are also filled with 0-bytes
data <- c(data, as.raw(integer(edf_block_size_limit - writenChars)));

fid <- file(filename, "wb") #0pen to write binary mode

if (fid == -1){
  stop('Error creating output EDF+ file')
  return(-1)
}

general_header_size <- 256#bytes
one_signal_header_size <- 256#bytes

# Write edf

# FIXED HEADER
header$version <- 0
header$local_patient_identification <- patId
header$local_recording_identification <- recId
header$startdate_recording <- startdate
header$starttime_recording <- starttime
header$num_signals <- 1
header$num_bytes_header <- general_header_size + one_signal_header_size
header$reserved <- 'EDF+C'
# Note we assume an 'Annotations only' EDF+ file
header$duration_data_record <- 0
# Each character "weights" 1-byte, thus the total number of necessary bytes
# equal the total number of characters to write
header$num_data_records <- num_necessary_blocks

modeUseBytes <- TRUE;

writeChar(sprintf('%-8s', as.character(header$version)), fid, nchars = 8, eos = NULL, useBytes = modeUseBytes);
writeChar(sprintf('%-80s', header$local_patient_identification), fid, nchars = 80, eos = NULL, useBytes = modeUseBytes);
writeChar(sprintf('%-80s', header$local_recording_identification), fid, nchars = 80, eos = NULL, useBytes = modeUseBytes);
writeChar(sprintf('%-8s', header$startdate_recording), fid, nchars = 8, eos = NULL, useBytes = modeUseBytes);
writeChar(sprintf('%-8s', header$starttime_recording), fid, nchars = 8, eos = NULL, useBytes = modeUseBytes);
actualHeaderBytes <- general_header_size + one_signal_header_size*header$num_signals
if (actualHeaderBytes != header$num_bytes_header){
    warning('EDFwriteSignal: Warning, num_bytes_header does not match the actual number of header bytes')
}
writeChar(sprintf('%-8s', as.character(actualHeaderBytes)), fid, nchars = 8, eos = NULL, useBytes = modeUseBytes);
writeChar(sprintf('%-44s', header$reserved), fid, nchars = 44, eos = NULL, useBytes = modeUseBytes);
writeChar(sprintf('%-8s', as.character(header$num_data_records)), fid, nchars = 8, eos = NULL, useBytes = modeUseBytes);
writeChar(sprintf('%-8s', as.character(header$duration_data_record)), fid, nchars = 8, eos = NULL, useBytes = modeUseBytes);
writeChar(sprintf('%-4s', as.character(header$num_signals)), fid, nchars = 4, eos = NULL, useBytes = modeUseBytes);

# SIGNAL DEPENDENT HEADER
header$signals_label <- 'EDF Annotations'
header$signals_transducer_type <- ''
header$signals_physical_dimension <- ''
header$signals_physical_min <- -32768
header$signals_physical_max <- 32767
header$signals_digital_min <- -32768
header$signals_digital_max <- 32767
header$signals_prefiltering <- ''
# One sample is a 2-byte (samples are encoded in EDF(+) using 16 bit pieces), and each character is using 1-byte, thus each character "weights" half a sample
header$signals_num_samples_datarecord <- edf_block_size_limit/2
header$signals_reserved <- ''

writeChar(sprintf('%-16s', header$signals_label), fid, nchars = 16, eos = NULL, useBytes = modeUseBytes);
writeChar(sprintf('%-80s', header$signals_transducer_type), fid, nchars = 80, eos = NULL, useBytes = modeUseBytes);
writeChar(sprintf('%-8s', header$signals_physical_dimension), fid, nchars = 8, eos = NULL, useBytes = modeUseBytes);
writeChar(sprintf('%-8s', as.character(header$signals_physical_min)), fid, nchars = 8, eos = NULL, useBytes = modeUseBytes);
writeChar(sprintf('%-8s', as.character(header$signals_physical_max)), fid, nchars = 8, eos = NULL, useBytes = modeUseBytes);
writeChar(sprintf('%-8s', as.character(header$signals_digital_min)), fid, nchars = 8, eos = NULL, useBytes = modeUseBytes);
writeChar(sprintf('%-8s', as.character(header$signals_digital_max)), fid, nchars = 8, eos = NULL, useBytes = modeUseBytes);
writeChar(sprintf('%-80s', header$signals_prefiltering), fid, nchars = 80, eos = NULL, useBytes = modeUseBytes);
writeChar(sprintf('%-8s', as.character(header$signals_num_samples_datarecord)), fid, nchars = 8, eos = NULL, useBytes = modeUseBytes);
writeChar(sprintf('%-32s', header$signals_reserved), fid, nchars = 32, eos = NULL, useBytes = modeUseBytes);


# DATA WRITING
header_length <- general_header_size + header$num_signals * one_signal_header_size
current_position <- seek(fid, where = NA, rw = "write"); # in bytes

if (header_length != current_position){
    warning('After writing header info: header length does not match current file pointer')
}

writeBin(data, fid, size = 1, endian = "little", useBytes = modeUseBytes);
statusok <- close(fid);

return(TRUE)
}
