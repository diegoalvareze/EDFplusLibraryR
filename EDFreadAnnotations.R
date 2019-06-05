# Reads annotations of an EDF+ file into a list structure
# Time-keeping annotations can be skipped by setting the parameter
# doSkipTimeKeep to TRUE
#
#   annotations <- EDFreadAnnotations(filename, doSkipTimeKeep)
#
# More info, see EDF+ specification at http://www.edfplus.info/
#
# See also: EDFreadHeader
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
EDFreadAnnotations = function (filename, doSkipTimeKeep = FALSE) {

annotations <- NULL

general_header_size <- 256 #bytes
one_signal_header_size <- 256 #bytes

header <- EDFreadHeader(filename)

if (length(header) == 0)
{
    stop('Unable to read file header')
    return(-1)
}

fid <- file(filename, "rb") #0pen to read in binary mode

if (fid == -1){
    stop('Error opening the file for writing')
    return(-1)
}

# Set the pointer to the annotations signal
signalIndex <- header$num_signals;

# Positioning on the begining of the selected file
header_length <- general_header_size + header$num_signals * one_signal_header_size;
seek(fid, where = header_length + header$signals_signalOffset[signalIndex], rw = "read", origin = "start");

record_size <- header$signals_num_samples_datarecord[signalIndex]; # In samples (1 sample = 2 bytes)
bytes_full_data_record = 2 * sum(header$signals_num_samples_datarecord);

# We will read all the annotations
numRecords2read <- header$num_data_records;

# Data reading
data <- as.integer(readBin(fid, raw(), 2*record_size, size = 1, signed = FALSE, endian = "little"));
for (i in 1:header$num_data_records-1) 
{
    seek(fid, where = bytes_full_data_record - 2*record_size, rw = "read", origin = "current");
    datatmp <- as.integer(readBin(fid, raw(), 2*record_size, size = 1, signed = FALSE, endian = "little"));
    data <- c(data, datatmp);
}
close(fid)

# ASCII '20' separates each annotation after the time stamp (onset+duration)
tmp <- which(data == 20);
# Detection of 'end of TAL': ASCII '20' followed by ASCII '0'. Each TAL may contain several annotations, and the last one finishes with this combination
# Note: Each TAL shares time_stamp and duration. Several TALS can be
#       included in the same block, but the first one contains an empty
#       annotation that indicates the time_stamp of the block itself.
numberTALends <- tmp[data[tmp+1] == 0];

countAnn <- 0;

for (k in 1:length(numberTALends)) 
{
  if (k == 1)
    tal <- data[1:numberTALends[1]]
  else
    tal <- data[(numberTALends[k-1]+2):numberTALends[k]];
  
  if (!identical(intToUtf8(tal[1]),'+') && !identical(intToUtf8(tal[1]),'-'))
  {
    # Remove previous characters til + or - (possibly by block transition, thus zeros at left are expected)
    strIdx = which(tal == 43); # Looking for '+' character
    if (length(strIdx) == 0)
      strIdx = which(tal == 43); # Looking for '-' character
    if (length(strIdx) == 0)
      stop('Incorrect start of TAL');
    tal = tal[strIdx:length(tal)];
  }
  
  # Find offset and duration
  sep = which(tal == 20);
  indxDur = which(tal == 21);
  if (length(indxDur) == 0)
  {
    dur = 0;
    offset = as.numeric(intToUtf8(tal[1:(sep[1]-1)]));
  }
  else
  {
    dur = as.numeric(intToUtf8(tal[(indxDur[1]+1):(sep[1]-1)]));
    offset = as.numeric(intToUtf8(tal[1:(indxDur[1]-1)]));
  }
  
  # Read annotations
  for (k1 in 1:(length(sep)-1))
  {
    idxS = sep[k1] + 1;
    idxE = sep[k1+1] - 1;
    
    if (idxS > idxE)
    {
      if (!doSkipTimeKeep)
      {
        countAnn <- countAnn + 1;
        annotations$label[countAnn] <- 'Time-keep';
        annotations$offset[countAnn] <- offset;
        annotations$duration[countAnn] <- dur;
      }
    }
    else
    {
      countAnn <- countAnn + 1;
      annotations$label[countAnn] <- intToUtf8(tal[idxS:idxE]);
      annotations$offset[countAnn] <- offset;
      annotations$duration[countAnn] <- dur;
    }
  }
}

return(annotations)
}