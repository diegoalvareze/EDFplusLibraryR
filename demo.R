source("EDFreadHeader.R")
source("EDFreadAnnotations.R")
source("annotations2edfplus.R")

edfPlusInputPath <-  "demoData/demoEDFplus.edf";
edfplusOutputPath <-  "demoOutputUsingR.EDF";

header <- EDFreadHeader(edfPlusInputPath);

#annotations <- EDFreadAnnotations(edfPlusInputPath, FALSE);
annotations <- EDFreadAnnotations(edfPlusInputPath, TRUE); # Skips EDF block's 'Time-keeping' annotations

# Writting annotations back
annotations2EDFplus(annotations, header$local_patient_identification, header$local_recording_identification, 
                    header$startdate_recording, header$starttime_recording, edfplusOutputPath);

