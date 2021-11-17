Place the mzXML file you are analyzing in this folder.

Name it toy.mzXML.



# toy.mzXML

```
cp /home/jia/Documents/Code/waterlooms/dia_benchmark/data/PXD005573/toy.mzXML /home/jia/Documents/Code/waterlooms/dia_data_reading/src/main/resources/toy.mzXML
```

toy.mzXML is a renamed Fig1_MP-DIA-120min-120kMS1-10W30k-14dppp_MHRM_R01.mzXML from PXD005573.
toy.mzXML_feature_all_z is a result of running Xiangyuan's FeatureDetect on toy.mzXML



# Code Execution
```
java -jar target/dia_data_reading-1.0-SNAPSHOT-jar-with-dependencies.jar
java -jar target/dia_data_reading-1.0-SNAPSHOT-jar-with-dependencies.jar -detectionParams src/main/resources/featuredetect.params -selectionParams src/main/resources/featureselect.params -mzXML src/main/resources/toy.mzXML
```

# TODOs
- Make it fail if param files are not passed in
- it will pass if the param files are not passed in and exit silently



# Notes
averagine_20160512.txt was sourced from [Online](https://zenodo.org/record/2652602).

Added the line to averagine file
"""
222.18	 1.000	 0.123
"""
- Probably doesn't make sense to have a precursor so small


toy_PeakCluster_xiangyuan.csv mimics toy_PeakCluster.csv 
	Generated from toy.mzXML which is a subset of dataset PXD005573
	Exported from running Xiangyuan's program on toy.mzXML

toy_peakCluster.csv was copied from /home/jia/Documents/Code/waterlooms/dia_benchmark/data/toy_data

