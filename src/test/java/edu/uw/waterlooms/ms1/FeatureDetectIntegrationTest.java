package edu.uw.waterlooms.ms1;
import edu.uw.waterlooms.msutil.OpenMzxml;
import org.junit.Test;

import java.io.FileReader;
import java.io.IOException;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Properties;

public class FeatureDetectIntegrationTest {

    @Test
    public void matchMS1TrailsWithFeatureDetectionResult() throws IOException {
        /* Given */
        String dataDirectory = FileSystems.getDefault().getPath("").toAbsolutePath().toString() + "/data/";
        String featureDetectParamFile = dataDirectory + "/featuredetect.params";
        String rawFileName = "toy.mzXML";
        String mzxmlFile = dataDirectory + rawFileName;
        int rtNextPeakTolSec = 5;  // 5sec
        double mzTolerancePPMStrict = 10e-6; // 10ppm
        int MIN_PEAKNUM_Strict = 2;
        double rtMaxRangeSec = 30; // 30sec
        double intensityNextPeakPercentageTol = 0.9; // 90%

        /* When */
        OpenMzxml openMzxml = new OpenMzxml(mzxmlFile);
        FeatureDetect featureDetect = new FeatureDetect(openMzxml, FeatureDetect.DetectionType.MS1);
//        MS1TrailDetect ms1TrailDetect = new MS1TrailDetect();

        /* Then */
//        featureDetect.detectMS1Trails(mzxmlFile, rtNextPeakTolSec, mzTolerancePPMStrict, MIN_PEAKNUM_Strict, rtMaxRangeSec, intensityNextPeakPercentageTol);
//        ms1TrailDetect.matchFeatureToTrails(mzxmlFile, dataDirectory, rawFileName, featureDetect.getMS1Trails());
//        featureDetect.writeMS1TrailsData(dataDirectory + rawFileName + "_ms1_trails.tsv", featureDetect.getMS1Trails());
    }

    /**
     * Given an Absolute Path to a FeatureDetect Params File, return the Parameters in a Property
     * @param featureDetectParamFile String with absolute value of path
     * @return Properties of parameters file
     */
    private Properties parseFeatureDetectParamsFile(String featureDetectParamFile){
        Properties featureDetectParams = new Properties();
        try {
            FileReader fin = new FileReader(featureDetectParamFile);
            featureDetectParams.load(fin);
        } catch (IOException $e) {
            System.err.printf($e.getMessage());
            System.exit(1);
        }
        return featureDetectParams;
    }

}
