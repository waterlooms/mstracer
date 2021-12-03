package edu.uw.waterlooms.service;

/**
 * This class is utilized to hold the constants & parameters for WaterlooMS. It has built-in methods
 * for parsing a JSON-encoded string that is typically passed in from DIA-Webapp. Additionally,
 * there are constant definitions for feature-detection/feature-selection parameters.
 */
public class ParameterService {
  // Parameters utilized for FeatureDetect
  private static double C12_MASS = 12;
  private static double C13_MASS = 13.003355;
  private static double PPM = 0.00001;
  private static double MZ_TOLERANCE_PPM = 0.00001;
  private static int MIN_CHECK = 2;
  private static int MASS_STEP = 10;
  private static double MIN_REL_HEIGHT = 0.01;
  private static int INVALID_VAL = -1;
  private static double PERCENTAGE_OF_MAX = 0.5;
  private static int INTENSITY_THRESHOLD = 1;
  private static int Z_RANGE = 5;
  private static int WINDOW_SIZE = 4;
  private static double INTENSITY_BOUND = 1.5;

  // Parameters utilized for FeatureSelect
  private static final int ID_INDEX = 0;
  private static final int MZ_INDEX = 1;
  private static final int RT_INDEX = 2;
  private static final int Z_INDEX = 3;
  private static final int ISONUM_INDEX = 4;
  private static final int INT_SHAPE_INDEX = 5;
  private static final int ISO_DISTR_INDEX = 6;
  private static final int INTENSITY_AREA_PERCENTAGE_INDEX = 7;
  private static final int RT_START_INDEX = 8;
  private static final int RT_END_INDEX = 9;
  private static final int SCAN_SUM_INDEX = 10;
  private static final int PEAKS_SUM_INDEX = 11;
  private static final int PEAKS_AREA_INDEX = 12;
  private static final int SVR_INDEX = 13;
  private static final int QUALITY_INDEX = 14;
  private static double MZ_ERROR = 0.000008;
  private static double RT_ERROR = 0.5;

  // Parameters Configured from DIA-Webapp
  private static int MS1_PPM = 30;
  private static int MS2_PPM = 10;
  private static double RT_TOL_MIN = 0.16;
  private static int MIN_PPT = 5;
  private static int MIN_PEPTIDE_LEN = 7;
  private static int MAX_PRECURSOR_Z = 2;

  public ParameterService(){
  }

  public static double getC12Mass() {
    return C12_MASS;
  }

  public static void setC12Mass(int c12Mass) {
    C12_MASS = c12Mass;
  }

  public static double getC13Mass() {
    return C13_MASS;
  }

  public static void setC13Mass(double c13Mass) {
    C13_MASS = c13Mass;
  }

  public static double getPPM() {
    return PPM;
  }

  public static void setPPM(double PPM) {
    ParameterService.PPM = PPM;
  }

  public static double getMzTolerancePpm() {
    return MZ_TOLERANCE_PPM;
  }

  public static void setMzTolerancePpm(double mzTolerancePpm) {
    MZ_TOLERANCE_PPM = mzTolerancePpm;
  }

  public static int getMinCheck() {
    return MIN_CHECK;
  }

  public static void setMinCheck(int minCheck) {
    MIN_CHECK = minCheck;
  }

  public static int getMassStep() {
    return MASS_STEP;
  }

  public static void setMassStep(int massStep) {
    MASS_STEP = massStep;
  }

  public static double getMinRelHeight() {
    return MIN_REL_HEIGHT;
  }

  public static void setMinRelHeight(double minRelHeight) {
    MIN_REL_HEIGHT = minRelHeight;
  }

  public static int getInvalidVal() {
    return INVALID_VAL;
  }

  public static void setInvalidVal(int invalidVal) {
    INVALID_VAL = invalidVal;
  }

  public static double getPercentageOfMax() {
    return PERCENTAGE_OF_MAX;
  }

  public static void setPercentageOfMax(double percentageOfMax) {
    PERCENTAGE_OF_MAX = percentageOfMax;
  }

  public static int getIntensityThreshold() {
    return INTENSITY_THRESHOLD;
  }

  public static void setIntensityThreshold(int intensityThreshold) {
    INTENSITY_THRESHOLD = intensityThreshold;
  }

  public static int getzRange() {
    return Z_RANGE;
  }

  public static void setzRange(int zRange) {
    Z_RANGE = zRange;
  }

  public static int getWindowSize() {
    return WINDOW_SIZE;
  }

  public static void setWindowSize(int windowSize) {
    WINDOW_SIZE = windowSize;
  }

  public static double getIntensityBound() {
    return INTENSITY_BOUND;
  }

  public static void setIntensityBound(double intensityBound) {
    INTENSITY_BOUND = intensityBound;
  }

  public static int getIdIndex() {
    return ID_INDEX;
  }

  public static int getMzIndex() {
    return MZ_INDEX;
  }

  public static int getRtIndex() {
    return RT_INDEX;
  }

  public static int getzIndex() {
    return Z_INDEX;
  }

  public static int getIsonumIndex() {
    return ISONUM_INDEX;
  }

  public static int getIntShapeIndex() {
    return INT_SHAPE_INDEX;
  }

  public static int getIsoDistrIndex() {
    return ISO_DISTR_INDEX;
  }

  public static int getIntensityAreaPercentageIndex() {
    return INTENSITY_AREA_PERCENTAGE_INDEX;
  }

  public static int getRtStartIndex() {
    return RT_START_INDEX;
  }

  public static int getRtEndIndex() {
    return RT_END_INDEX;
  }

  public static int getScanSumIndex() { return SCAN_SUM_INDEX; }

  public static int getPeaksSumIndex() {
    return PEAKS_SUM_INDEX;
  }

  public static int getPeaksAreaIndex() {
    return PEAKS_AREA_INDEX;
  }

  public static int getSvrIndex() {
    return SVR_INDEX;
  }

  public static int getQualityIndex() {
    return QUALITY_INDEX;
  }

  public static double getMzError() {
    return MZ_ERROR;
  }

  public static void setMzError(double mzError) {
    MZ_ERROR = mzError;
  }

  public static double getRtError() {
    return RT_ERROR;
  }

  public static void setRtError(double rtError) {
    RT_ERROR = rtError;
  }

  public static int getMs1Ppm() {
    return MS1_PPM;
  }

  public static void setMs1Ppm(int ms1Ppm) {
    MS1_PPM = ms1Ppm;
  }

  public static int getMs2Ppm() {
    return MS2_PPM;
  }

  public static void setMs2Ppm(int ms2Ppm) {
    MS2_PPM = ms2Ppm;
  }

  public static double getRtTolMin() {
    return RT_TOL_MIN;
  }

  public static void setRtTolMin(double rtTolMin) {
    RT_TOL_MIN = rtTolMin;
  }

  public static int getMinPpt() {
    return MIN_PPT;
  }

  public static void setMinPpt(int minPpt) {
    MIN_PPT = minPpt;
  }

  public static int getMinPeptideLen() {
    return MIN_PEPTIDE_LEN;
  }

  public static void setMinPeptideLen(int minPeptideLen) {
    MIN_PEPTIDE_LEN = minPeptideLen;
  }

  public static int getMaxPrecursorZ() {
    return MAX_PRECURSOR_Z;
  }

  public static void setMaxPrecursorZ(int maxPrecursorZ) {
    MAX_PRECURSOR_Z = maxPrecursorZ;
  }

}
