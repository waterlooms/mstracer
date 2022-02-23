package edu.uw.waterlooms.ms1;

import edu.uw.waterlooms.entity.Pair;
import edu.uw.waterlooms.entity.XIC;
import edu.uw.waterlooms.msutil.OpenMzxml;
import java.io.*;
import java.util.*;

import edu.uw.waterlooms.service.ParameterService;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.spectrum.ISpectrum;
import umich.ms.fileio.exceptions.FileParsingException;

import javax.crypto.spec.PSource;

public class FeatureDetect {
  public enum DetectionType {
    MS1, MS2
  }
  private Set<Map.Entry<Integer, IScan>> scanEntries;
  private Set<Map.Entry<Integer, IScan>> scanEntries2;

  // MS1 params
  private double C12;
  private double C13;
  private double ppm;
  private double mzSearchRangePPM; // max tolerant m/z difference as to the same m/z
  private int MIN_CHECK; // number of consequent isotopic features to be identified
  private int MassStep;
  private double MinRelHeight;
  private int InvalidVal;
  private double PercentageOfMax;
  private double intensityThreshold;
  private int z_range;
  private int window_size;
  private int windowHi;
  private double IntensityBound;
  // MS1 features
  private double[][] intensityMap;
  private double[][] mzMap;
  private List<Integer>[] localpepMax;
  private int IDX;
  private double[] RT; // the array of real RT
  private double mzLow = Double.MAX_VALUE; // the smallest mz overall
  private double mzHi = Double.MIN_VALUE; // the largest mz overall
  private List<Double> scoreIntShape;
  private List<Double> scoreIsoDistr;
  private List<List<Pair<Double, Double>>> scorePepAll;
  private List<List<Pair<Integer, Integer>>>  scorePepAllPos;
  private List<Integer> scoreZVal;
  private List<Integer> scoreScanNum;
  private List<List<Double>> intCluster;
  private List<Double> intSum;
  private List<Double> intensityPercentage;
  private List<List<XIC>> trailsAll;

  /**
   * Read in all mzXML file information
   *
   * @param of new mzXML file opened by OpenMzxml
   */
  public FeatureDetect(OpenMzxml of, DetectionType detectionType) {
    switch (detectionType) {
      case MS1:
        setParameters();
        scanEntries = of.scanEntries;
        IDX = scanEntries.size();
        intensityMap = new double[IDX][];
        mzMap = new double[IDX][];
        RT = new double[IDX];
        scanEntries2 = of.scanEntries2;
        break;
      case MS2:
        scanEntries2 = of.scanEntries2;
    }
  }
  /** Helper function for constructor. Sets parameter values defined above. */
  private void setParameters() {
    // MS1 parameter
    C12 = ParameterService.getC12Mass();
    C13 = ParameterService.getC13Mass();
    ppm = ParameterService.getPPM();
    mzSearchRangePPM = ParameterService.getMzTolerancePpm();
    MIN_CHECK = ParameterService.getMinCheck();
    MassStep = ParameterService.getMassStep();
    MinRelHeight = ParameterService.getMinRelHeight();
    InvalidVal = ParameterService.getInvalidVal();
    PercentageOfMax = ParameterService.getPercentageOfMax();
    intensityThreshold = ParameterService.getIntensityThreshold();
    z_range = ParameterService.getzRange();
    window_size = ParameterService.getWindowSize();
    IntensityBound = ParameterService.getIntensityBound();
  }
  /**
   * Call this function to detect features from LC-MS, write to a file, followed by machine learning
   * scoring methods
   */
  public void detectFeatures(String oldFilePath) {
    String filepath = oldFilePath.replaceFirst("[.][^.]+$", "");
    init_MS1();
    findLocalMax();
//    System.out.println("findLocalMax completed");
    searchIsotope(window_size, z_range);
//    System.out.println("searchIsotope completed");
    isoDistributionScore();
//    System.out.println("isoDistributionScore completed");
    detectTrailofPrecursor();
//    System.out.println("detectTrailofPrecursor completed");
    try {
      writeFeatures(
              filepath
                      + "_featureAllZ.tsv"); // since a feature is searched with all possible charge states
    } catch (IOException e) {
      System.out.println("Write to file failed.");
    }
  }

  /** Read in for each Spectrum, including intensities of recorded m/z values */
  /** Read in for each Spectrum, including intensities of recorded m/z values */
  private void init_MS1() {
    int count = 0;
    double rtHi = 0;
    for (Map.Entry<Integer, IScan> scanEntry : scanEntries) {
      IScan scan = scanEntry.getValue();
      ISpectrum spectrum = null;
      try {
        spectrum = scan.fetchSpectrum();
      } catch (FileParsingException e) {
        e.printStackTrace();
      }
      double[] mz = spectrum.getMZs();
      double[] intensities = spectrum.getIntensities();
      if (mz[mz.length - 1] > mzHi) {
        mzHi = mz[mz.length - 1];
      }
      if (mz[0] < mzLow) {
        mzLow = mz[0];
      }

      int intAboveThresholdNum = 0;
      for (double intensity: intensities) {
        if (intensity > intensityThreshold) {
          intAboveThresholdNum++;
        }
      }

      double[] mzsNew = new double[intAboveThresholdNum];
      double[] intsNew = new double[intAboveThresholdNum];
      int i = 0;
      for (int j = 0; j < intensities.length; j++) {
        double intensity = intensities[j];
        if (intensity > intensityThreshold) {
          mzsNew[i] = mz[j];
          intsNew[i] = intensities[j];
          i++;
        }
      }
      RT[count] = scan.getRt();
      intensityMap[count] = intsNew;
      mzMap[count] = mzsNew;
      count++;

      if (rtHi < scan.getRt()) {
        rtHi = scan.getRt();
      }
    }
    // Choose window size; must be in 3-15 scans
    int windowCalculated = (int) Math.ceil(rtHi/15); // scaling function
    if (windowCalculated > window_size) {
      window_size = windowCalculated;
    }
    for (Map.Entry<Integer, IScan> scanEntry2 : scanEntries2) {
      IScan scan = scanEntry2.getValue();
      Object scanPrecursorCharge = scan.getPrecursor().getCharge();
      if (scanPrecursorCharge != null && z_range < (int)scanPrecursorCharge) {
        z_range = scan.getPrecursor().getCharge();
      }else {
        z_range=8;
      }
    }
//    System.out.println("window_size: " + window_size);
    System.out.println("Charge state searched: 1 ... " + z_range);
  }

  /** Initial detection on local maximal intensities on features */
  private void findLocalMax(){
    // Copy of intensityMap
    // Records if a m/z is reached. If reached, set to -1
    double[][] intMap = new double[IDX][];
    for (int i = 0; i < IDX; i++) {
      double[] copy = new double[intensityMap[i].length];
      for (int j = 0; j < intensityMap[i].length; j++) {
        copy[j] = intensityMap[i][j];
      }
      intMap[i] = copy;
    }

    localpepMax = new List[IDX];
    for(int i = 0; i < IDX; i++){
      localpepMax[i] = new ArrayList<>(); // initiate local max map
    }

    // Search max near one m/z horizontally
//    System.out.println("idx = " + IDX);
    for(int width = 0; width < IDX; width++){
      for(int height = 0; height < intMap[width].length; height++){
        if (intMap[width][height] == InvalidVal) {
          continue;
        }
        int count = 1;
        int col_for_max = width; // Column for max
        int pos_for_max = height; // Position for max
        int prev_col_for_max = InvalidVal; // Column for max
        int prev_pos_for_max = InvalidVal; // Position for max
        double pit_intval = InvalidVal;
        int last_col = width;
        int last_pos = height;
        double value_int_max = intMap[col_for_max][pos_for_max]; // The max value of intensity of one peptide
        double prev_value_int_max = InvalidVal;
        double prev_val = intMap[last_col][last_pos];
        intMap[col_for_max][pos_for_max] = InvalidVal; // Change the visited position to invalid

        // Search at the right column
        int rightPos = 0;
        boolean isDecrease = false;
        for (int rightCol = width + 1; rightCol < IDX; rightCol++){
          rightPos = searchRange(mzMap[rightCol], mzMap[last_col][last_pos] * (1-mzSearchRangePPM), mzMap[last_col][last_pos] * (1+mzSearchRangePPM));
          // multi-local-max
          //start
          if (rightPos == InvalidVal){
            if (count >= MIN_CHECK) {
              localpepMax[col_for_max].add(pos_for_max);
              if (prev_value_int_max!=InvalidVal && pit_intval < value_int_max * PercentageOfMax) {
                localpepMax[prev_col_for_max].add(prev_pos_for_max);
              }
            }
            count = 1;
            break;
          }
          else if (intMap[rightCol][rightPos] < prev_val) {
            if (value_int_max <= prev_value_int_max) {
              value_int_max = prev_value_int_max;
              col_for_max = prev_col_for_max;
              pos_for_max = prev_pos_for_max;
            }
            else if (prev_value_int_max!=InvalidVal && pit_intval < value_int_max * PercentageOfMax) { // the second max, value_int_max > prev_value_int_max
              localpepMax[prev_col_for_max].add(prev_pos_for_max);
            }
            prev_value_int_max = InvalidVal; // turned off when not necessary
            isDecrease = true;
          }
          else if (isDecrease) { // cur_val >= prev_val
            if (prev_val < value_int_max * PercentageOfMax) { // prev_val is lowest intensity point
              localpepMax[col_for_max].add(pos_for_max);
            } else {
              prev_col_for_max = col_for_max;
              prev_pos_for_max = pos_for_max;
              prev_value_int_max = value_int_max;
            }
            value_int_max = intMap[rightCol][rightPos];
            col_for_max = rightCol;
            pos_for_max = rightPos;
            pit_intval = prev_val;
            isDecrease = false;
          }
          if (intMap[rightCol][rightPos] > value_int_max) {
            value_int_max = intMap[rightCol][rightPos];
            col_for_max = rightCol;
            pos_for_max = rightPos;
          }
          // update
          last_col = rightCol;
          last_pos = rightPos;
          prev_val = intMap[rightCol][rightPos];
          count++;
          //end

          // turn the visited to invalid
          intMap[rightCol][rightPos] = InvalidVal;

        }
        // add last local max position;
        if (rightPos == IDX && count >= MIN_CHECK) {
          localpepMax[col_for_max].add(pos_for_max);
        }
      }
    }
    // sort
    for (int i = 0; i < IDX; i++) {
      Collections.sort(localpepMax[i]);
    }
  }
  /**
   * Search for isotopes of all charge states based on local maximal intensities; calculate
   * intensity shape function to measure feature quality, record only when isotopes >= 2
   *
   * @param window_size the max distance to search isotopes
   * @param num_of_z z values: 1 ... z
   */
  private void searchIsotope(int window_size, int num_of_z) {
    scoreIntShape = new ArrayList<>();
    scorePepAll = new ArrayList<>();
    scorePepAllPos = new ArrayList<>();
    scoreZVal = new ArrayList<>();
    intCluster = new ArrayList<>(); // Cluster of intensities corresponding to feature cluster
    intSum = new ArrayList<>();
    scoreScanNum = new ArrayList<>();
    intensityPercentage = new ArrayList<>();
    int groupID = 0; // Corresponding to scorePepAll index

    // Search for feature all z
    for (int z = num_of_z; z >= 1; z--) {
      List<Integer>[] maxMap = new List[IDX];
      for (int i = 0; i < IDX; i++) {
        List<Integer> copy = new ArrayList<>(localpepMax[i]);
        maxMap[i] = copy;
      }
      for (int i = 0; i < IDX - window_size + 1; i++) {
        for (int k = 0; k < maxMap[i].size(); k++) {
          if (maxMap[i].get(k) == InvalidVal) {
            continue;
          }
          double cur_mz = mzMap[i][maxMap[i].get(k)]; // M/Z value of current point
          double center_mz = cur_mz;
          List<Pair<Integer, Integer>> tempPos =
                  new ArrayList<>(); // All positions found in terms of localPepMax. {i, k}
          tempPos.add(new Pair<>(i, k));

          // Intensities of current peak and peaks within window size
          List<Double> curInt = new ArrayList<>();
          List<Double> maxInt = new ArrayList<>(); // max intensity for an isotope

          // Intensities for a whole cluster in increasing order of mz
          List<List<Double>> clusterInt = new ArrayList<>(); // intensities for a whole cluster
          curInt.add(intensityMap[i][maxMap[i].get(k)]);
          maxInt.add(intensityMap[i][maxMap[i].get(k)]);

          // Looking for intensity
          // Recorded as 0, if no reads
          for (int p = 1; p < window_size; p++) {
            int pos_search =
                    searchRange(
                            mzMap[i + p], cur_mz * (1 - mzSearchRangePPM), cur_mz * (1 + mzSearchRangePPM));
            if (pos_search != InvalidVal) {
              curInt.add(intensityMap[i + p][pos_search]);
            } else {
              curInt.add(0.0);
            }
          }
          clusterInt.add(curInt);

          double lowrange = cur_mz * (1 - ppm) - (C13 - C12) / z;
          double hirange = lowrange + 2 * cur_mz * ppm;
          int lowerPos = InvalidVal;
          while (true) {
            int j = 1;
            List<Double> newInt = new ArrayList<>();
            for (; j < window_size; j++) {
              lowerPos = searchPos(maxMap[i + j], mzMap[i + j], lowrange, hirange);
              if (lowerPos != InvalidVal) {
                cur_mz = mzMap[i + j][maxMap[i + j].get(lowerPos)];
                lowrange = cur_mz * (1 - ppm) - (C13 - C12) / z;
                hirange = lowrange + 2 * cur_mz * ppm;
                tempPos.add(0, new Pair<>(i + j, lowerPos));
                break;
              }
            }
            // If no point found below, break
            // If found, record intensity for the feature
            if (lowerPos == InvalidVal) {
              break;
            } else {
              for (int u = 0; u < window_size; u++) {
                if (u == j) {
                  maxInt.add(0, intensityMap[i + j][maxMap[i + j].get(lowerPos)]);
                  newInt.add(intensityMap[i + j][maxMap[i + j].get(lowerPos)]);
                  continue;
                }
                int pos_search =
                        searchRange(
                                mzMap[i+u], cur_mz * (1-mzSearchRangePPM), cur_mz * (1+mzSearchRangePPM));
                if (pos_search != -1) {
                  newInt.add(intensityMap[i+u][pos_search]);
                } else  {
                  newInt.add(0.00);
                }
              }
            }
            clusterInt.add(0, newInt);
          }

          // Check above the center peptide for isotopes
          lowrange = center_mz * (1 - ppm) + (C13 - C12) / z;
          hirange = lowrange + 2 * center_mz * ppm;
          int higherPos = InvalidVal;
          while (true) {
            int j = 0;
            List<Double> newInt = new ArrayList<>();
            for (; j < window_size; j++) {
              higherPos = searchPos(maxMap[i + j], mzMap[i + j], lowrange, hirange);
              if (higherPos != InvalidVal) {
                cur_mz = mzMap[i + j][maxMap[i + j].get(higherPos)];
                lowrange = cur_mz * (1 - ppm) + (C13 - C12) / z;
                hirange = lowrange + 2 * center_mz * ppm;
                tempPos.add(new Pair<>(i + j, higherPos));
                break;
              }
            }
            // If no point found above, break
            // If found, record intensity for the feature
            if (higherPos == InvalidVal) {
              break;
            } else {
              for (int u = 0; u < window_size; u++) {
                if (u == j) {
                  maxInt.add(intensityMap[i + j][maxMap[i + j].get(higherPos)]);
                  newInt.add(intensityMap[i + j][maxMap[i + j].get(higherPos)]);
                  continue;
                }
                int pos_search =
                        searchRange(
                                mzMap[i + u], cur_mz * (1 - mzSearchRangePPM), cur_mz * (1 + mzSearchRangePPM));
                if (pos_search != InvalidVal) {
                  newInt.add(intensityMap[i + u][pos_search]);
                } else {
                  newInt.add(0.0);
                }
              }
            }
            clusterInt.add(newInt);
          }
          int num_pep = tempPos.size();

          // Check if there are >=MIN_CHECK isotopes and intensity of major isotope is >=intensityThreshold
          if (num_pep >= MIN_CHECK && maxInt.get(0) >= intensityThreshold) {
            double intShape = 0;

            // Calculate intensity shape score
            for (int t = 0; t < tempPos.size() - 1; t++) {
              intShape += addScore(clusterInt.get(t), clusterInt.get(t + 1));
            }
            intShape /= tempPos.size() - 1;
            // Calculate sum of intensity for all lines
            List<Double> intSumLine = new ArrayList<>();
            for (int t = 0; t < clusterInt.size(); t++) {
              double intSum = 0;
              for (int v = 0; v < clusterInt.get(t).size(); v++) {
                intSum += clusterInt.get(t).get(v);
              }
              intSumLine.add(intSum);
            }
            // Turn local max into invalid
            // Turn tempPos into real position
            List<Pair<Double, Double>> realPos = new ArrayList<>();
            List<Pair<Integer, Integer>> virtualPos = new ArrayList<>();
            for (int l = 0; l < tempPos.size(); l++) {
              // Turn local max into invalid
              maxMap[tempPos.get(l).getL()].set(tempPos.get(l).getR(), InvalidVal);
              // Turn tempPos into real position
              realPos.add(
                      new Pair<>(
                              RT[tempPos.get(l).getL()],
                              mzMap[tempPos.get(l).getL()][
                                      localpepMax[tempPos.get(l).getL()].get(tempPos.get(l).getR())]));
              virtualPos.add(
                      new Pair<>(
                              tempPos.get(l).getL(),
                              localpepMax[tempPos.get(l).getL()].get(tempPos.get(l).getR())));

            }

            // calculate intensity percentage of local area
            double lowMZ = realPos.get(0).getR() - IntensityBound;
            double hiMZ = realPos.get(realPos.size()-1).getR() + IntensityBound;
            double areaIntSum = 0;
            double isoIntSum = 0;
            for (int m = 0; m < window_size; m++) {
              areaIntSum += intensityColSum(mzMap[i+m], intensityMap[i+m], lowMZ, hiMZ);
            }
            for (int m = 0; m < intSumLine.size(); m++) {
              isoIntSum += intSumLine.get(m);
            }
            double isoIntPercentage = isoIntSum/areaIntSum;

            groupID++;
            scorePepAll.add(realPos);
            scorePepAllPos.add(virtualPos);
            scoreIntShape.add(intShape);
            scoreZVal.add(z);
            scoreScanNum.add(tempPos.size());
            intCluster.add(intSumLine);
            intSum.add(isoIntSum/window_size);
            intensityPercentage.add(isoIntPercentage);
          }
        }
      }
    }
//    System.out.println("All peptide found before deleting duplicates: " + groupID);
  }

  private void isoDistributionScore() {
    scoreIsoDistr = new ArrayList<>();
    for (int cur_pep = 0; cur_pep < intCluster.size(); cur_pep++) {
      // Calculate vector score of experimental and theoretical isotope intensity
      float m_over_z = scorePepAll.get(cur_pep).get(0).getR().floatValue();
      float m = m_over_z * scoreZVal.get(cur_pep);
      int maxMass = (int) (m + (C13 - C12) * scoreScanNum.get(cur_pep) / scoreZVal.get(cur_pep) + 1);
      TheoreticalIsotope iso =
              new TheoreticalIsotope(scoreScanNum.get(cur_pep), maxMass, MassStep, MinRelHeight);
      float[] iso_shape = iso.computeShape(m);
      int IsoNum = iso_shape.length;
      double xi_sqr_sum = 0;
      double yi_sqr_sum = 0;
      double[] xi = new double[IsoNum];

      // Calculate isotope distribution score
      for (int i = 0; i < IsoNum; i++) {
        yi_sqr_sum += Math.pow(iso_shape[i], 2);
        xi[i] = intCluster.get(cur_pep).get(i);
        xi_sqr_sum += Math.pow(xi[i], 2);
      }

      double xi_yi_dot_sum = 0;
      for (int i = 0; i < IsoNum; i++) {
        xi_yi_dot_sum += xi[i] * iso_shape[i];
      }
      double iso_score = xi_yi_dot_sum / (Math.sqrt(xi_sqr_sum) * Math.sqrt(yi_sqr_sum));
      scoreIsoDistr.add(iso_score);
    }
  }

  private void detectTrailofPrecursor(){
    trailsAll = new ArrayList<>();

    for (List<Pair<Integer, Integer>> cluster: scorePepAllPos) {
      List<XIC> trails = new ArrayList<>();

      for (Pair<Integer, Integer> isotope: cluster) {
        ArrayList<Double>
                ints = new ArrayList<>(),
                mzs = new ArrayList<>(),
                rts = new ArrayList<>();

        ints.add(intensityMap[isotope.getL()][isotope.getR()]);
        mzs.add(mzMap[isotope.getL()][isotope.getR()]);
        rts.add(RT[isotope.getL()]);

        double lowMZ = mzMap[isotope.getL()][isotope.getR()] * (1-mzSearchRangePPM);
        double hiMZ = mzMap[isotope.getL()][isotope.getR()] * (1+mzSearchRangePPM);
        int leftScan = isotope.getL() - 1;
        int rightScan = isotope.getL() + 1;
        double leftMaxIntensity = intensityMap[isotope.getL()][isotope.getR()];
        double rightMaxIntensity = intensityMap[isotope.getL()][isotope.getR()];
        for (;leftScan >= 0 ; leftScan--){
          int leftPos = searchRange(mzMap[leftScan], lowMZ, hiMZ);
          if (leftPos == InvalidVal || intensityMap[leftScan][leftPos] > leftMaxIntensity) {
            break;
          }
          lowMZ = mzMap[leftScan][leftPos] * (1-mzSearchRangePPM);
          hiMZ = mzMap[leftScan][leftPos] * (1+mzSearchRangePPM);
          leftMaxIntensity = intensityMap[leftScan][leftPos];

          ints.add(0, intensityMap[leftScan][leftPos]);
          mzs.add(0, mzMap[leftScan][leftPos]);
          rts.add(0, RT[leftScan]);

        }
        for (;rightScan < IDX ; rightScan++){
          int rightPos = searchRange(mzMap[rightScan], lowMZ, hiMZ);
          if (rightPos == InvalidVal || intensityMap[rightScan][rightPos] > rightMaxIntensity) {
            break;
          }
          lowMZ = mzMap[rightScan][rightPos] * (1-mzSearchRangePPM);
          hiMZ = mzMap[rightScan][rightPos] * (1+mzSearchRangePPM);
          rightMaxIntensity = intensityMap[rightScan][rightPos];

          ints.add(intensityMap[rightScan][rightPos]);
          mzs.add(mzMap[rightScan][rightPos]);
          rts.add( RT[rightScan]);
        }
        if (leftScan < 0) { leftScan = 0; }
        if (rightScan == IDX) { rightScan = IDX - 1; }

        XIC trail = new XIC(ints, mzs, rts);
        trails.add(trail);
      }
      trailsAll.add(trails);
    }
  }

  public int searchRange(double[] array, double lowrange, double highrange) {
    int pos = InvalidVal;
    for (int i = 0; i < array.length; i++) {
      if (array[i] >= lowrange && array[i] <= highrange) {
        pos = i;
        break;
      }
      if (array[i] > highrange) {
        break;
      }
    }
    return pos;
  }
  // Given array of localmax positions, check in range
  // Return position at posArray
  public int searchPos(
          List<Integer> posArray, double[] mzArray, double lowrange, double highrange) {
    int pos = InvalidVal;
    for (int i = 0; i < posArray.size(); i++) {
      if (posArray.get(i) == InvalidVal) continue;
      if (mzArray[posArray.get(i)] >= lowrange && mzArray[posArray.get(i)] <= highrange) {
        pos = i;
        break;
      }
      if (mzArray[posArray.get(i)] > highrange) {
        break;
      }
    }
    return pos;
  }

  double addScore(List<Double> a, List<Double> b) {
    double dot_sum = 0;
    double a_sqr_sum = 0;
    double b_sqr_sum = 0;
    for (int i = 0; i < a.size(); i++) {
      dot_sum += a.get(i) * b.get(i);
      a_sqr_sum += a.get(i) * a.get(i);
      b_sqr_sum += b.get(i) * b.get(i);
    }
    return dot_sum / (Math.sqrt(a_sqr_sum) * Math.sqrt(b_sqr_sum));
  }

  public double intensityColSum(double[] array, double[] intensity, double lowrange, double highrange){
    double intSum = 0;
    for (int i = 0; i < array.length; i++) {
      if (array[i] >= lowrange && array[i] <= highrange) {
        intSum += intensity[i];
      }
      if (array[i] > highrange) {
        break;
      }
    }
    return intSum;
  }

  void writeFeatures(String filename) throws IOException {
    File file = new File(filename);
    try {
      // create FileWriter object with file as parameter
      FileWriter outputfile = new FileWriter(file);
      // PrintWriter
      PrintWriter printWriter = new PrintWriter(outputfile);
      // add header
      String header =
              "id"
                      + '\t'
                      + "mz"
                      + '\t'
                      + "rt"
                      + '\t'
                      + "z"
                      + '\t'
                      + "isotope_num"
                      + '\t'
                      + "intensity_shape_score"
                      + '\t'
                      + "isotope_distribution_score"
                      + '\t'
                      + "intensity_area_percentage"
                      + '\t'
                      + "scan_num"
                      + '\t'
                      + "quantification_peak_sum"
                      + '\t'
                      + "quantification_peak_area"
                      + '\t'
                      + "svr_score"
                      + '\t'
                      + "quality_score"
                      + '\t'
                      + "mzs,rts,ints"
                      + '\n';
      printWriter.print(header);

      for (int i = 0; i < scorePepAll.size(); i++) {
        int size = scorePepAll.get(i).size();
        for (int j = 0; j < size; j++) {
          int id = i + 1;
          XIC trail = trailsAll.get(i).get(j);
          int scannum = trail.getScanNum();
          double peak_sum = trail.getPeakSum();
          double peak_area = trail.getPeakArea();

          ArrayList<Double> ints = trail.getIntensities();
          ArrayList<Double> mzs = trail.getMassChargeRatios();
          ArrayList<Double> rts = trail.getRetentionTimes();

          String mzsStr =  mzs.get(0).toString();
          String rtsStr =  rts.get(0).toString();
          String intsStr =  ints.get(0).toString();
          for (int k = 1; k < scannum; k++) {
            mzsStr += '\t' + mzs.get(k).toString();
            rtsStr += '\t' + rts.get(k).toString();
            intsStr += '\t' + ints.get(k).toString();
          }

          String data = "";
          data += id;
          data += '\t' + scorePepAll.get(i).get(j).getR().toString();
          data += '\t' + scorePepAll.get(i).get(j).getL().toString();
          data += '\t' + scoreZVal.get(i).toString();
          data += '\t' + scoreScanNum.get(i).toString();
          data += '\t' + scoreIntShape.get(i).toString();
          data += '\t' + scoreIsoDistr.get(i).toString();
          data += '\t' + intensityPercentage.get(i).toString();
          data += '\t' + Integer.toString(scannum);
          data += '\t' + Double.toString(peak_sum);
          data += '\t' + Double.toString(peak_area);
          data += '\t' + "0"; // placeholder for svr_score
          data += '\t' + "0"; // placeholder for quality score
          data += '\t' + mzsStr;
          data += '\t' + rtsStr;
          data += '\t' + intsStr;
          data += '\n';
          printWriter.print(data);
        }
      }
      // closing writer connection
      printWriter.close();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
}
