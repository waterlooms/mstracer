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
  private double IntensityThreshold;
  // MS1 features
  private double[][] intensityMap;
  private double[][] mzMap;
  private List<Pair<Integer, Integer>>[] localpepMax;
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
  private List<Double> intensityPercentage;
  private List<List<Integer>> posInMS1Trail;
  private ArrayList<XIC> MS1Trails;

  // TODO: refactor parameters
  private int SEC_IN_MIN = 60;

  private double rtNextPeakTolMin = 0.1;  // 0.1min
  private double mzTolerancePPM = 10e-6; // 10ppm
  private int MIN_PEAKNUM = 2;
  private double rtMaxRangeMin = 0.5; // 0.5min
  private double intensityNextPeakPercentageTol = 0.9; // 90%
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

        // TODO: Added by jia to prevent init() from crashing
        scanEntries2 = of.scanEntries2;
        break;
      case MS2:
        scanEntries2 = of.scanEntries2;
    }
  }
  /**
   * Call this function to detect features from LC-MS, write to a file, followed by machine learning
   * scoring methods
   */
  public void detectFeatures(String filepath) {
    init_MS1();
    searchInfoList(rtNextPeakTolMin, mzTolerancePPM, MIN_PEAKNUM, rtMaxRangeMin, intensityNextPeakPercentageTol);
    System.out.println("getLocalMaxFromMS1Trails completed");
    searchIsotope(window_size, z_range);
    System.out.println("searchIsotope completed");
    calculateIsoDistributionScore();
    System.out.println("calculateIsoDistributionScore completed");
    try {
      writeFeatures(
              filepath
                      + "_feature_all_z"); // since a feature is searched with all possible charge states
    } catch (IOException e) {
      System.out.println("Write to file failed.");
    }
  }

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
    System.out.println("window_size: " + window_size);
    System.out.println("z_range: " + z_range);
  }

  private void searchInfoList(double rtNextPeakTol, double mzTolerancePPM, int MIN_PEAKNUM, double rtMaxRangeTol, double intensityNextPeakPercentageTol){
    MS1Trails = new ArrayList<>();
    ArrayList<Pair<Integer, Integer>> resultIndex; // A list of <index, rt> wrt the XIC.
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
    for(int i = 0; i < IDX; i++) {
      if (i % 100 == 0) {
        System.out.println("At: " + i + "/" + IDX);
      }
      int IDXRight = i;
      for (int k = i + 1; k < IDX; k++) {
        if (RT[k] - RT[i] > rtNextPeakTol) {
          IDXRight = k - 1;
        }
      }
      for (int j = 0; j < mzMap[i].length; j++) {
        if (intMap[i][j] == InvalidVal) {
          continue;
        }
        double prev_mz = mzMap[i][j];
        double prev_rt = RT[i];
        double prev_intensity = intMap[i][j];
        boolean ifDecrease = false;
        double max_intensity = prev_intensity;
        double max_peak_rt = prev_rt;
        resultIndex = new ArrayList<>();
        resultIndex.add(new Pair<>(i, j));
        for (int k = i + 1; k < IDXRight; k++) {
          if (RT[k] > max_peak_rt + rtMaxRangeTol / 2) {
            break;
          }
          int pos = searchTarget(mzMap[k], intMap[k], prev_mz);
          if (pos == InvalidVal) {
            continue;
          }
          double mz2 = mzMap[k][pos];
          double rt2 = RT[k];
          double intensity2 = intMap[k][pos];
          if (intensity2 != InvalidVal
                  && mz2 >= prev_mz * (1 - mzTolerancePPM)
                  && mz2 <= prev_mz * (1 + mzTolerancePPM)) {
            if (intensity2 > prev_intensity && ifDecrease) {
              break;
            } else ifDecrease = intensity2 < prev_intensity * intensityNextPeakPercentageTol;
            prev_intensity = intensity2;
            resultIndex.add(new Pair<>(k, pos));
            prev_mz = mz2;
            prev_rt = rt2;
            if (intensity2 > max_intensity) {
              max_intensity = intensity2;
              max_peak_rt = rt2;
            }
          }
        }
        // Must contain at least MIN_PEAKNUM peaks
        if (resultIndex.size() < MIN_PEAKNUM) {
          continue;
        }
        // Add XIC trail
        ArrayList<Double> mzs = new ArrayList<>();
        ArrayList<Double> rts = new ArrayList<>();
        ArrayList<Double> ints = new ArrayList<>();
        ArrayList<Integer> cols = new ArrayList<>();
        ArrayList<Integer> rows = new ArrayList<>();
        for (Pair<Integer, Integer> item : resultIndex) { // Add to trails
          int col = item.getL();
          int row = item.getR();
          mzs.add(mzMap[col][row]);
          rts.add(RT[col]);
          cols.add(col);
          rows.add(row);
          ints.add(intMap[col][row]);
          intMap[col][row] = InvalidVal;
        }
        XIC trail = new XIC(ints, mzs, rts, cols, rows);
        MS1Trails.add(trail);

        // transfer to local max map
        int columnAtMaxIntensity = trail.getColumeInMapAtMaxIntensity();
        int rowAtMaxIntensity = trail.getRowInMapAtMaxIntensity();
        localpepMax[columnAtMaxIntensity].add(new Pair<>(rowAtMaxIntensity, MS1Trails.size() - 1));
      }
    }
    int count = 0;
    int total = 0;
    for (int i = 0; i < IDX; i++) {
      total += intMap[i].length;
      for (int j = 0; j < intMap[i].length; j++)
      if (intMap[i][j] == InvalidVal) {
        count++;
      }
    }
    double percentage = (double)count/total;
    System.out.println("Number of peaks covered: " + count + "/" + total);
    System.out.println("Percentage: " + percentage);
    System.out.println("Number of trails detected: " + MS1Trails.size());
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
    scoreScanNum = new ArrayList<>();
    intensityPercentage = new ArrayList<>();
    posInMS1Trail = new ArrayList<>();
    int groupID = 0; // Corresponding to scorePepAll index

    // Search for feature all z
    for (int z = num_of_z; z >= 1; z--) {
      List<Pair<Integer, Integer>>[] maxMap = new List[IDX];
      for (int i = 0; i < IDX; i++) {
        List<Pair<Integer, Integer>> copy = new ArrayList<>(localpepMax[i]);
        maxMap[i] = copy;
      }
      for (int i = 0; i < IDX - window_size + 1; i++) {
        for (int k = 0; k < maxMap[i].size(); k++) {
          if (maxMap[i].get(k).getL() == InvalidVal) {
            continue;
          }
          double cur_mz = mzMap[i][maxMap[i].get(k).getL()]; // M/Z value of current point
          double center_mz = cur_mz;
          List<Pair<Integer, Integer>> tempPos =
              new ArrayList<>(); // All positions found in terms of localPepMax. {i, k}
          tempPos.add(new Pair<>(i, k));

          // Intensities of current peak and peaks within window size
          List<Double> curInt = new ArrayList<>();
          List<Double> maxInt = new ArrayList<>(); // max intensity for an isotope

          // Intensities for a whole cluster in increasing order of mz
          List<List<Double>> clusterInt = new ArrayList<>(); // intensities for a whole cluster
          curInt.add(intensityMap[i][maxMap[i].get(k).getL()]);
          maxInt.add(intensityMap[i][maxMap[i].get(k).getL()]);

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
                cur_mz = mzMap[i + j][maxMap[i + j].get(lowerPos).getL()];
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
                  maxInt.add(0, intensityMap[i + j][maxMap[i + j].get(lowerPos).getL()]);
                  newInt.add(intensityMap[i + j][maxMap[i + j].get(lowerPos).getL()]);
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
                cur_mz = mzMap[i + j][maxMap[i + j].get(higherPos).getL()];
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
                  maxInt.add(intensityMap[i + j][maxMap[i + j].get(higherPos).getL()]);
                  newInt.add(intensityMap[i + j][maxMap[i + j].get(higherPos).getL()]);
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

          // Check if there are >=MIN_CHECK isotopes
          if (num_pep >= MIN_CHECK) {
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
            List<Integer> ms1trailPos = new ArrayList<>();
            for (int l = 0; l < tempPos.size(); l++) {
              // Turn local max into invalid
              maxMap[tempPos.get(l).getL()].set(tempPos.get(l).getR(), new Pair<>(InvalidVal, InvalidVal));
              // Turn tempPos into real position
              realPos.add(
                  new Pair<>(
                      RT[tempPos.get(l).getL()],
                      mzMap[tempPos.get(l).getL()][
                          localpepMax[tempPos.get(l).getL()].get(tempPos.get(l).getR()).getL()]));
              virtualPos.add(
                  new Pair<>(
                      tempPos.get(l).getL(),
                      localpepMax[tempPos.get(l).getL()].get(tempPos.get(l).getR()).getL()));
              // get position in MS1Trails
              ms1trailPos.add(localpepMax[tempPos.get(l).getL()].get(tempPos.get(l).getR()).getR());

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
            intensityPercentage.add(isoIntPercentage);
            posInMS1Trail.add(ms1trailPos);
          }
        }
      }
    }
    System.out.println("All peptide found before deleting duplicates: " + groupID);
  }

  private void calculateIsoDistributionScore() {
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

  public int searchTarget(double[] array, double[] validify, double target) {
    int pos = InvalidVal;
    double diff = 2 * target * (1 + mzSearchRangePPM); // dumb number
    for (int i = 0 ; i < array.length; i++) {
      if (validify[i] == InvalidVal) {
      } else {
        double cur_diff = Math.abs(array[i] - target);
        if (cur_diff < diff) {
          diff = cur_diff;
          pos = i;
        } else {
          break;
        }
      }
    }
    return pos;
  }
  // Given array of localmax positions, check in range
  // Return position at posArray
  public int searchPos(
          List<Pair<Integer, Integer>> posArray, double[] mzArray, double lowrange, double highrange) {
    int pos = InvalidVal;
    for (int i = 0; i < posArray.size(); i++) {
      if (posArray.get(i).getL() == InvalidVal) continue;
      if (mzArray[posArray.get(i).getL()] >= lowrange && mzArray[posArray.get(i).getL()] <= highrange) {
        pos = i;
        break;
      }
      if (mzArray[posArray.get(i).getL()] > highrange) {
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
    IntensityThreshold = ParameterService.getIntensityThreshold();
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
                  + "intensity_percentage"
                  + '\t'
              + "rt_start"
              + '\t'
              + "rt_end"
              + '\t'
              + "quantification_peaks_sum"
              + '\t'
              + "quantification_peaks_area"
              + '\n';
      printWriter.print(header);

      for (int i = 0; i < scorePepAll.size(); i++) {
        int size = scorePepAll.get(i).size();
        for (int j = 0; j < size; j++) {
          Integer id = i + 1;
          String data = "";
          data += id.toString();
          data += '\t' + scorePepAll.get(i).get(j).getR().toString();
          data += '\t' + scorePepAll.get(i).get(j).getL().toString();
          data += '\t' + scoreZVal.get(i).toString();
          data += '\t' + scoreScanNum.get(i).toString();
          data += '\t' + scoreIntShape.get(i).toString();
          data += '\t' + scoreIsoDistr.get(i).toString();
          data += '\t' + intensityPercentage.get(i).toString();
          data += '\t' + MS1Trails.get(posInMS1Trail.get(i).get(j)).getStartRT();
          data += '\t' + MS1Trails.get(posInMS1Trail.get(i).get(j)).getEndRT();
          data += '\t' + MS1Trails.get(posInMS1Trail.get(i).get(j)).getPeakSum();
          data += '\t' + MS1Trails.get(posInMS1Trail.get(i).get(j)).getPeakArea();
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
