/** Copyright Rapid Novor Inc. */
package edu.uw.waterlooms.ms1;

import java.util.Arrays;

/**
 * Compute the theoretical isotope shape using the average abundance of C, H, O, N, S elements.
 *
 * @author binma
 */
public class TheoreticalIsotope {
  private static TheoreticalIsotope defaultInstance;

  /**
   * The elements' properties
   *
   * @author binma
   */
  enum Element {
    // source https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
    C(12, 0.317f, new double[] {0.9893, 0.0107}, new double[] {12.0, 13.003355}),
    H(1, 0.494f, new double[] {0.999885, 0.000115}, new double[] {1.007825, 2.014102}),
    O(
        16,
        0.097f,
        new double[] {0.99757, 0.00038, 0.00205},
        new double[] {15.994915, 16.99131, 17.999160}),
    N(14, 0.088f, new double[] {0.99636, 0.00364}, new double[] {14.003074, 15.001090}),
    S(
        32,
        0.003f,
        new double[] {0.9499, 0.0075, 0.0425, 1e-12, 0.0001},
        new double[] {31.972071, 32.971459, 33.967867, 34.96903216, 35.96708076});
    private int nominalMass;
    private double frequency; // the frequency of the element in an average peptide
    private double[] abundance; // the natural abundance of the isotopes of the element.
    private double[] mass;

    Element(int nominalmass, float frequency, double[] abundance, double[] mass) {
      this.nominalMass = nominalmass;
      this.frequency = frequency;
      this.abundance = abundance;
      this.mass = mass;
      float total = 0;
      for (int i = 0; i < abundance.length; i++) {
        total += this.abundance[i];
      }
      assert (Math.abs(total - 1) < 1e-7);
      assert mass.length == abundance.length;
      for (int i = 0; i < abundance.length; i++) {
        this.abundance[i] /= total;
      }
    }

    int getNominalMass() {
      return nominalMass;
    }

    double getFrequency() {
      return frequency;
    }

    double[] getAbundance() {
      return abundance;
    }

    double[] getDeltaMass() {
      double[] deltaMass = new double[mass.length];
      for (int i = 0; i < mass.length; i++) {
        deltaMass[i] = mass[i] - mass[0];
      }
      return deltaMass;
    }
  }

  private int npeaks; // the maximum number of isotope peaks to compute.
  private int maxM;
  private int massStep;
  private double
      minRelativeHeight; // only keep the isotope peaks if they are taller than 0.01*the maximum
  // isotope peak height
  private float[][] shapes; // shapes[i] is the shapes for mass i*massStep daltons

  /**
   * Create a IsotopeShape calculator.
   *
   * @param nIsotopes maximum number of isotopes to compute (including the monoisotope)
   * @param maxMass maximum mass value to use.
   * @param massStep store the values for every massStep Dalton. isotope shapes do not change
   *     significantly for changing of massStep Daltons
   * @param minRelativeHeight the minimum relative height of maximum isotope peak to be reported.
   *     The isotope peaks at left of the maximum isotope peak are always reported.
   */
  public TheoreticalIsotope(int nIsotopes, int maxMass, int massStep, double minRelativeHeight) {
    this.npeaks = nIsotopes;
    this.maxM = maxMass;
    this.massStep = massStep;
    this.minRelativeHeight = minRelativeHeight;
    double[][] isoshapes = new double[maxMass + 1][nIsotopes];
    // init
    isoshapes[0][0] = 1;
    // dynamic programming
    for (int m = 1; m <= maxMass; m++) {
      double[] current = new double[nIsotopes];
      for (Element e : Element.values()) {
        double abundance[] = e.getAbundance();
        int mass = e.getNominalMass();
        if (mass <= m) {
          for (int i = 0; i < nIsotopes; i++) {
            for (int j = 0; j < Math.min(nIsotopes - i, abundance.length); j++) {
              current[i + j] += isoshapes[m - mass][i] * abundance[j] * e.getFrequency();
            }
          }
        }
      }
      // normalize and copy current
      float total = 0;
      for (int i = 0; i < nIsotopes; i++) {
        total += current[i];
      }
      if (total > 0) {
        for (int i = 0; i < nIsotopes; i++) {
          isoshapes[m][i] = current[i] / total;
        }
      }
    }
    // copy to shapes for later use, reduce memory by copying only the multiple of massStep.
    shapes = new float[maxMass / massStep + 1][]; // keep the shapes for every massStep daltons
    for (int m = 0; m < shapes.length; m++) {
      shapes[m] = trim(isoshapes[m * massStep]);
    }
  }

  // trim an isotope shape to only keep the first several isotope peaks with significant abundance
  private float[] trim(double shape[]) {
    double max = shape[0];
    for (int i = 1; i < shape.length; i++) {
      if (shape[i] > max) max = shape[i];
    }
    double min = max * minRelativeHeight;
    int n = shape.length;
    for (; n > 1; n--) {
      if (shape[n - 1] >= min) break;
    }
    float[] trimmed = new float[n];
    for (int i = 0; i < n; i++) {
      trimmed[i] = (float) shape[i];
    }
    return trimmed;
  }

  /**
   * create a default isotope shape calculator
   *
   * @return
   */
  public static TheoreticalIsotope getDefaultInstance() {
    if (defaultInstance == null) {
      defaultInstance = new TheoreticalIsotope(50, 20000, 30, 0.01);
    }
    return defaultInstance;
  }

  public float[] computeShape(float m) {
    int nm = (int) (m * 0.999506242f / massStep); // mass to nominal mass
    if (nm >= shapes.length) {
      throw new RuntimeException(
          "The computeShape(" + m + ") exceeds the maximum mass TheoreticalIsotope supports");
    } else {
      return Arrays.copyOf(shapes[nm], shapes[nm].length);
    }
  }
}
