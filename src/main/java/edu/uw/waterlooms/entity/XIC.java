package edu.uw.waterlooms.entity;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;



public class XIC implements Comparable<XIC>, Serializable {
  private double rtAtMaxIntensity;
  private double mzAtMaxIntensity;
  private double maxIntensity;
  private int columeInMapAtMaxIntensity;
  private int rowInMapAtMaxIntensity;
  private ArrayList<Double> intensities;
  private ArrayList<Double> massChargeRatios;
  private ArrayList<Double> retentionTimes;
  private double precursorMZ;
  private double precursorRT;
  private double peakSum;
  private double peakArea;

  /**
   * Lightweight constructor to allow for mz searching.
   * TODO: REMOVE THIS
   * @param mz double for MZ of a given XIC.
   */
  public XIC(double mz) {
    this.mzAtMaxIntensity = mz;
  }

  /**
   * Lightweight constructor to allow for MZ or RT searching.
   * @param mz double for MZ of a given XIC.
   * @param rt double for RT of a given XIC.
   *
   */
  public XIC(double mz, double rt){
    this.mzAtMaxIntensity = mz;
    this.rtAtMaxIntensity = rt;
  }

  /**
   * Constructor for XIC. Also called a "Trail".
   *
   * @param ints ArrayList of intensities
   * @param mzs ArrayList of mass-charge ratios
   * @param rts ArrayList of retention times
   */
  public XIC(ArrayList<Double> ints, ArrayList<Double> mzs, ArrayList<Double> rts) {
    this.intensities = ints;
    this.massChargeRatios = mzs;
    this.retentionTimes = rts;
    // Set the MZ, Intensity, and RetentionTime values at the maximum intensity of the trail
    setLocalMaxValues();
    // Compute the quantification of peaks
    quantifyPeaks();
  }
  public XIC(ArrayList<Double> ints, ArrayList<Double> mzs, ArrayList<Double> rts, double precursorMZ, double precursorRT) {
    this.intensities = ints;
    this.massChargeRatios = mzs;
    this.retentionTimes = rts;
    this.precursorMZ = precursorMZ;
    this.precursorRT = precursorRT;
    // Set the MZ, Intensity, and RetentionTime values at the maximum intensity of the trail
    setLocalMaxValues();
    // Compute the quantification of peaks
    quantifyPeaks();
  }

  public void quantifyPeaks() {
    ArrayList<Double> ints = this.intensities;
    ArrayList<Double> rts = this.retentionTimes;
    double peaks_sum = ints.get(0);;
    double peaks_area = 0;
    for (int j = 1; j < ints.size(); j++) {
      peaks_sum += ints.get(j);
      peaks_area += (ints.get(j) + ints.get(j-1)) * (rts.get(j) - rts.get(j-1)) / 2;
    }
    this.peakSum = peaks_sum;
    this.peakArea = peaks_area;
  }

  /**
   * Helper function to compute the Maximum Intensity value and set the respective mz and rt fields.
   */
  private void setLocalMaxValues() {
    double maxIntensity = Collections.max(this.intensities);
    int maxIntensityIdx = this.intensities.indexOf(maxIntensity);
    this.maxIntensity = maxIntensity;
    this.mzAtMaxIntensity = this.massChargeRatios.get(maxIntensityIdx);
    this.rtAtMaxIntensity = this.retentionTimes.get(maxIntensityIdx);
  }

  /**
   * Getter for XIC/Trail intensities.
   *
   * @return ArrayList<Double>
   */
  public ArrayList<Double> getIntensities() {
    return this.intensities;
  }

  /**
   * Setter for XIC/Trail intensities.
   *
   * @param intensities ArrayList of intensities to set.
   */
  public void setIntensities(ArrayList<Double> intensities) {
    this.intensities = intensities;
  }

  /**
   * Getter for XIC/Trail mass-charge ratios.
   *
   * @return ArrayList<Double>
   */
  public ArrayList<Double> getMassChargeRatios() {
    return this.massChargeRatios;
  }

  /**
   * Setter for XIC/Trail mass-charge ratios.
   *
   * @param massChargeRatios ArrayList of m/z to set.
   */
  public void setMassChargeRatios(ArrayList<Double> massChargeRatios) {
    this.massChargeRatios = massChargeRatios;
  }

  /**
   * Getter for XIC/Trail retention times.
   *
   * @return ArrayList<Double>
   */
  public ArrayList<Double> getRetentionTimes() {
    return this.retentionTimes;
  }

  /**
   * Setter for XIC/Trail retention times.
   *
   * @param retentionTimes ArrayList of retention times to set.
   */
  public void setRetentionTimes(ArrayList<Double> retentionTimes) {
    this.retentionTimes = retentionTimes;
  }

  /**
   * Getter for RetentionTime at Max Intensity.
   *
   * @return double RetentionTime at the max intensity.
   */
  public double getRtAtMaxIntensity() {
    return rtAtMaxIntensity;
  }

  /**
   * Set the RetentionTime at the Max Intensity.
   *
   * @param rtAtMaxIntensity sets the RetentionTime at MaxIntensity
   */
  public void setRtAtMaxIntensity(double rtAtMaxIntensity) {
    this.rtAtMaxIntensity = rtAtMaxIntensity;
  }

  /**
   * Getter for MZ at Max Intensity.
   *
   * @return double MZ at the Max Intensity.
   */
  public double getMZAtMaxIntensity() {
    return mzAtMaxIntensity;
  }

  /**
   * Set the MZ at the Max Intensity.
   *
   * @param mzAtMaxIntensity sets the MZ at MaxIntensity
   */
  public void setMZAtMaxIntensity(double mzAtMaxIntensity) {
    this.mzAtMaxIntensity = mzAtMaxIntensity;
  }

  /**
   * Getter for Max Intensity value.
   *
   * @return double Max Intensity value.
   */
  public double getMaxIntensity() {
    return maxIntensity;
  }

  /**
   * Set the Maximum Intensity.
   *
   * @param maxIntensity Max Intensity.
   */
  public void setMaxIntensity(double maxIntensity) {
    this.maxIntensity = maxIntensity;
  }

  @Override
  public int compareTo(XIC xic){
    double x = this.getRtAtMaxIntensity();
    double y = xic.getRtAtMaxIntensity();

    return x > y ? 1: -1;
  }

  public int getColumeInMapAtMaxIntensity() {
    return columeInMapAtMaxIntensity;
  }

  public void setColumeInMapAtMaxIntensity(int columeInMapAtMaxIntensity) {
    this.columeInMapAtMaxIntensity = columeInMapAtMaxIntensity;
  }

  public int getRowInMapAtMaxIntensity() {
    return rowInMapAtMaxIntensity;
  }

  public void setRowInMapAtMaxIntensity(int rowInMapAtMaxIntensity) {
    this.rowInMapAtMaxIntensity = rowInMapAtMaxIntensity;
  }

  public double getPrecursorMZ () { return precursorMZ; }

  public void setPrecursorMZ(double precursorMZ) {
    this.precursorMZ = precursorMZ;
  }

  public double getPrecursorRT () { return precursorRT; }

  public void setPrecursorRT(double precursorRT) {
    this.precursorRT = precursorRT;
  }

  public double getPeakArea() { return peakArea; }

  public void setPeakArea(double peakArea) { this.peakArea = peakArea; }

  public double getPeakSum() { return peakSum; }

  public void setPeakSum(double peakSum) { this.peakSum = peakSum; }

  public int getScanNum() {return this.intensities.size(); }

  public int compareRtXIC(XIC xic){
    double x = this.getRtAtMaxIntensity();
    double y = xic.getRtAtMaxIntensity();
    if (x == y){
      return 0;
    }
    return x > y ? 1: -1;
  }

  public int compareMzXIC(XIC xic){
    double x = this.getMZAtMaxIntensity();
    double y = xic.getMZAtMaxIntensity();
    if (x == y){
      return 0;
    }
    return x > y ? 1: -1;
  }

  public int compareIntXIC(XIC xic){
    double x = this.getMaxIntensity();
    double y = xic.getMaxIntensity();
    if (x == y){
      return 0;
    }
    return x > y ? 1: -1;
  }

}
