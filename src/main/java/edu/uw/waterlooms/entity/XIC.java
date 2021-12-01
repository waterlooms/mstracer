package edu.uw.waterlooms.entity;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;



public class XIC implements Comparable<XIC>, Serializable {
  final double RT_MARGIN = 0.0000001;

  private double rtAtMaxIntensity;
  private double mzAtMaxIntensity;
  private double maxIntensity;
  private int columeInMapAtMaxIntensity;
  private int rowInMapAtMaxIntensity;
  private double peakSum;
  private double peakArea;
  private ArrayList<Double> intensities;
  private ArrayList<Double> massChargeRatios;
  private ArrayList<Double> retentionTimes;
  private ArrayList<Integer> columeInMap;
  private ArrayList<Integer> rowInMap;

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
  public XIC(ArrayList<Double> ints, ArrayList<Double> mzs, ArrayList<Double> rts, ArrayList<Integer> columeInMap, ArrayList<Integer> rowInMap) {
    this.intensities = ints;
    this.massChargeRatios = mzs;
    this.retentionTimes = rts;
    this.columeInMap = columeInMap;
    this.rowInMap = rowInMap;
    // Set the MZ, Intensity, and RetentionTime values at the maximum intensity of the trail
    setLocalMaxValues();
    quantify();
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
    this.columeInMapAtMaxIntensity = this.columeInMap.get(maxIntensityIdx);
    this.rowInMapAtMaxIntensity = this.rowInMap.get(maxIntensityIdx);
  }

  private void quantify(){
    this.peakSum = this.intensities.get(0);
    this.peakArea = 0;
    for (int i = 1; i < this.intensities.size(); i++) {
      double intensity = this.intensities.get(i);
      double prev_intensity = this.intensities.get(i-1);
      double rt = this.retentionTimes.get(i);
      double prev_rt = this.retentionTimes.get(i-1);
      this.peakSum += intensity;
      this.peakArea += (intensity + prev_intensity) * (rt - prev_rt) / 2;
    }
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

  public double getPeakSum() {
    return peakSum;
  }

  public double getPeakArea() {
    return peakArea;
  }

  public double getStartRT() {
    return retentionTimes.get(0);
  }

  public double getEndRT() {
    return retentionTimes.get(retentionTimes.size()-1);
  }

  @Override
  public int compareTo(XIC xic){
    double x = this.getRtAtMaxIntensity();
    double y = xic.getRtAtMaxIntensity();

    return x > y ? 1: -1;
    /*
    double diff = Math.abs(this.rtAtMaxIntensity - xic.rtAtMaxIntensity);
    if (diff < RT_MARGIN){
      return 1;
    } else if (diff > RT_MARGIN){
      return -1;
    }
    return 0;*/
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
