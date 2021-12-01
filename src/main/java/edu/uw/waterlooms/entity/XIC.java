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
  private double startRT;
  private double endRT;

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
    // Set the MZ, Intensity, and RetentionTime values at the maximum intensity of the trail
    setLocalValues(ints, mzs, rts, columeInMap, rowInMap);
    quantify(ints, rts);
  }
  /**
   * Helper function to compute the Maximum Intensity value and set the respective mz and rt fields.
   */
  private void setLocalValues(ArrayList<Double> ints, ArrayList<Double> mzs, ArrayList<Double> rts, ArrayList<Integer> columeInMap, ArrayList<Integer> rowInMap) {
    double maxIntensity = Collections.max(ints);
    int maxIntensityIdx = ints.indexOf(maxIntensity);
    this.maxIntensity = maxIntensity;
    this.mzAtMaxIntensity = mzs.get(maxIntensityIdx);
    this.rtAtMaxIntensity = rts.get(maxIntensityIdx);
    this.columeInMapAtMaxIntensity = columeInMap.get(maxIntensityIdx);
    this.rowInMapAtMaxIntensity = rowInMap.get(maxIntensityIdx);
    this.startRT = rts.get(0);
    this.endRT = rts.get(rts.size() - 1);
  }

  private void quantify(ArrayList<Double> ints, ArrayList<Double> rts){
    this.peakSum = ints.get(0);
    this.peakArea = 0;
    for (int i = 1; i < ints.size(); i++) {
      double intensity = ints.get(i);
      double prev_intensity = ints.get(i-1);
      double rt = rts.get(i);
      double prev_rt = rts.get(i-1);
      this.peakSum += intensity;
      this.peakArea += (intensity + prev_intensity) * (rt - prev_rt) / 2;
    }
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

  public double getStartRT() {
    return startRT;
  }

  public void setStartRT(double startRT) {
    this.startRT = startRT;
  }

  public double getEndRT() {
    return endRT;
  }

  public void setEndRT(double endRT) {
    this.endRT = endRT;
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
