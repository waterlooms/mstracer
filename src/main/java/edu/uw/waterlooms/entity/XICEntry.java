package edu.uw.waterlooms.entity;

public class XICEntry {
  private double mz;
  private double rt;
  private double intensity;

  public XICEntry(double mz, double rt, double intensity) {
    this.mz = mz;
    this.rt = rt;
    this.intensity = intensity;
  }
  public double getMz(){
    return this.mz;
  }
  public double getRt(){
    return this.rt;
  }
  public double getIntensity(){
    return this.intensity;
  }

}
