package edu.uw.waterlooms.msutil;

public class ScanEntry implements Comparable<ScanEntry> {
  public double mz;
  public double rt;
  public double intensity;

  public ScanEntry(double d1, double d2, double d3) {
    mz = d1;
    rt = d2;
    intensity = d3;
  }

  public boolean equals(Object obj) {
    if (obj == null) return false;
    if (obj == this) return true;
    if (!(obj instanceof ScanEntry)) return false;
    ScanEntry o = (ScanEntry) obj;
    return o.mz == this.mz && o.rt == this.rt && o.intensity == this.intensity;
  }

  @Override
  public int compareTo(ScanEntry other) {
    if (mz > other.mz) {
      return 1;
    } else if (mz < other.mz) {
      return -1;
    } else {
      return 0;
    }
  }

  @Override
  public int hashCode() {
    return (Double.valueOf(this.mz).hashCode() + Double.valueOf(this.intensity).hashCode() + Double.valueOf(this.rt).hashCode());
  }
}
