package edu.uw.waterlooms.msutil;

public class SpecEntry implements Comparable<SpecEntry> {
  public double rt;
  public double mz;
  public double spec_int;
  public double cor;

  public SpecEntry(double d1, double d2, double d3) {
    mz = d1;
    rt = d2;
    spec_int = d3;
  }

  public boolean equals(Object obj) {
    if (obj == null) return false;
    if (obj == this) return true;
    if (!(obj instanceof SpecEntry)) return false;
    SpecEntry o = (SpecEntry) obj;
    return o.mz == this.mz && o.rt == this.rt;
  }

  @Override
  public int compareTo(SpecEntry other) {
    return Double.compare(cor, other.cor);
  }

  @Override
  public int hashCode() {
    return (Double.valueOf(this.mz).hashCode() + Double.valueOf(this.rt).hashCode());
  }
}
