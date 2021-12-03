package edu.uw.waterlooms.entity;

public class SVRScore {
  private double id;
  private double massChargeRatio;
  private double retentionTime;
  private double charge;
  private double isotopeNumber;
  private double intensityShapeScore;
  private double isotopeDistributionScore;
  private double intensityAreaPercentage;
  private double retentionTimeStart;
  private double retentionTimeEnd;
  private double scanNum;
  private double peaksSum;
  private double peaksArea;
  private double svrScore;
  private double qualityScore;
  /**
   * Constructor for the model_result object returned from FeatureSelect.java. Represents a single precursor detected.
   *
   * @param id ID of Match
   * @param massChargeRatio m/z of precursor
   * @param retentionTime rt of precursor
   * @param charge z of precursor
   * @param isotopeNumber
   * @param intensityShapeScore
   * @param isotopeDistributionScore
   * @param intensityAreaPercentage
   * @param retentionTimeStart
   * @param retentionTimeEnd
   * @param scanNum
   * @param svrScore
   * @param qualityScore
   */
  public SVRScore(
      double id,
      double massChargeRatio,
      double retentionTime,
      double charge,
      double isotopeNumber,
      double intensityShapeScore,
      double isotopeDistributionScore,
      double intensityAreaPercentage,
      double retentionTimeStart,
      double retentionTimeEnd,
      double scanNum,
      double peaksArea,
      double peaksSum,
      double svrScore,
      double qualityScore) {
    this.id = id;
    this.massChargeRatio = massChargeRatio;
    this.retentionTime = retentionTime;
    this.charge = charge;
    this.isotopeNumber = isotopeNumber;
    this.intensityShapeScore = intensityShapeScore;
    this.isotopeDistributionScore = isotopeDistributionScore;
    this.intensityAreaPercentage = intensityAreaPercentage;
    this.retentionTimeStart = retentionTimeStart;
    this.retentionTimeEnd = retentionTimeEnd;
    this.scanNum = scanNum;
    this.peaksSum = peaksSum;
    this.peaksArea = peaksArea;
    this.svrScore = svrScore;
    this.qualityScore = qualityScore;
  }

  public double getId() {
    return id;
  }

  public void setId(double id) {
    this.id = id;
  }

  public double getMassChargeRatio() {
    return massChargeRatio;
  }

  public void setMassChargeRatio(double massChargeRatio) {
    this.massChargeRatio = massChargeRatio;
  }

  public double getRetentionTime() {
    return retentionTime;
  }

  public void setRetentionTime(double retentionTime) {
    this.retentionTime = retentionTime;
  }

  public double getCharge() {
    return charge;
  }

  public void setCharge(double charge) {
    this.charge = charge;
  }

  public double getIsotopeNumber() {
    return isotopeNumber;
  }

  public void setIsotopeNumber(double isotopeNumber) {
    this.isotopeNumber = isotopeNumber;
  }

  public double getIntensityShapeScore() {
    return intensityShapeScore;
  }

  public void setIntensityShapeScore(double intensityShapeScore) {
    this.intensityShapeScore = intensityShapeScore;
  }

  public double getIsotopeDistributionScore() {
    return isotopeDistributionScore;
  }

  public void setIsotopeDistributionScore(double isotopeDistributionScore) {
    this.isotopeDistributionScore = isotopeDistributionScore;
  }

  public double getIntensityAreaPercentage() {
    return intensityAreaPercentage;
  }

  public void setIntensityAreaPercentage(double intensityAreaPercentage) {
    this.intensityAreaPercentage = intensityAreaPercentage;
  }

  public double getRetentionTimeStart() {
    return retentionTimeStart;
  }

  public void setRetentionTimeStart(double retentionTimeStart) {
    this.retentionTimeStart = retentionTimeStart;
  }

  public double getRetentionTimeEnd() {
    return retentionTimeEnd;
  }

  public void setRetentionTimeEnd(double retentionTimeEnd) {
    this.retentionTimeEnd = retentionTimeEnd;
  }

  public double getSvrScore() {
    return svrScore;
  }

  public void setSvrScore(double svrScore) {
    this.svrScore = svrScore;
  }

  public double getQualityScore() {
    return qualityScore;
  }

  public void setQualityScore(double qualityScore) {
    this.qualityScore = qualityScore;
  }
}
