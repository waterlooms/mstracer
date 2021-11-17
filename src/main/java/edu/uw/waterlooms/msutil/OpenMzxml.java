package edu.uw.waterlooms.msutil;

import umich.ms.datatypes.LCMSData;
import umich.ms.datatypes.LCMSDataSubset;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scan.StorageStrategy;
import umich.ms.datatypes.scancollection.IScanCollection;
import umich.ms.datatypes.scancollection.ScanIndex;
import umich.ms.datatypes.scancollection.impl.ScanCollectionDefault;
import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.mzml.MZMLFile;
import umich.ms.fileio.filetypes.mzxml.MZXMLFile;
import umich.ms.fileio.filetypes.mzxml.MZXMLIndex;

import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

public class OpenMzxml {
  String filePath;
  public Set<Map.Entry<Integer, IScan>> scanEntries = null;
  public Set<Map.Entry<Integer, IScan>> scanEntries2 = null;
  public TreeMap<Integer, IScan> num2scan = null;
  public TreeMap<Integer, IScan> num2scan2 = null;
  public IScanCollection scans = null;

  public OpenMzxml(String str) {
    filePath = str;
    try {
      init();
    } catch (FileParsingException e) {
      e.printStackTrace();
    }
  }

  // different way to load data
  void init2() throws FileParsingException {
    MZXMLFile source = new MZXMLFile(filePath);
    // if a scan has zero peaks in its spectrum it will still be parsed
    source.setExcludeEmptyScans(false);
    // null means use as many cores as reported by Runtime.getRuntime().availableProcessors()
    source.setNumThreadsForParsing(null);
    // 30 sec timeout for worker threads - each worker must parse its chunk of spectra within that
    // time
    source.setParsingTimeout(30L);

    MZXMLIndex mzxmlIndex = source.fetchIndex();

    // this is a data structure used to store scans and to navigate around the run
    scans = new ScanCollectionDefault();
    // softly reference spectral data, make it reclaimable by GC
    scans.setDefaultStorageStrategy(StorageStrategy.SOFT);
    // set it to automatically re-parse spectra from the file if spectra were not yet parsed or were
    // reclaimed
    // to make auto-loading work you'll need to use IScan#fetchSpectrum() method instead of
    // IScan#getSpectrum()
    scans.isAutoloadSpectra(true);

    // set the MZXML file as the data source for this scan collection
    scans.setDataSource(source);
    // load the whole run, with forced parsing of MS2 spectra, using default StorageStrategy.
    scans.loadData(LCMSDataSubset.MS2_WITH_SPECTRA);

    //        scans.getMapMsLevel2index()
    //        num2scan2 = scans.getMapNum2scan();
    //        scanEntries2 = num2scan2.entrySet();

  }

  // load data to cache
  void init() throws FileParsingException {
    LCMSData data = null;
    if (filePath.toLowerCase().endsWith(".mzml")) {
      MZMLFile source = new MZMLFile(filePath);
      data = new LCMSData(source);
    } else if (filePath.toLowerCase().endsWith(".mzxml")) {
      MZXMLFile source = new MZXMLFile(filePath);
      data = new LCMSData(source);
    }
    // Get the index (fetchXXX() methods will parse data from the file if it has not yet been
    // parsed) and
    // cache it in the object for reuse.
    // You'll only need the index if you want to convert between internal scan numbers and raw scan
    // numbers
    // in the file. Some files might have non-consecutive scan numbers, for example, but internally
    // they'll be
    // renumbered to start from 1 and increment by one for each next scan.
    //        MZXMLIndex idx = source.fetchIndex();
    // info about the run
    //        LCMSRunInfo runInfo = source.fetchRunInfo();

    // To parse a single scan from the file (or a range of scans) we first create a predicate
    // matching the
    // scan to be parsed.
    // For example, parse scans from 1 to 3 at MS level 2.
    //        Set<Integer> msLevel = Collections.singleton(1);
    //        LCMSDataSubset subset = new LCMSDataSubset(1, 3, msLevel, null);
    //        List<IScan> parsedScans = source.parse(subset);

    // If you want higher level access to data, create an LCMSData object
    // necessary
    //        LCMSData data = new LCMSData(source);
    // load the whole structure of the run, and parse all spectra for MS1 scans
    data.load(LCMSDataSubset.WHOLE_RUN);
    //        data.releaseMemory();

    // or load the whole structure, but only get m/z-intensity info at MS level 2
    //        data.load(new LCMSDataSubset(null, null, msLevel, null));
    //        data.releaseMemory();
    // alternatively, use this shortcut
    //        data.load(LCMSDataSubset.MS1_WITH_SPECTRA);
    //        data.releaseMemory();

    // If you need memory management, you can also pass an instance of an object, which will be
    // considered
    // the owner of prased data. When this object is garbage collected, this will be detected
    // automatically
    // and corresponding spectra released.
    //        Object dataUser = new Object();
    //        data.load(LCMSDataSubset.WHOLE_RUN, dataUser);
    //        System.out.printf("The data is loaded and used by [%s] object.\n",
    // System.identityHashCode(dataUser));
    // at this point dataUser might be garbage collected as it's not referenced anymore, and the
    // data might
    // get unloaded automatically
    //        dataUser = null; // just to be sure that we don't have a strong reference

    // If you don't want to fiddle around with memory management at all, but still want it to play
    // nicely
    // there's one more feature - auto-loading of spectra.
    // You can parse the whole structure of the file and keep it in memory (it's rather small), and
    // just magically get the spectra whenever you need them.
    // Also set referenceing type to soft, so that garbage collector could reclaim unused spectra.
    //        data.load(LCMSDataSubset.STRUCTURE_ONLY);

    scans = data.getScans();
    scans.isAutoloadSpectra(true); // set automatic spectra loading
    //        scans.setDefaultStorageStrategy(StorageStrategy.SOFT); // mz-intensity data will be
    // softly referenced
    TreeMap<Integer, ScanIndex> msLevel2index = scans.getMapMsLevel2index();
    ScanIndex ms2idx = msLevel2index.get(1); // get the index at MS level 1

    // we'll iterate by scan numbers
    // num2scan is a tree map; key: scan number, value: scan object (ms1)
    num2scan = ms2idx.getNum2scan();
    // scanEntries is a set view of the mappings contained in this map (ms1)
    scanEntries = num2scan.entrySet();

    ScanIndex ms2idx2 = msLevel2index.get(2); // get the index at MS level 2

    // we'll iterate by scan numbers
    // same num2scan and scanEntries for ms2
    num2scan2 = ms2idx2.getNum2scan();
    scanEntries2 = num2scan2.entrySet();
  }
}
