package uk.ac.ox.well.hts.reliablegenome;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.Stack;
import java.util.TreeSet;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

import at.ac.univie.cs.mis.lds.index.itree.IncompatibleIntervalsException;
import at.cibiv.ngs.tools.bed.AbstractBedIterator;
import at.cibiv.ngs.tools.bed.BedIterator;
import at.cibiv.ngs.tools.bed.MemoryBedIterator;
import at.cibiv.ngs.tools.bed.SimpleBEDFile;
import at.cibiv.ngs.tools.lds.GenomicITree;
import at.cibiv.ngs.tools.lds.GenomicInterval;
import at.cibiv.ngs.tools.util.CanonicalChromsomeComparator;
import at.cibiv.ngs.tools.util.GenomicPosition;
import at.cibiv.ngs.tools.util.Histogram;
import at.cibiv.ngs.tools.util.MathUtil;
import at.cibiv.ngs.tools.util.Statistics;
import at.cibiv.ngs.tools.util.StringUtils;
import at.cibiv.ngs.tools.util.TabIterator;
import at.cibiv.ngs.tools.vcf.HeterozygousityClassification.GENOTYPE;
import at.cibiv.ngs.tools.vcf.SimpleVCFFile;
import at.cibiv.ngs.tools.vcf.SimpleVCFVariant;
import at.cibiv.ngs.tools.vcf.VCFContextWinIterator;
import at.cibiv.ngs.tools.vcf.VCFIterator;
import at.cibiv.ngs.tools.vcf.VCFSampleIterator;
import at.cibiv.ngs.tools.vcf.VCFWriter;
import at.cibiv.ngs.tools.wig.WigOutputStream;
import at.cibiv.ngs.tools.wig.WigTools;

/**
 * Main RG class.
 * 
 * @author niko@well.ox.ac.uk
 * 
 */
public class CalcreliabilitySignals {

    public static final String CMD = "CalcreliabilitySignals";
    public static final String CMD_INFO = "Calculate signals from vcf files.";
    public static boolean debug = false;

    public final static String FILTER_STRING_UNRELIABLE = "Unreliable";
    public final static String FILTER_STRING_NO_CLASS = "Unclassified";

    /**
     * Go to next entry in a VCF iterator that is not filtered and in the set of
     * tested chromosomes.
     * 
     * @param it
     * @param testedChroms
     * @param excludeNonPass
     * @return
     */
    private static SimpleVCFVariant incVcfIt(Iterator<SimpleVCFVariant> it, Set<String> testedChroms, boolean excludeNonPass) {
	SimpleVCFVariant nxt = null;
	do {
	    if (!it.hasNext())
		return null;
	    nxt = it.next();
	} while ((excludeNonPass && !(nxt.getFilter().equalsIgnoreCase("PASS") || nxt.getFilter().equals(".")))
		|| (testedChroms != null && !testedChroms.contains(nxt.getChromosomeOrig())));
	return nxt;
    }

    /**
     * Go to next entry in a BED iterator that is in the set of tested
     * chromosomes.
     * 
     * @param bi
     * @param testedChroms
     * @return
     */
    private static GenomicInterval incBedIt(AbstractBedIterator bi, Set<String> testedChroms) {
	GenomicInterval nxt = null;
	do {
	    if (!bi.hasNext())
		return null;
	    nxt = bi.next();
	} while ((testedChroms != null && !testedChroms.contains(nxt.getOriginalChrom())));
	return nxt;
    }

    /**
     * Check whether the passed genomic position is in the set of passed
     * chromosomes and inculded regions and not in the passed set of excluded
     * regions.
     * 
     * @param pos
     * @param testedChroms
     * @param included
     * @param excluded
     * @return true if the position is included or false otherwise.
     */
    private static boolean checkPositionIncluded(GenomicPosition pos, Set<String> testedChroms, List<AbstractBedIterator> included,
	    List<AbstractBedIterator> excluded) {
	// check included/excluded lists
	boolean isIncluded = true;

	for (AbstractBedIterator bi : excluded) {
	    if (bi == null) {
		// this means: exclude everything!
		isIncluded = false;
		break;
	    }
	    // get current interval
	    while (bi.hasNext() && bi.getCurrentInterval() != null && bi.getCurrentInterval().getRightPosition().compareTo(pos) < 0) {
		incBedIt(bi, testedChroms); // load next
		// interval

	    }
	    if (bi.getCurrentInterval() != null && bi.getCurrentInterval().contains(pos))
		isIncluded = false;
	}

	for (AbstractBedIterator bi : included) {
	    // get current interval
	    while (bi.hasNext() && bi.getCurrentInterval() != null && bi.getCurrentInterval().getRightPosition().compareTo(pos) < 0) {
		incBedIt(bi, testedChroms); // load next
		// interval
	    }
	    if (bi.getCurrentInterval() != null && bi.getCurrentInterval().contains(pos))
		isIncluded = true;
	}

	// System.out.println(pos + " " + isIncluded);
	return isIncluded;
    }

    /**
     * Main method for calculating a genomic partition.
     * 
     * @param vcfFiles
     * @param excludeFiles
     * @param includeFiles
     * @param testedChromsArray
     * @param scoreConcordantHit
     * @param scoreDiscordantHit
     * @param dontCheckSort
     * @param excludeNonPass
     * @param dropAllFiltered
     * @param createWigs
     * @param winSize
     * @param reliableMin
     * @param unreliableMax
     * @param smallWinSizeFactor
     * @param smallWindowSizeRegions
     * @param outdir
     * @return
     * @throws IOException
     */
    static Statistics calc(String[] vcfFiles, String[] excludeFiles, String[] includeFiles, String[] testedChromsArray, int scoreConcordantHit,
	    int scoreDiscordantHit, boolean dontCheckSort, boolean excludeNonPass, boolean dropAllFiltered, boolean createWigs, int winSize,
	    double reliableMin, double unreliableMax, double smallWinSizeFactor, String smallWindowSizeRegionsF, File outdir) throws IOException {

	long start = System.currentTimeMillis();
	Statistics stats = new Statistics();
	stats.setString("Num of Input VCF Files", vcfFiles.length + "");
	stats.setString("Input VCF Files", Arrays.toString(vcfFiles));
	stats.setString("Excluded Regions", Arrays.toString(excludeFiles));
	stats.setString("included Regions", Arrays.toString(includeFiles));
	stats.setString("Tested Chromosomes", Arrays.toString(testedChromsArray));
	stats.setString("Score for concordant position", scoreConcordantHit + "");
	stats.setString("Score for discordant position", scoreDiscordantHit + "");
	stats.setString("ExcludeNonPass", excludeNonPass + "");
	stats.setString("DropAllFiltered", dropAllFiltered + "");
	stats.setString("DontCheckSort", dontCheckSort + "");
	stats.setString("Window size", winSize + "");
	stats.setString("Reliable Min", reliableMin + "");
	stats.setString("Unreliable Max", unreliableMax + "");
	stats.setString("Small Window Size Factor", smallWinSizeFactor + "");
	stats.setString("Regions with small window size", smallWindowSizeRegionsF + "");
	stats.setString("Outdir", outdir.toString());
	if (debug)
	    System.out.println(stats);

	// prepare chromosome scope.
	Set<String> testedChroms = null;
	if (testedChromsArray != null) {
	    testedChroms = new HashSet<String>();
	    for (String tc : testedChromsArray)
		testedChroms.add(tc);
	}

	// prepare smallWindowRegions
	GenomicITree smallWindowSizeRegions = null;
	if (smallWindowSizeRegionsF != null) {
	    if (debug)
		System.out.println("loading regions with smaller (" + smallWinSizeFactor + ") window size from " + smallWindowSizeRegionsF);
	    smallWindowSizeRegions = (new SimpleBEDFile(new File(smallWindowSizeRegionsF))).getGenomicITree();
	}

	// sorted list of current pointers
	SortedSet<GenomicPosition> pointers = new TreeSet<GenomicPosition>();

	// prepare include/exclude BEDs
	if (debug)
	    System.out.println("Prepare included/excluded regions");
	List<AbstractBedIterator> included = new ArrayList<AbstractBedIterator>();
	if (includeFiles != null)
	    for (String inc : includeFiles) {
		if (inc.contains(":")) {
		    MemoryBedIterator it = new MemoryBedIterator();
		    it.setName("includes");
		    it.add(inc);
		    GenomicInterval firstInterval = incBedIt(it, testedChroms);
		    if (firstInterval != null) { // there is at least one
			// interval
			included.add(it);
		    }
		} else {
		    BedIterator it = new BedIterator(new File(inc));
		    it.setName("includes");
		    GenomicInterval firstInterval = incBedIt(it, testedChroms);
		    if (firstInterval != null) { // there is at least one
			// interval
			included.add(it);
		    }
		}
	    }
	List<AbstractBedIterator> excluded = new ArrayList<AbstractBedIterator>();
	if (excludeFiles != null)
	    for (String ex : excludeFiles) {
		if (ex.equals("ALL")) {
		    excluded.clear();
		    excluded.add(null);
		    break;
		} else if (ex.contains(":")) {
		    MemoryBedIterator it = new MemoryBedIterator();
		    it.setName("excludes");
		    it.add(ex);
		    GenomicInterval firstInterval = incBedIt(it, testedChroms);
		    if (firstInterval != null) { // there is at least one
			// interval
			excluded.add(it);
		    }
		} else {
		    BedIterator it = new BedIterator(new File(ex));
		    it.setName("excludes");
		    GenomicInterval firstInterval = incBedIt(it, testedChroms);
		    if (firstInterval != null) { // there is at least one
			// interval
			excluded.add(it);
		    }
		}
	    }

	// prepare iterators
	if (debug)
	    System.out.println("Prepare input files");
	List<VCFIterator> vits = new ArrayList<VCFIterator>();
	for (String vcffn : vcfFiles) {

	    String sample = null;
	    if (vcffn.contains("@")) {
		String[] tmp = vcffn.split("@");
		sample = tmp[0];
		vcffn = tmp[1];
	    }
	    File vcf = new File(vcffn);
	    // check VCF sorting?
	    if (!dontCheckSort) {
		if (debug)
		    System.out.println("Check sorting of " + vcf);
		if (!SimpleVCFFile.checkSortedVCF(vcf, testedChroms))
		    throw new IOException("File " + vcf + " is not sorted by coords");
	    }
	    VCFIterator it = null;
	    if (sample == null)
		it = new VCFIterator(vcf);
	    else
		it = new VCFSampleIterator(vcf, sample);

	    SimpleVCFVariant v = incVcfIt(it, testedChroms, excludeNonPass);
	    if (v != null)
		pointers.add(v.getGenomicPosition());
	    vits.add(it);
	}

	if (pointers.size() == 0)
	    throw new IOException("All iterators were empty - do chrom names match the tested ones (" + testedChroms
		    + ") ? Does the file contain GT information!?");

	// create output dir
	if (!outdir.exists())
	    outdir.mkdirs();

	// filename prefix
	String prefix = "RG-win" + winSize + "-score" + scoreConcordantHit + "_" + scoreDiscordantHit + "-";
	if (smallWindowSizeRegions != null)
	    prefix += "sw" + new File(smallWindowSizeRegionsF).getName() + "-";

	// prepare output streams
	WigOutputStream powerOut = createWigs ? new WigOutputStream(new File(outdir, prefix + "power.wig")) : null;
	WigOutputStream concordanceOut = createWigs ? new WigOutputStream(new File(outdir, prefix + "concordance.wig")) : null;
	// the concordance.win file is always written as it is needed to calc
	// the reliable/unreliable regions. It will be deleted if createWigs=F
	File concordanceOutWinFile = new File(outdir, prefix + "concordance.win" + winSize + ".wig");
	WigOutputStream concordanceOutWin = new WigOutputStream(concordanceOutWinFile);
	BACKUPSignalInterpolator concordanceInterpol = new BACKUPSignalInterpolator(winSize, smallWinSizeFactor, smallWindowSizeRegions);

	File dataFile = new File(outdir, prefix + "data.tsv");
	PrintStream dataOut = new PrintStream(dataFile);
	dataOut.println("## pos: genomic position");
	dataOut.println("## distance: genomic distance to last matched position; is zero for the first entry per chrom.");
	dataOut.println("## power: Number of datasets supporting this position.");
	dataOut.println("## score: reliability score.");
	dataOut.println("## score for concordant position" + scoreConcordantHit);
	dataOut.println("## score for discordant position" + scoreDiscordantHit);
	dataOut.println("pos\tdistance\tpower\tscore");
	stats.setString("DataFile", dataFile.toString());

	File unionFile = new File(outdir, prefix + "union.vcf");
	VCFWriter unionout = new VCFWriter(new PrintStream(unionFile));
	List<SimpleVCFVariant> toEmit = new ArrayList<SimpleVCFVariant>();
	stats.setString("UnionFile", unionFile.toString());

	long lastpos = 0L;
	String lastChr = "";
	int count = 0;
	String lastCtx = "?"; // indicates whether the last considered position
	// was reliable or unreliable.
	if (debug)
	    System.out.println("Start");
	// loop over variants
	do {

	    // filled with all matches at pos <pointers.first()>
	    List<SimpleVCFVariant> matchVariants = new ArrayList<SimpleVCFVariant>();
	    int filteredVariantsCount = 0; // keeps track of how many filtered
	    // variants are in the current hit.

	    // check iterators
	    for (int i = 0; i < vits.size(); i++) {
		VCFIterator it = vits.get(i);

		if (it.getCurrentVariant().getGenomicPosition().equals(pointers.first())) {

		    // match at current position
		    matchVariants.add(it.getCurrentVariant());
		    if (it.getCurrentVariant().isFiltered())
			filteredVariantsCount++;

		    if (it.hasNext()) {
			// handle multiple variants at same position!
			SimpleVCFVariant var = incVcfIt(it, testedChroms, excludeNonPass);
			if (var != null) {
			    GenomicPosition nxt = var.getGenomicPosition();
			    while (nxt.equals(pointers.first())) {
				var = incVcfIt(it, testedChroms, excludeNonPass);
				if (var == null) {
				    break;
				}
				nxt = var.getGenomicPosition();
			    }
			    // add next position
			    pointers.add(nxt);
			}
		    }

		}
	    }

	    stats.inc("Raw-matched-positions");

	    boolean positionIncluded = checkPositionIncluded(pointers.first(), testedChroms, included, excluded);

	    // drop matches where all variants are filtered (applies only if
	    // filtered vars are included!)
	    if (dropAllFiltered && matchVariants.size() == filteredVariantsCount)
		positionIncluded = false;

	    // calc distance to last position
	    long distance = 0;
	    if (pointers.first().getChromosomeOriginal().equals(lastChr))
		distance = pointers.first().get1Position() - lastpos;
	    lastpos = pointers.first().get1Position();
	    lastChr = pointers.first().getChromosomeOriginal();

	    if (positionIncluded) {
		stats.inc("Matched-positions");

		int power = matchVariants.size();
		ScoreResults scoreResults = score(matchVariants, scoreConcordantHit, scoreDiscordantHit);

		if (powerOut != null)
		    powerOut.pushAndFillZero(pointers.first(), power);
		if (concordanceOut != null)
		    concordanceOut.pushAndFillZero(pointers.first(), scoreResults.score);

		concordanceInterpol.addDataPoint(pointers.first(), scoreResults.score, concordanceOutWin);
		// powerInterpol.addDataPoint(pointers.first(), power,
		// powerOutWin);

		dataOut.println(pointers.first().toString1basedCanonicalHuman() + "\t" + distance + "\t" + power + "\t" + scoreResults.score);

		// create UNION file. Create one variant per refAllele
		Map<String, SimpleVCFVariant> outVariants = new HashMap<String, SimpleVCFVariant>();
		Map<String, Set<String>> outVariantAltAlleles = new HashMap<String, Set<String>>();
		for (SimpleVCFVariant x : matchVariants) {
		    SimpleVCFVariant v = outVariants.get(x.getRefString());
		    if (v == null) {
			v = new SimpleVCFVariant();
			v.setChromosomeOrig(matchVariants.get(0).getChromosomeOrig());
			v.setPosition(matchVariants.get(0).getPosition());
			v.setID(matchVariants.get(0).getID());
			v.setRefString(x.getRefString());
			v.setFilter(scoreToFilter(scoreResults.score, reliableMin, unreliableMax));
			v.addField("power", power + "");
			v.addField("score", scoreResults.score + "");
			// call percentage ordered by caller id.
			v.addField("CP", StringUtils.concat(scoreResults.CP, ","));
			v.addField("AVGQUAL", StringUtils.concat(scoreResults.AVGQUAL, ","));
			if (scoreResults.AVGMDI != null)
			    v.addField("AVGMDI", scoreResults.AVGMDI + "");
			if (scoreResults.AVGDPRAW != null)
			    v.addField("AVGDPRAW", scoreResults.AVGDPRAW + "");
			v.addField("GTHOM", MathUtil.fmt(scoreResults.getHomFrac()));
			v.addField("GTHET", MathUtil.fmt(scoreResults.getHetFrac()));
			v.addField("GTOTHER", MathUtil.fmt(scoreResults.getOtherFrac()));
			// AF = 2 * HOM + HET / 2 * ALL
			v.addField("AF", MathUtil.fmt( (2*scoreResults.getHomFrac() + scoreResults.getHetFrac() ) / ( 2d * (double) vcfFiles.length ) ) );
		    }

		    Set<String> altAlleles = outVariantAltAlleles.get(x.getRefString());
		    if (altAlleles == null)
			altAlleles = new HashSet<String>();
		    for (String a : x.getAltAlleles())
			altAlleles.add(a);
		    outVariants.put(x.getRefString(), v);
		    outVariantAltAlleles.put(x.getRefString(), altAlleles);
		}

		// set context for variants from previous position
		for (int i = toEmit.size() - 1; i >= 0; i--) {

		    // previous emitted variant
		    String ctx = "";
		    if (i > 0)
			ctx += filterToContextChar(toEmit.get(i - 1).getFilter());
		    else
			ctx += lastCtx;

		    // current variant
		    ctx += filterToContextChar(toEmit.get(i).getFilter());

		    // next emitted variant
		    if (i < toEmit.size() - 1)
			ctx += filterToContextChar(toEmit.get(i + 1).getFilter());
		    else {
			if (outVariants.size() == 0)
			    ctx += "0";
			else
			    ctx += filterToContextChar(outVariants.values().iterator().next().getFilter());
		    }

		    toEmit.get(i).addField("Context", ctx);
		}
		if (toEmit.size() > 0)
		    lastCtx = filterToContextChar(toEmit.get(toEmit.size() - 1).getFilter());
		else
		    lastCtx = "?";

		// write variants from previous position
		unionout.add(toEmit);
		toEmit = new ArrayList<SimpleVCFVariant>();

		// emit new variants
		for (SimpleVCFVariant x : outVariants.values()) {
		    x.setAltString(StringUtils.concat(outVariantAltAlleles.get(x.getRefString()), ","));
		    toEmit.add(x);
		}

	    } // included

	    if (debug)
		if (++count % 10000 == 0) {
		    System.out.println("> " + pointers.first().toString1based());
		}

	    // go to next position
	    pointers.remove(pointers.first());

	} while (!pointers.isEmpty());

	// write last variants!
	unionout.add(toEmit);

	// write stats
	stats.setString("Execution time [ms]", (System.currentTimeMillis() - start) + "");
	stats.toFile(new PrintStream(new File(outdir, "README")));

	if (powerOut != null)
	    powerOut.close();
	if (concordanceOut != null)
	    concordanceOut.close();
	concordanceOutWin.close();
	dataOut.close();
	unionout.close();

	// extract reliable/unreliable regions
	File reliableOut = new File(outdir, prefix + "RELIABLE-above" + reliableMin + ".bed");
	WigTools.wig2bed(concordanceOutWinFile, reliableMin, null, "Con_", reliableOut);
	stats.setString("ReliableRegionsFile", reliableOut.toString());

	File unreliableOut = new File(outdir, prefix + "UNRELIABLE-below" + unreliableMax + ".bed");
	WigTools.wig2bed(concordanceOutWinFile, null, unreliableMax, "Dis_", unreliableOut);
	stats.setString("UnreliableRegionsFile", unreliableOut.toString());

	// delete (large) wig files?
	if (!createWigs) {
	    concordanceOutWinFile.delete();
	}

	return stats;
    }

    /**
     * Convert a score into a filter string
     * 
     * @param score
     * @param reliableMin
     * @param unreliableMax
     * @return
     */
    private static String scoreToFilter(double score, double reliableMin, double unreliableMax) {
	if (score >= reliableMin) {
	    return (".");
	} else if (score < unreliableMax) {
	    return (FILTER_STRING_UNRELIABLE);
	} else {
	    return (FILTER_STRING_NO_CLASS);
	}
    }

    /**
     * Convert a filter string into a context character.
     * 
     * @param filterString
     * @return
     */
    private static String filterToContextChar(String filterString) {
	if (filterString == null)
	    return "?";
	if (filterString.equals("."))
	    return "R";
	else if (filterString.equals(FILTER_STRING_UNRELIABLE))
	    return "U";
	else
	    return "?";
    }

    /**
     * Scoring function.
     * 
     * @param matchVariants
     *            all variants that match at the same position.
     * @param scoreConcordantHit
     *            weight w^c
     * @param scoreDiscordantHit
     *            weight w^d
     * @return
     */
    private static ScoreResults score(List<SimpleVCFVariant> matchVariants, int scoreConcordantHit, int scoreDiscordantHit) {
	int dim = 0;
	if (matchVariants.size() > 0)
	    dim = matchVariants.get(0).getGenotypes().size();
	ScoreResults scoreResults = new ScoreResults(dim);

	double score = 0d;
	double count = 0d;
	double countMdi = 0d;
	double countDP = 0d;
	for (SimpleVCFVariant v : matchVariants) {
	    if (v == null || v.getGenotypes() == null)
		throw new RuntimeException("Internal error. Matched variant or its genotypes were null: " + v);
	    // double nocall = 0;
	    int hom = 0;
	    int het = 0;
	    int other = 0;
	    for (int i = 0; i < v.getGenotypes().size(); i++) {
		GENOTYPE g = v.getGenotypes().get(i);
		if (g == GENOTYPE.HETEROZYGOUS) {
		    het++;
		    scoreResults.CP[i]++;
		} else if (g == GENOTYPE.HOMOZYGOUS) {
		    hom++;
		    scoreResults.CP[i]++;
		} else
		    other++;
		// if (g == GENOTYPE.UNKNOWN)
		// nocall++;
	    }

	    // positive score only if genotypes match
	    if (Math.max(het, hom) == v.getGenotypes().size()) {
		score += scoreConcordantHit;
		count += Math.abs(scoreConcordantHit);
	    } else {
		score += scoreDiscordantHit;
		count += Math.abs(scoreDiscordantHit);
	    }

	    // increase overall het/hom/other count
	    scoreResults.calledHet += het;
	    scoreResults.calledHom += hom;
	    scoreResults.calledOther += other;

	    // add variant qualities
	    String[] tmp = v.getInfoField("QUALS").split(",");
	    for (int i = 0; i < dim; i++)
		scoreResults.AVGQUAL[i] += Float.parseFloat(tmp[i]);

	    // add minimum distance to indel (if any)
	    if (v.getInfoField("MDI") != null) {
		if (scoreResults.AVGMDI == null)
		    scoreResults.AVGMDI = 0d;
		scoreResults.AVGMDI += Integer.parseInt(v.getInfoField("MDI"));
		countMdi++;
	    }

	    // add raw read depth
	    if (v.getInfoField("DPRAW") != null) {
		if (scoreResults.AVGDPRAW == null)
		    scoreResults.AVGDPRAW = 0d;
		scoreResults.AVGDPRAW += Integer.parseInt(v.getInfoField("DPRAW"));
		countDP++;
	    }
	}

	scoreResults.score = score / count;
	for (int i = 0; i < dim; i++) {
	    scoreResults.CP[i] /= (double) matchVariants.size();
	    scoreResults.AVGQUAL[i] /= (double) matchVariants.size();
	}
	if (scoreResults.AVGMDI != null)
	    scoreResults.AVGMDI /= countMdi;

	if (scoreResults.AVGDPRAW != null)
	    scoreResults.AVGDPRAW /= countDP;

	return scoreResults;
    }

    /**
     * Helper class for storing scoring function results.
     * 
     * @author niko
     * 
     */
    private static class ScoreResults {

	public ScoreResults(int dim) {
	    this.dim = dim;
	    CP = new double[dim];
	    AVGQUAL = new double[dim];
	}

	public int dim;
	public double score;
	public double[] CP;
	public double[] AVGQUAL;
	// avg read depth from all samples
	public Double AVGDPRAW = null;
	// average MDI from all samples
	public Double AVGMDI = null;
	public int calledHom;
	public int calledHet;
	public int calledOther;

	public double getHomFrac() {
	    return (double) calledHom / (double) dim;
	}

	public double getHetFrac() {
	    return (double) calledHet / (double) dim;
	}

	public double getOtherFrac() {
	    return (double) calledOther / (double) dim;
	}
    }

    /**
     * Simple performance measuring method used in the evaluation.
     * 
     * @param union
     * @param reliable
     * @param unreliable
     * @param out
     * @throws IOException
     */
    static Statistics measure(List<String> goldStandardJoinFiles, File reliableF, File unreliableF, File outF) throws IOException {
	if (reliableF == null && unreliableF == null)
	    throw new IOException("A reliable and/or an unreliable region file needs to be provided!;");
	PrintStream out = null;
	if (outF != null)
	    out = new PrintStream(outF);
	GenomicITree reliable = reliableF != null ? (new SimpleBEDFile(reliableF)).getGenomicITree() : null;
	GenomicITree unreliable = unreliableF != null ? (new SimpleBEDFile(unreliableF)).getGenomicITree() : null;

	Statistics stats = new Statistics();
	stats.setString("GoldStandardJoinFiles", StringUtils.concat(goldStandardJoinFiles, ","));
	stats.setString("reliableF", (reliableF != null ? reliableF.toString() : "NA"));
	stats.setString("unreliableF", (unreliableF != null ? unreliableF.toString() : "NA"));

	for (String f : goldStandardJoinFiles) {
	    VCFIterator vi = new VCFIterator(new File(f));
	    while (vi.hasNext()) {
		SimpleVCFVariant v = vi.next();

		boolean inReliable = false;
		if (reliable != null) {
		    SortedSet<? extends GenomicInterval> res = reliable.query(v.getGenomicPosition());
		    inReliable = (res != null && res.size() > 0);
		}

		boolean inUnreliable = false;
		if (unreliable != null) {
		    SortedSet<? extends GenomicInterval> res = unreliable.query(v.getGenomicPosition());
		    inUnreliable = (res != null && res.size() > 0);
		}

		if (reliable == null && unreliable != null)
		    inReliable = !inUnreliable;
		else if (reliable != null && unreliable == null)
		    inUnreliable = !inReliable;
		else if (reliable == null && unreliable == null) {
		    if (out != null)
			out.close();
		    throw new IOException("internal error.");
		}

		String fil = v.getFilter();

		if (fil.equals(FILTER_STRING_UNRELIABLE)) {
		    stats.inc("ALL");

		    // variant should be classified as unreliable
		    if (inReliable) {
			// FALSE-classify-unreliable-as-reliable
			stats.inc("FP");
			stats.setString("EXAMPLE FP", f + "," + v.getGenomicPosition().toString1basedCanonicalHuman());
		    } else if (inUnreliable) {
			// TRUE-classify-unreliable-as-unreliable
			stats.inc("TN");
			stats.setString("EXAMPLE TN", f + "," + v.getGenomicPosition().toString1basedCanonicalHuman());
		    } else {
			// NEUTRAL-no-classify-unreliable
			stats.inc("UCA");
			stats.setString("EXAMPLE unclassified unreliable", f + "," + v.getGenomicPosition().toString1basedCanonicalHuman());
		    }
		} else if (fil.equals(".") || fil.equals("PASS")) {
		    stats.inc("ALL");

		    // variant should be classified as reliable
		    if (inReliable) {
			// TRUE-classify-reliable-as-reliable
			stats.inc("TP");
			stats.setString("EXAMPLE TP", f + "," + v.getGenomicPosition().toString1basedCanonicalHuman());
		    } else if (inUnreliable) {
			// FALSE-classify-reliable-as-unreliable
			stats.inc("FN");
			stats.setString("EXAMPLE FN", f + "," + v.getGenomicPosition().toString1basedCanonicalHuman());
		    } else {
			// NEUTRAL-no-classify-reliable
			stats.inc("UCB");
			stats.setString("EXAMPLE unclassified reliable", f + "," + v.getGenomicPosition().toString1basedCanonicalHuman());
		    }
		} else {
		    // ignore entries that are unclassified in the gold standard
		    stats.inc("Ignored");
		}
	    }
	}

	// ensure counts
	stats.set("ALL", stats.get("ALL", 0));
	stats.set("TP", stats.get("TP", 0));
	stats.set("FP", stats.get("FP", 0));
	stats.set("TN", stats.get("TN", 0));
	stats.set("FN", stats.get("FN", 0));
	stats.set("UCA", stats.get("UCA", 0));
	stats.set("UCB", stats.get("UCB", 0));
	double ALL = stats.get("ALL");
	double TP = stats.get("TP");
	double FP = stats.get("FP");
	double TN = stats.get("TN");
	double FN = stats.get("FN");
	double UCA = stats.get("UCA");
	double UCB = stats.get("UCB");

	// calc measures
	stats.set("Unreliable", (FP + TN + UCA));
	stats.set("Reliable", (TP + FN + UCB));
	stats.set("Precision", TP / (TP + FP)); // aka positive predictive value
						// (PPV)
	stats.set("Recall", TP / (TP + FN)); // aka sensitivity aka true
					     // positive rate TPR
	stats.set("Specificity", TN / (TN + FP)); // aka true negative rate TPR
	stats.set("Accuracy", (TP + TN) / (TP + TN + FP + FN));
	stats.set("NPV", TN / (TN + FN));
	stats.set("FDR", FP / (FP + TP));
	stats.set("FNR", FN / (FN + FP));
	stats.set("F1", (2 * TP) / (2 * TP + FP + FN));

	stats.set("Classified", (TP + TN + FP + FN) / (ALL));
	stats.set("Unclassified", (UCA + UCB) / (ALL));
	if (out != null) {
	    stats.toFile(out);
	    out.close();
	}
	return stats;
    }

    /**
     * Evaluation function used in evaluation experiment 1.
     * 
     * @param datasetsFile
     * @param outDir
     * @param tsc
     * @param rep
     * @throws IOException
     */
    static void evaluate(File datasetsFile, File outDir, int tsc, int rep, int scoreConcordantHit, int scoreDiscordantHit, boolean excludeNonPass,
	    boolean dropAllFiltered, int winSize, double reliableMin, double unreliableMax, double smallWinSizeFactor, String smallWindowSizeRegionsF,
	    boolean useExactUnreliableClassification) throws IOException {

	boolean dontCheckSort = true;
	boolean createWigs = false;

	TabIterator ti = new TabIterator(datasetsFile, "#");
	List<String> datasets = new ArrayList<>();
	while (ti.hasNext()) {
	    String[] t = ti.next();
	    datasets.add(t[0]);
	}

	if (!outDir.exists())
	    throw new FileNotFoundException("Outdir " + outDir + " does not exist.");

	// random shuffle of datasets
	Collections.shuffle(datasets);
	List<String> train = new ArrayList<>();
	List<String> goldStandard = new ArrayList<>();
	for (int i = 0; i < datasets.size(); i++)
	    if (i < tsc)
		train.add(datasets.get(i));
	    else
		goldStandard.add(datasets.get(i));

	System.out.println("------------ TSC: " + tsc + ", REP: " + rep + "------------");
	System.out.println("training size: " + train.size());
	System.out.println("gold standard size: " + goldStandard.size());

	File outDirTrain = new File(outDir, "Train-tsc-" + tsc + "-rep" + rep);
	Statistics trainStats = BACKUPCalcreliabilitySignals.calc(train.toArray(new String[train.size()]), null, null, null, scoreConcordantHit, scoreDiscordantHit,
		dontCheckSort, excludeNonPass, dropAllFiltered, createWigs, winSize, reliableMin, unreliableMax, smallWinSizeFactor, smallWindowSizeRegionsF,
		outDirTrain);
	File classifiedReliable = new File(trainStats.getString("ReliableRegionsFile"));
	File classifiedUnreliable = null;
	if (useExactUnreliableClassification)
	    classifiedUnreliable = new File(trainStats.getString("UnreliableRegionsFile"));

	// MEASURE
	Statistics results = BACKUPCalcreliabilitySignals.measure(goldStandard, classifiedReliable, classifiedUnreliable, null);

	File goldStandardResults = new File(outDir, "Results-tsc-" + tsc + "-rep" + rep);
	PrintStream out = new PrintStream(goldStandardResults);
	// r.toFileAsDataRow(out, true);

	String[] columns = new String[] { "Accuracy", "Classified", "Unclassified", "Precision", "Recall", "Specificity", "ALL", "TP", "FP", "TN", "FN", "UCA",
		"UCB", "Unreliable", "Reliable", "EXAMPLE TP", "EXAMPLE FP", "EXAMPLE TN", "EXAMPLE FN", "EXAMPLE unclassified unreliable",
		"EXAMPLE unclassified reliable" };
	StringBuffer header = new StringBuffer();
	StringBuffer body = new StringBuffer();

	boolean first = true;
	for (String k : columns) {
	    if (first)
		first = false;
	    else {
		header.append("\t");
		body.append("\t");
	    }
	    header.append(k);
	    String d = results.get(k) == null ? null : results.get(k) + "";
	    if (d == null)
		d = results.getString(k) == null ? "NA" : results.getString(k);
	    body.append(d);
	}

	out.println(header);
	out.println(body);
	out.close();
    }

    /**
     * Join method.
     * 
     * This will call a position concordant iff the variant is called in all
     * input VCFs and if their GT fields match.
     * 
     * @see test data at src/test/resources/RG
     * 
     * @param vcfs
     * @param scopes
     * @param testedChromsArray
     * @param dontCheckSort
     * @param prefix
     * @throws IOException
     */
    static void join(String[] vcfFiles, String[] labels, String[] excludeFiles, String[] includeFiles, String[] testedChromsArray, boolean dontCheckSort,
	    boolean excludeNonPass, boolean dropAllFiltered, boolean dropAnyFiltered, int indelMergeWin, File outFileSNV, File outFileINDEL) throws IOException {

	long start = System.currentTimeMillis();
	Statistics stats = new Statistics();
	stats.setString("Num of Input VCF Files", vcfFiles.length + "");
	stats.setString("Input VCF Files", Arrays.toString(vcfFiles));
	stats.setString("Input VCF Labels", Arrays.toString(labels));
	stats.setString("Excluded Regions", Arrays.toString(excludeFiles));
	stats.setString("included Regions", Arrays.toString(includeFiles));
	stats.setString("Tested Chromosomes", Arrays.toString(testedChromsArray));
	stats.setString("ExcludeNonPass", excludeNonPass + "");
	stats.setString("DropAllFiltered", dropAllFiltered + "");
	stats.setString("DropAnyFiltered", dropAnyFiltered + "");
	stats.setString("DontCheckSort", dontCheckSort + "");
	stats.setString("indelMergeWin", indelMergeWin + "");
	stats.setString("OutfileSNV", outFileSNV.toString());
	stats.setString("OutfileINDEL", outFileINDEL.toString());
	if (debug)
	    System.out.println(stats);

	if (labels == null) {
	    labels = vcfFiles;
	} else {
	    if (labels.length != vcfFiles.length)
		throw new IOException("Number of vcffiles and labels have to match! ");
	}

	// prepare chromosome scope.
	Set<String> testedChroms = null;
	if (testedChromsArray != null) {
	    testedChroms = new HashSet<String>();
	    for (String tc : testedChromsArray)
		testedChroms.add(tc);
	}

	// sorted list of current pointers
	SortedSet<GenomicPosition> pointers = new TreeSet<GenomicPosition>();

	// prepare include/exclude BEDs
	List<AbstractBedIterator> included = new ArrayList<AbstractBedIterator>();
	if (includeFiles != null)
	    for (String inc : includeFiles) {
		if (inc.contains(":")) {
		    MemoryBedIterator it = new MemoryBedIterator();
		    it.setName("includes");
		    it.add(inc);
		    GenomicInterval firstInterval = incBedIt(it, testedChroms);
		    if (firstInterval != null) { // there is at least one
			// interval
			included.add(it);
		    }
		} else {
		    BedIterator it = new BedIterator(new File(inc));
		    it.setName("includes");
		    GenomicInterval firstInterval = incBedIt(it, testedChroms);
		    if (firstInterval != null) { // there is at least one
			// interval
			included.add(it);
		    }
		}
	    }
	List<AbstractBedIterator> excluded = new ArrayList<AbstractBedIterator>();
	if (excludeFiles != null)
	    for (String ex : excludeFiles) {
		if (ex.equals("ALL")) {
		    excluded.clear();
		    excluded.add(null);
		    break;
		} else if (ex.contains(":")) {
		    MemoryBedIterator it = new MemoryBedIterator();
		    it.setName("excludes");
		    it.add(ex);
		    GenomicInterval firstInterval = incBedIt(it, testedChroms);
		    if (firstInterval != null) { // there is at least one
			// interval
			excluded.add(it);
		    }
		} else {
		    BedIterator it = new BedIterator(new File(ex));
		    it.setName("excludes");
		    GenomicInterval firstInterval = incBedIt(it, testedChroms);
		    if (firstInterval != null) { // there is at least one
			// interval
			excluded.add(it);
		    }
		}
	    }

	// prepare iterators
	List<VCFIterator> vits = new ArrayList<VCFIterator>();
	for (String vcffn : vcfFiles) {

	    String sample = null;
	    if (vcffn.contains("@")) {
		String[] tmp = vcffn.split("@");
		sample = tmp[0];
		vcffn = tmp[1];
	    }
	    File vcf = new File(vcffn);
	    // check VCF sorting?
	    if (!dontCheckSort)
		if (!SimpleVCFFile.checkSortedVCF(vcf, testedChroms))
		    throw new IOException("File " + vcf + " is not sorted by coords");
	    VCFIterator it = null;
	    if (sample == null)
		it = new VCFIterator(vcf);
	    else
		it = new VCFSampleIterator(vcf, sample);

	    SimpleVCFVariant v = incVcfIt(it, testedChroms, excludeNonPass);
	    if (v != null)
		pointers.add(v.getGenomicPosition());
	    vits.add(it);
	}

	if (pointers.size() == 0)
	    throw new IOException("All iterators were empty - do chrom names match the tested ones (" + testedChroms
		    + ") ? Does the file contain GT information!?");

	// create output streams
	@SuppressWarnings("resource")
	PrintStream joinOut = new PrintStream(outFileSNV);
	StringBuffer buf = new StringBuffer();
	buf.append("##fileformat=VCFv4.1\n");
	buf.append("##fileDate=" + new Date() + "\n");
	buf.append("##source=" + BACKUPCalcreliabilitySignals.class + "; SNVs only\n");
	buf.append("##INFO=<ID=QUALS,Number=.,Type=Float,Description=\"Comma-separated list of the QUAL fields of the individual callers. QUAL "
		+ "is the Phred-scaled quality score for the assertion made in ALT. i.e. -10log10 prob (call in ALT is wrong). "
		+ "If ALT is \'.\' (no variant) then this is -10log10 prob(variant), and if ALT is not \'.\' this is -10 log10 prob(no variant)."
		+ "If unknown, the missing value should be specified.\">\n");
	buf.append("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Maximum found read-depth for this variant.\">\n");
	buf.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");

	buf.append("##FILTER=<ID=PASS,Description=\"A reliable position\">\n");
	buf.append("##FILTER=<ID=" + FILTER_STRING_UNRELIABLE + ",Description=\"An unreliable position\">\n");
	buf.append("##FILTER=<ID=" + FILTER_STRING_NO_CLASS + ",Description=\"An unclassified position\">\n");

	buf.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (String l : labels) {
	    buf.append("\t" + l);
	}
	buf.append("\n");
	PrintStream joinOutIndels = new PrintStream(outFileINDEL);
	StringBuffer bufIndel = new StringBuffer();
	bufIndel.append("##fileformat=VCFv4.1\n");
	bufIndel.append("##fileDate=" + new Date() + "\n");
	bufIndel.append("##source=" + BACKUPCalcreliabilitySignals.class + "; INDELs only\n");
	bufIndel.append("##INFO=<ID=QUALS,Number=.,Type=Float,Description=\"Comma-separated list of the QUAL fields of the individual callers. QUAL "
		+ "is the Phred-scaled quality score for the assertion made in ALT. i.e. -10log10 prob (call in ALT is wrong). "
		+ "If ALT is \'.\' (no variant) then this is -10log10 prob(variant), and if ALT is not \'.\' this is -10 log10 prob(no variant). "
		+ "If unknown, the missing value should be specified.\">\n");
	bufIndel.append("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Maximum found read-depth for this variant (usually the samtools raw reads field!).\">\n");
	bufIndel.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");

	bufIndel.append("##FILTER=<ID=PASS,Description=\"A reliable position\">\n");
	bufIndel.append("##FILTER=<ID=" + FILTER_STRING_UNRELIABLE + ",Description=\"An unreliable position\">\n");
	bufIndel.append("##FILTER=<ID=" + FILTER_STRING_NO_CLASS + ",Description=\"An unclassified position\">\n");

	bufIndel.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (String l : labels) {
	    bufIndel.append("\t" + l);
	}
	bufIndel.append("\n");

	int count = 0;

	// for storing indels.
	List<GenomicInterval> indelList = new ArrayList<GenomicInterval>();

	// loop over variants
	do {
	    // filled with all matches at pos <pointers.first()>
	    List<SimpleVCFVariant> matchVariants = new ArrayList<SimpleVCFVariant>();
	    int filteredVariantsCount = 0; // keeps track of how many filtered
	    Boolean allINDELs = null;

	    // check iterators
	    for (int i = 0; i < vits.size(); i++) {
		VCFIterator it = vits.get(i);
		if (it.getCurrentVariant().getGenomicPosition().equals(pointers.first())) {
		    // match at current position
		    it.getCurrentVariant().addField("vcf-id", i + "");
		    matchVariants.add(it.getCurrentVariant());
		    if (it.getCurrentVariant().isFiltered())
			filteredVariantsCount++;
		    if (it.getCurrentVariant().isINDEL()) {
			if (allINDELs == null)
			    allINDELs = true;
		    } else
			allINDELs = false;
		    if (it.hasNext()) {
			// handle multiple variants at same position!
			SimpleVCFVariant var = incVcfIt(it, testedChroms, excludeNonPass);
			if (var != null) {
			    GenomicPosition nxt = var.getGenomicPosition();
			    while (nxt.equals(pointers.first())) {
				var = incVcfIt(it, testedChroms, excludeNonPass);
				if (var == null) {
				    break;
				}
				nxt = var.getGenomicPosition();
			    }
			    // add next position
			    pointers.add(nxt);
			}
		    }

		}
	    }

	    boolean positionIncluded = checkPositionIncluded(pointers.first(), testedChroms, included, excluded);

	    // drop matches where any variant was filtered
	    if (dropAnyFiltered && filteredVariantsCount > 0)
		positionIncluded = false;

	    // drop matches where all variants are filtered (applies only if
	    // filtered vars are included!)
	    if (dropAllFiltered && matchVariants.size() == filteredVariantsCount)
		positionIncluded = false;

	    // handle INDEL positions!
	    if (allINDELs != null && allINDELs) {
		positionIncluded = false;
		for (SimpleVCFVariant x : matchVariants) {
		    GenomicInterval g = new GenomicInterval(x.getChromosome(), x.getPosition() - indelMergeWin, x.getPosition() + x.getWidth() + indelMergeWin,
			    null);
		    g.setOriginalChrom(x.getChromosomeOrig());
		    g.setAnnotation("var", x);
		    indelList.add(g);
		}
	    }

	    if (positionIncluded) {
		SimpleVCFVariant v = new SimpleVCFVariant();

		v.setChromosomeOrig(matchVariants.get(0).getChromosomeOrig());
		v.setPosition(matchVariants.get(0).getPosition());
		v.setID(matchVariants.get(0).getID());
		Set<String> refAlleles = new HashSet<String>();
		Set<String> altAlleles = new HashSet<String>();
		for (SimpleVCFVariant x : matchVariants) {
		    refAlleles.add(x.getRefString());
		    for (String a : x.getAltAlleles())
			altAlleles.add(a);
		}
		boolean ignoreVariant = false;
		if (refAlleles.size() > 1) {
		    stats.inc("Ignored-Due-To-Ref-Mismatch");
		    ignoreVariant = true;
		}

		if (!ignoreVariant) {
		    v.setRefString(StringUtils.concat(refAlleles, ","));
		    v.setAltString(StringUtils.concat(altAlleles, ","));
		    List<Float> QUAL = new ArrayList<Float>();
		    Integer maxEstimatedDepth = 0;

		    Map<String, List<String>> gtdata = new HashMap<>();
		    List<String> gts = new ArrayList<>();
		    // get genotypes
		    String countMatchGT = null;
		    for (int i = 0; i < vits.size(); i++) {
			String GT = "./.";
			Float qual = 0f;
			for (SimpleVCFVariant x : matchVariants) {
			    String idstr = x.getInfoField("vcf-id");
			    if (idstr != null && Integer.parseInt(idstr) == i) {
				// get quality
				qual = x.getQuality();

				// get max depth
				Integer dp = x.estimateCoverage();
				if (dp != null) {
				    if (maxEstimatedDepth == null)
					maxEstimatedDepth = dp;
				    else
					maxEstimatedDepth = Math.max(dp, maxEstimatedDepth);
				}

				// get genotype
				if (x.getGenotypes() == null || x.getGenotypes().size() != 1) {
				    throw new IOException("Invalid genotype spec at " + x);
				}
				GENOTYPE g = x.getGenotypes().get(0);
				if (g == GENOTYPE.HETEROZYGOUS)
				    GT = "0/1";
				else if (g == GENOTYPE.HOMOZYGOUS)
				    GT = "1/1";
				else if (g == GENOTYPE.UNKNOWN)
				    GT = "./.";
				else
				    throw new IOException("Invalid genotype spec at " + x);
				break;
			    }
			}
			if (countMatchGT == null)
			    countMatchGT = GT;
			else if (!GT.equals(countMatchGT)) {
			    stats.inc("positions-with-discordant-GT");
			    v.setFilter(FILTER_STRING_UNRELIABLE);
			}
			gts.add(GT);
			QUAL.add(qual);
		    }
		    gtdata.put("GT", gts);
		    v.setGtdata(gtdata);
		    v.addField("QUALS", StringUtils.concat(QUAL, ","));
		    if (maxEstimatedDepth != null)
			v.addField("DP", maxEstimatedDepth + "");

		    buf.append(v.toStringAsIs() + "\n");
		    stats.inc("Emitted SNVs");
		    if (buf.length() > 100000) {
			joinOut.print(buf.toString());
			buf = new StringBuffer();
		    }
		}

	    } // included

	    // can we resolve some indels? its safe at chr-changes
	    if (indelList.size() > 0 && !indelList.get(indelList.size() - 1).getOriginalChrom().equals(pointers.first().getChromosomeOriginal())) {
		resolveIndels(indelList, bufIndel, stats, dropAllFiltered, dropAnyFiltered, joinOutIndels, vits);
		indelList.clear();
		if (bufIndel.length() > 100000) {
		    joinOutIndels.print(bufIndel.toString());
		    bufIndel = new StringBuffer();
		}
	    }

	    if (debug)
		if (++count % 10000 == 0) {
		    System.out.println("> " + pointers.first().toString1based() + ", INDELs: " + indelList.size());
		}

	    // go to next position
	    pointers.remove(pointers.first());

	} while (!pointers.isEmpty());
	joinOut.print(buf.toString());
	joinOut.close();

	resolveIndels(indelList, bufIndel, stats, dropAllFiltered, dropAnyFiltered, joinOutIndels, vits);

	joinOutIndels.print(bufIndel.toString());
	joinOutIndels.close();

	// write stats
	stats.setString("Execution time [ms]", (System.currentTimeMillis() - start) + "");
	stats.toFile(new PrintStream(new File(outFileSNV.getAbsolutePath() + ".stats")));

    }

    /**
     * Helper method for the join() method. Resolves the current list of INDELs.
     * 
     * @param indelList
     * @param bufIndel
     * @param stats
     * @param dropAllFiltered
     * @param joinOutIndels
     * @param vits
     * @throws IOException
     */
    private static void resolveIndels(List<GenomicInterval> indelList, StringBuffer bufIndel, Statistics stats, boolean dropAllFiltered,
	    boolean dropAnyFiltered, PrintStream joinOutIndels, List<VCFIterator> vits) throws IOException {

	if (indelList.size() == 0)
	    return;

	if (debug)
	    System.out.println("Resolve " + indelList.size() + " IDNELs such as " + indelList.get(0));

	// ----------------------------------------------------
	// resolve indels
	// ----------------------------------------------------
	Stack<GenomicInterval> stack = new Stack<>();
	for (int idx = 0; idx < indelList.size(); idx++) {
	    GenomicInterval gi = indelList.get(idx);
	    if (!stack.isEmpty() && (stack.lastElement().getRightPosition().compareTo(gi.getLeftPosition()) < 0 || idx == indelList.size() - 1)) {
		boolean positionIncluded = true;
		boolean allFiltered = true;
		boolean anyFiltered = false;

		Set<String> ids = new HashSet<>();
		List<SimpleVCFVariant> matchVariants = new ArrayList<SimpleVCFVariant>();
		for (GenomicInterval x : stack) {
		    SimpleVCFVariant v = (SimpleVCFVariant) x.getAnnotation("var");
		    ids.add(v.getInfoField("vcf-id"));
		    matchVariants.add(v);
		    if (!v.isFiltered())
			allFiltered = false;
		    else
			anyFiltered = true;
		}
		// System.out.println("IDs: " + ids.size());
		SimpleVCFVariant matchVariant = matchVariants.get(0);
		SimpleVCFVariant v = new SimpleVCFVariant();
		v.setChromosomeOrig(matchVariant.getChromosomeOrig());
		v.setPosition(matchVariant.getPosition());
		v.setID(matchVariant.getID());
		v.setRefString(matchVariant.getRefString());
		v.setAltString(matchVariant.getAltString());
		List<Float> QUAL = new ArrayList<Float>();
		Integer maxEstimatedDepth = 0;
		Map<String, List<String>> gtdata = new HashMap<>();
		List<String> gts = new ArrayList<>();
		// get genotypes
		String countMatchGT = null;
		Set<String> indelStartPositions = new HashSet<String>();
		for (int i = 0; i < vits.size(); i++) {
		    String GT = "./.";
		    Float qual = 0f;
		    for (SimpleVCFVariant x : matchVariants) {
			String idstr = x.getInfoField("vcf-id");
			if (idstr != null && Integer.parseInt(idstr) == i) {
			    // get quality
			    qual = x.getQuality();

			    // get max depth
			    Integer dp = x.estimateCoverage();
			    if (dp != null) {
				if (maxEstimatedDepth == null)
				    maxEstimatedDepth = dp;
				else
				    maxEstimatedDepth = Math.max(dp, maxEstimatedDepth);
			    }

			    // count number of different start positions
			    // (debugging only)
			    indelStartPositions.add(x.getGenomicPosition().toString());

			    // get GT
			    if (x.getGenotypes() == null || x.getGenotypes().size() != 1) {
				throw new IOException("Invalid genotype spec at " + x);
			    }
			    GENOTYPE g = x.getGenotypes().get(0);
			    if (g == GENOTYPE.HETEROZYGOUS)
				GT = "0/1";
			    else if (g == GENOTYPE.HOMOZYGOUS)
				GT = "1/1";
			    else if (g == GENOTYPE.UNKNOWN)
				GT = "./.";
			    else
				throw new IOException("Invalid genotype spec at " + x);
			    break;
			}
		    }
		    v.addField("matchedINDELs", matchVariants.size() + "");
		    v.addField("indelStartPositions", indelStartPositions.size() + "");

		    if (countMatchGT == null)
			countMatchGT = GT;
		    else if (!GT.equals(countMatchGT)) {
			stats.inc("positions-with-discordant-GT");
			v.setFilter(FILTER_STRING_UNRELIABLE);
		    }
		    gts.add(GT);
		    QUAL.add(qual);
		}
		gtdata.put("GT", gts);
		v.setGtdata(gtdata);
		v.addField("QUALS", StringUtils.concat(QUAL, ","));
		if (maxEstimatedDepth != null)
		    v.addField("DP", maxEstimatedDepth + "");

		// drop matches where any variant was filtered
		if (dropAnyFiltered && anyFiltered)
		    positionIncluded = false;

		// drop matches where all variants are filtered (applies only if
		// filtered vars are included!)
		if (dropAllFiltered && allFiltered)
		    positionIncluded = false;

		if (positionIncluded) {
		    stats.inc("Emitted INDELs");
		    bufIndel.append(v.toStringAsIs() + "\n");
		}
		stack.clear();
	    }
	    stack.push(gi);
	}
    }

    /**
     * Calculates a histogram of windowed unreliability density values.
     * 
     * @param file
     * @throws IOException
     * @throws IncompatibleIntervalsException
     */
    private static void calcDensityHistogram(File vcf, int winSize, int cutoff, File bedOut) throws IOException, IncompatibleIntervalsException {
	VCFContextWinIterator it = new VCFContextWinIterator(vcf, winSize);
	Histogram<Integer> hist = new Histogram<>("filtered/unfiltered variant changes");

	int intcount = 0;
	int c = 0;
	GenomicITree gi = new GenomicITree(null);
	while (it.hasNext()) {
	    LinkedList<SimpleVCFVariant> ctx = it.next();
	    SimpleVCFVariant v = it.getCurrentVariant();

	    c = 0;
	    boolean wasFiltered = false;
	    // count # of filter/non-filter changes in context
	    for (SimpleVCFVariant x : ctx) {
		if (x.isFiltered()) {
		    if (!wasFiltered) {
			c++;
			wasFiltered = true;
		    }
		} else {
		    if (!x.isFiltered()) {
			if (wasFiltered) {
			    c++;
			    wasFiltered = false;
			}
		    }
		}
	    }

	    if (c >= cutoff) {
		gi.insert(new GenomicInterval(v.getChromosomeOrig(), v.getPosition() - winSize, v.getPosition() + winSize, ("int" + (++intcount) + "_c" + c)));
	    }
	    hist.push(c);
	}
	gi.buildTree();

	// create bed
	PrintStream out = new PrintStream(bedOut);
	for (GenomicInterval i : gi.getMergedSortedIntervals()) {
	    out.println(i.getOriginalChrom() + "\t" + i.getMin() + "\t" + i.getMax() + "\tmerged_" + StringUtils.countMatches(i.getUri(), '+'));
	}
	out.close();

	// create stats
	out = new PrintStream(bedOut + ".hist");
	out.println("# Found " + intcount + " intervals for cutoff " + cutoff);
	out.println(hist);
	out.close();

    }

    /**
     * Print usage information.
     * 
     * @param options
     */
    private static void usage(Options options, String subcommand, String e) {

	if (subcommand == null) {
	    System.out.println("Usage:\t\tjava -jar x.jar " + CMD + " <command> [options]:\t");
	    System.out.println();
	    System.out.println("Command:\tjoin\tJoins multiple VCFs from different VC pipelines");
	    System.out.println("Command:\tcalc\tCalculates a genomic partition from joined VCF files");
	    System.out.println("Command:\tmeasure\tMeasures statistics from a union vcf");
	    System.out.println("Command:\teval\tEvaluates a genomic partition by using a gold-standard union vcf");
	    System.out.println();
	} else {

	    HelpFormatter hf = new HelpFormatter();
	    hf.setLeftPadding(10);
	    hf.setDescPadding(2);
	    hf.setWidth(160);
	    hf.setSyntaxPrefix("Usage:    ");
	    hf.printHelp("java -jar x.jar " + CMD + " " + subcommand, "Params:", options, "", true);
	}
	if (e != null) {
	    System.out.println("\nError: " + e);
	}
	System.exit(1);
    }

    public static void main(String[] args) {

	// args = new String[] { "calcDensityHistogram",
	//
	// "-i", "/temp/Niko/RG/test.vcf",
	//
	// "-w", "1000",
	//
	// "-c", "10",
	//
	// "-o", "/temp/Niko/RG/test.vcf.HIGHDENS.bed" };

	// args = new String[] { "calc", "-o", "c:/data/reliableGenome/calc",
	// "-w", "1000", "-scoringSchema", "1,-3", "-thresholds", "0.5,0.5",
	// "-dontCheckSort",
	// "-dropAllFiltered", "-v", "-d",
	// "c:/data/reliableGenome/vcf/joinSNVs.vcf" };

	// args = new String[] { "join",
	//
	// "-d", "src/test/resources/RG/C1.vcf", "-dl", "C1",
	//
	// "-d", "src/test/resources/RG/C2.vcf", "-dl", "C2",
	//
	// "-d", "src/test/resources/RG/C3.vcf", "-dl", "C3",
	//
	// "-o", "src/test/resources/RG/C1+2+3.vcf",
	//
	// "-oi", "src/test/resources/RG/C1+2+3.INDEL.vcf",
	//
	// "-indelMergeWin", "3",
	//
	// "-dropAllFiltered",
	//
	// "-dontCheckSort",
	//
	// "-v"
	//
	// };

	// String[] extensions = new String[] { "vcf.gz" };
	// List<File> files = (List<File>) FileUtils.listFiles(new
	// File("/temp/Niko/RG/test/"), extensions, false);
	// List<String> argsl = new ArrayList<String>();
	// argsl.add("calc");
	// for (File f : files) {
	// argsl.add("-d");
	// argsl.add(f.toString());
	// }
	// argsl.add("-o");
	// argsl.add("/temp/Niko/RG/test/test0");
	// argsl.add("-v");
	// argsl.add("-dontCheckSort");
	// argsl.add("-w");
	// argsl.add("1000");
	// args = argsl.toArray(new String[argsl.size()]);

	CommandLineParser parser = new PosixParser();

	// create the Options
	Options options = new Options();
	options.addOption("h", "help", false, "Print this usage information.");
	String subcommand = null;

	try {
	    // parse the command line arguments
	    CommandLine line = parser.parse(options, args, true);

	    // validate that block-size has been set
	    if (line.getArgs().length == 0 || line.hasOption("h")) {
		usage(options, subcommand, null);
	    }

	    subcommand = line.getArgs()[0];

	    if (subcommand.equalsIgnoreCase("join")) {
		Option o = new Option("d", "data", true,
			"Input data VCF files. You may use the Sample@VCFFile notation to iterate over individual samples in a VCF file.");
		o.setRequired(true);
		options.addOption(o);

		o = new Option("dl", "datalabel", true,
			"Lables for the analyzed data samples. Has to match the number of input vcf files/samples. If omitted, filenames are used.");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("e", "excluded", true, "Excluded regions file(s): (gziped) BED(s). Use ALL to exclude all.");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("i", "included", true, "Included regions file(s): (gziped) BED(s).");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("chrom", true, "Tested (compared) chromosomes (comma-separated). ");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("excludeNonPass", false, "Do not included filtered SNPs. Use with care.");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("dropAllFiltered", false,
			"If set, all matching positions where *all* variants are filtered will be ignored. Only useful w/o excludeNonPass.");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("dropAnyFiltered", false,
			"If set, all matching positions where *any* variants is filtered will be ignored. Only useful w/o excludeNonPass.");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("dontCheckSort", false, "Do not check whether files are sorted. Use with care.");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("indelMergeWin", true,
			"INDEL intervals will be extended by this number of nts up/downstream and overlapping INDELs will be merged! (default: 0bp).");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("o", "outSNV", true, "Output file for SNVs.");
		o.setRequired(true);
		options.addOption(o);

		o = new Option("oi", "outINDEL", true, "Output file for INDELs.");
		o.setRequired(true);
		options.addOption(o);

		options.addOption("v", "verbose", false, "be verbose.");

		line = parser.parse(options, args);
		if (line.hasOption("v"))
		    debug = true;
		else
		    debug = false;

		if (line.hasOption("excludeNonPass") && line.hasOption("dropAllFiltered"))
		    System.err.println("Warning: dropAllFilter cannot be applied when non-pass variants are excluded.");

		join(line.getOptionValues("d"), line.getOptionValues("dl"), line.getOptionValues("e"), line.getOptionValues("i"),
			line.hasOption("chrom") ? line.getOptionValues("chrom") : CanonicalChromsomeComparator.canonicalChromosomes,
			line.hasOption("dontCheckSort"), line.hasOption("excludeNonPass"), line.hasOption("dropAllFiltered"),
			line.hasOption("dropAnyFiltered"), line.hasOption("indelMergeWin") ? Integer.parseInt(line.getOptionValue("indelMergeWin")) : 0,
			new File(line.getOptionValue("o")), new File(line.getOptionValue("oi")));

		System.err.println("Finished.");

		System.exit(0);

	    } else if (subcommand.equalsIgnoreCase("eval")) {
		Option o = new Option("i", true, "Input dataset file.");
		o.setRequired(true);
		options.addOption(o);

		o = new Option("o", "out", true, "Output file.");
		o.setRequired(true);
		options.addOption(o);

		o = new Option("tsc", "trainSetCount", true,
			"Number of (random) datasets that will be used for classifier training. All remaining datasets will form the 'gold standard' group");
		o.setRequired(true);
		options.addOption(o);

		o = new Option("rep", "iteration", true, "Iteration count. Will solely be used in the resulting output filenames.");
		o.setRequired(true);
		options.addOption(o);

		o = new Option("excludeNonPass", false, "Do not included filtered SNPs. Use with care.");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("dropAllFiltered", false,
			"If set, all matching positions where *all* variants are filtered will be ignored. Only useful w/o excludeNonPass.");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("s", "scoringSchema", true,
			"Scoring schema consisting of two comma-separated weights, a positive one for concordant calls and a negative one for discordant calls. (def: '1,-3').");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("t", "thresholds", true, "Tresholds for reliable/unreliable regions (def: '0.5,-0.001').");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("w", "win", true, "Window size in bp (def: 1000bp).");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("swf", "smallWindowFactor", true,
			"The configured window size will be multiplied by this factor to calculate the smaller window size (def: 0.1).");
		o.setRequired(false);
		options.addOption(o);

		o = new Option(
			"swr",
			"smallWindowFactorRegions",
			true,
			"Optional BED file that contains regions that are considered rather unreliable (e.g., low-complexity regions). "
				+ "If configured, the signal interpolator will switch to a smaller window size (see smallWindowFactor) in such regions (def: N/A).");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("useExactUnreliableClassification", false, "If unset, all positions that are not reliable will be called unreliable (default).");
		o.setRequired(false);
		options.addOption(o);

		line = parser.parse(options, args);
		if (line.hasOption("v"))
		    debug = true;
		else
		    debug = false;

		File datasetsFile = new File(line.getOptionValue("i"));
		File outDir = new File(line.getOptionValue("o"));
		int tsc = Integer.parseInt(line.getOptionValue("tsc"));
		int rep = Integer.parseInt(line.getOptionValue("rep"));

		String scoringschema = line.hasOption("scoringSchema") ? line.getOptionValue("scoringSchema") : "1,-3";
		if (!scoringschema.contains(","))
		    throw new IOException("Invalid scoring schema. Cannot parse " + scoringschema);
		int scoreConcordantHit = Integer.parseInt(scoringschema.split(",")[0]);
		int scoreDiscordantHit = Integer.parseInt(scoringschema.split(",")[1]);
		if (!(scoreConcordantHit >= 0 && scoreDiscordantHit <= 0))
		    throw new IOException("Invalid scoring schema. weight for concordant calls must be positive, weight for discordant calls negative.");

		String thresholds = line.hasOption("thresholds") ? line.getOptionValue("thresholds") : "0.5,-0.001";
		if (!thresholds.contains(","))
		    throw new IOException("Invalid thresholds specification. Cannot parse " + thresholds);
		double reliableMin = Double.parseDouble(thresholds.split(",")[0]);
		double unreliableMax = Double.parseDouble(thresholds.split(",")[1]);
		if (reliableMin > 1 || reliableMin < -1 || unreliableMax > 1 || unreliableMax < -1 || reliableMin < unreliableMax)
		    throw new IOException("reliableMin & unreliableMax have to be from [-1, 1] and reliableMin>=unreliableMax");

		double smallWindowFactor = line.hasOption("smallWindowFactor") ? Double.parseDouble(line.getOptionValue("smallWindowFactor")) : 0.1d;

		evaluate(datasetsFile, outDir, tsc, rep, scoreConcordantHit, scoreDiscordantHit, line.hasOption("excludeNonPass"),
			line.hasOption("dropAllFiltered"), line.hasOption("win") ? Integer.parseInt(line.getOptionValue("win")) : 1000, reliableMin,
			unreliableMax, smallWindowFactor, line.getOptionValue("swr"), line.hasOption("useExactUnreliableClassification"));

		System.err.println("Finished.");

		System.exit(0);
	    } else if (subcommand.equalsIgnoreCase("measure")) {
		Option o = new Option("i", "union", true,
			"Input union VCF file(s) as computed by the calc method for some tested set of VCFs. This will be used as gold standard.");
		o.setRequired(true);
		options.addOption(o);

		o = new Option("r", "reliable", true,
			"Input BED file of reliable regions. If not set, everything that is not in an unreliable region will be considered reliable.");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("u", "unreliable", true,
			"Input BED file of unreliable regions. If not set, everything that is not in a reliable region will be considered unreliable.");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("o", "out", true, "Output file.");
		o.setRequired(true);
		options.addOption(o);

		options.addOption("v", "verbose", false, "be verbose.");

		line = parser.parse(options, args);
		if (line.hasOption("v"))
		    debug = true;
		else
		    debug = false;

		if (!line.hasOption("r") && !line.hasOption("u"))
		    throw new IOException("A reliable and/or an unreliable region file needs to be provided!");

		List<String> goldStandardJoinFiles = new ArrayList<String>();
		for (String f : line.getOptionValues("i"))
		    goldStandardJoinFiles.add(f);
		measure(goldStandardJoinFiles, line.hasOption("r") ? new File(line.getOptionValue("r")) : null,
			line.hasOption("u") ? new File(line.getOptionValue("u")) : null, new File(line.getOptionValue("o")));

		System.err.println("Finished.");

		System.exit(0);

	    } else if (subcommand.equalsIgnoreCase("calc")) {

		Option o = new Option("d", "data", true,
			"Input data VCF files. You may use the Sample@VCFFile notation to iterate over individual samples in a VCF file.");
		o.setRequired(true);
		options.addOption(o);

		o = new Option("e", "excluded", true, "Excluded regions file(s): (gziped) BED(s). Use ALL to exclude all.");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("i", "included", true, "Included regions file(s): (gziped) BED(s).");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("chrom", true, "Tested (compared) chromosomes (comma-separated). ");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("excludeNonPass", false, "Do not included filtered SNPs. Use with care.");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("dropAllFiltered", false,
			"If set, all matching positions where *all* variants are filtered will be ignored. Only useful w/o excludeNonPass.");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("dontCheckSort", false, "Do not check whether files are sorted. Use with care.");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("createWigs", false, "If set, WIG files will be created.");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("o", "out", true, "Output dir (will be created if not existing).");
		o.setRequired(true);
		options.addOption(o);

		o = new Option("s", "scoringSchema", true,
			"Scoring schema consisting of two comma-separated weights, a positive one for concordant calls and a negative one for discordant calls. (def: '1,-3').");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("t", "thresholds", true, "Tresholds for reliable/unreliable regions (def: '0.5,-0.001').");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("w", "win", true, "Window size in bp (def: 1000bp).");
		o.setRequired(false);
		options.addOption(o);

		o = new Option("swf", "smallWindowFactor", true,
			"The configured window size will be multiplied by this factor to calculate the smaller window size (def: 0.1).");
		o.setRequired(false);
		options.addOption(o);

		o = new Option(
			"swr",
			"smallWindowFactorRegions",
			true,
			"Optional BED file that contains regions that are considered rather unreliable (e.g., low-complexity regions). "
				+ "If configured, the signal interpolator will switch to a smaller window size (see smallWindowFactor) in such regions (def: N/A).");
		o.setRequired(false);
		options.addOption(o);

		options.addOption("v", "verbose", false, "be verbose.");

		line = parser.parse(options, args);
		if (line.hasOption("v"))
		    debug = true;
		else
		    debug = false;

		if (line.hasOption("excludeNonPass") && line.hasOption("dropAllFiltered"))
		    System.err.println("Warning: dropAllFilter cannot be applied when non-pass variants are excluded.");

		String scoringschema = line.hasOption("scoringSchema") ? line.getOptionValue("scoringSchema") : "1,-3";
		if (!scoringschema.contains(","))
		    throw new IOException("Invalid scoring schema. Cannot parse " + scoringschema);
		int scoreConcordantHit = Integer.parseInt(scoringschema.split(",")[0]);
		int scoreDiscordantHit = Integer.parseInt(scoringschema.split(",")[1]);
		if (!(scoreConcordantHit >= 0 && scoreDiscordantHit <= 0))
		    throw new IOException("Invalid scoring schema. weight for concordant calls must be positive, weight for discordant calls negative.");

		String thresholds = line.hasOption("thresholds") ? line.getOptionValue("thresholds") : "0.5,-0.001";
		if (!thresholds.contains(","))
		    throw new IOException("Invalid thresholds specification. Cannot parse " + thresholds);
		double reliableMin = Double.parseDouble(thresholds.split(",")[0]);
		double unreliableMax = Double.parseDouble(thresholds.split(",")[1]);
		if (reliableMin > 1 || reliableMin < -1 || unreliableMax > 1 || unreliableMax < -1 || reliableMin < unreliableMax)
		    throw new IOException("reliableMin & unreliableMax have to be from [-1, 1] and reliableMin>=unreliableMax");

		double smallWindowFactor = line.hasOption("smallWindowFactor") ? Double.parseDouble(line.getOptionValue("smallWindowFactor")) : 0.1d;

		calc(line.getOptionValues("d"), line.getOptionValues("e"), line.getOptionValues("i"), line.hasOption("chrom") ? line.getOptionValues("chrom")
			: CanonicalChromsomeComparator.canonicalChromosomes, scoreConcordantHit, scoreDiscordantHit, line.hasOption("dontCheckSort"),
			line.hasOption("excludeNonPass"), line.hasOption("dropAllFiltered"), line.hasOption("createWigs"),
			line.hasOption("win") ? Integer.parseInt(line.getOptionValue("win")) : 1000, reliableMin, unreliableMax, smallWindowFactor,
			line.getOptionValue("swr"), new File(line.getOptionValue("o")));

		System.err.println("Finished.");

		System.exit(0);
	    } else if (subcommand.equalsIgnoreCase("calcDensityHistogram")) {
		Option o = new Option("i", "in", true, "Input VCF file (usually the output of the calc() operation).");
		o.setRequired(true);
		options.addOption(o);

		o = new Option("w", "win", true, "Window size.");
		o.setRequired(true);
		options.addOption(o);

		o = new Option("c", "cutoff", true, "Density cutoff.");
		o.setRequired(true);
		options.addOption(o);

		o = new Option("o", "out", true, "Output BED file.");
		o.setRequired(true);
		options.addOption(o);

		options.addOption("v", "verbose", false, "be verbose.");

		line = parser.parse(options, args);
		if (line.hasOption("v"))
		    debug = true;
		else
		    debug = false;

		calcDensityHistogram(new File(line.getOptionValue("i")), Integer.parseInt(line.getOptionValue("w")),
			Integer.parseInt(line.getOptionValue("c")), new File(line.getOptionValue("o")));

		System.exit(0);

	    } else {
		usage(options, null, null);
	    }
	} catch (Exception e) {
	    e.printStackTrace();
	    usage(options, subcommand, e.toString());
	}
    }

}
