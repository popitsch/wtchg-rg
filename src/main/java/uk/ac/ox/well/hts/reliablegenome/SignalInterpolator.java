package uk.ac.ox.well.hts.reliablegenome;

import at.cibiv.ngs.tools.lds.GenomicITree;
import at.cibiv.ngs.tools.util.GenomicPosition;
import at.cibiv.ngs.tools.util.GenomicPosition.COORD_TYPE;
import at.cibiv.ngs.tools.wig.WigOutputStream;

/**
 * Interpolates score values for non-polymorphic genomic positions.
 * 
 * @author niko@well.ox.ac.uk
 * 
 */
public class SignalInterpolator {

    private long largeWin;
    private long smallWin;
    private long w;
    GenomicITree smallWindowSizeRegions;
    boolean inSmallWinRegion = false;

    GenomicPosition lastpos = null;
    double lastValue = 0;

    /**
     * Constructor.
     * 
     * @param largeWin
     * @param smallWinSizeFactor
     * @param smallWindowSizeRegions
     */
    public SignalInterpolator(long largeWin, double smallWinSizeFactor, GenomicITree smallWindowSizeRegions) {
	this.largeWin = largeWin;
	this.smallWin = (long) ((double) largeWin * smallWinSizeFactor);
	this.smallWindowSizeRegions = smallWindowSizeRegions;
	// we start outside of a small-win region
	this.inSmallWinRegion = false;
	this.w = largeWin;
    }

    /**
     * Add a data point to a WIG
     * 
     * @param pos
     * @param value
     * @param out
     */
    public void addDataPoint(GenomicPosition pos, double value, WigOutputStream out) {
	if (lastpos == null || !lastpos.onSameChromosome(pos))
	    lastpos = pos;
	long diff = pos.get1Position() - lastpos.get1Position();

	if (diff <= w) {
	    // special - windows overlap. Lin. interpos between points.
	    for (long x = lastpos.get1Position(); x < pos.get1Position(); x++) {
		GenomicPosition p = new GenomicPosition(pos.getChromosomeOriginal(), x, COORD_TYPE.ONEBASED);
		double perc = (double) (x - lastpos.get1Position()) / diff;
		double signal = lastValue + perc * (value - lastValue);
		if (smallWindowSizeRegions != null) {
		    if (inSmallWinRegion && !smallWindowSizeRegions.contains(p)) {
			// we just left a small-win-region.
			w = largeWin;
			lastpos = p;
			lastValue = signal;
			inSmallWinRegion = false;
			addDataPoint(pos, value, out);
			return;
		    } else if (!inSmallWinRegion && smallWindowSizeRegions.contains(p)) {
			// we just entered a small-win-region
			w = smallWin;
			lastpos = p;
			lastValue = signal;
			inSmallWinRegion = true;
			addDataPoint(pos, value, out);
			return;
		    }
		}
		out.push(p, signal);
	    }
	} else {
	    long w2 = w / 2;
	    long breakpoint1 = lastpos.get1Position() + w2;
	    long breakpoint2 = pos.get1Position() - w2;

	    for (long x = lastpos.get1Position(); x < pos.get1Position(); x++) {
		GenomicPosition p = new GenomicPosition(pos.getChromosomeOriginal(), x, COORD_TYPE.ONEBASED);

		double signal = 0;
		if (x <= breakpoint1) {

		    double perc = (1d - (double) (x - lastpos.get1Position()) / (double) w2);
		    signal = lastValue * perc;
		} else if (x >= breakpoint2) {
		    double perc = (1d - (double) (pos.get1Position() - x) / (double) w2);
		    signal = value * perc;
		} else
		    signal = 0;

		if (smallWindowSizeRegions != null) {
		    if (inSmallWinRegion && !smallWindowSizeRegions.contains(p)) {
			// we just left a small-win-region.
			w = largeWin;
			lastpos = p;
			lastValue = signal;
			inSmallWinRegion = false;
			addDataPoint(pos, value, out);
			return;
		    } else if (!inSmallWinRegion && smallWindowSizeRegions.contains(p)) {
			// we just entered a small-win-region
			w = smallWin;
			lastpos = p;
			lastValue = signal;
			inSmallWinRegion = true;
			addDataPoint(pos, value, out);
			return;
		    }
		}
		out.push(p, signal);
	    }

	}

	lastpos = pos;
	lastValue = value;
    }

}
