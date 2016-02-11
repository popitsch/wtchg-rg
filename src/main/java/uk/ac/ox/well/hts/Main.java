package uk.ac.ox.well.hts;

import java.security.NoSuchAlgorithmException;
import java.sql.SQLException;
import java.util.Arrays;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

import uk.ac.ox.well.hts.reliablegenome.BACKUPCalcreliabilitySignals;

/**
 * Main class. Implements CLI.
 * 
 * @author niko.popitsch@well.ox.ac.uk
 * 
 */
public class Main {

    private static final String VERSION = "v1.0";

    private static void usage(Options options) {

	System.out.println("Welcome to the WTCHG HTS package!");
	System.out.println("Version:\t" + VERSION);
	System.out.println();
	System.out.println("Usage:\t\tjava -jar x.jar <command> [options]:\t");
	System.out.println();
	System.out.println("Command:\t" + BACKUPCalcreliabilitySignals.CMD + "\t\t" + BACKUPCalcreliabilitySignals.CMD_INFO);
	System.exit(1);
    }

    /**
     * Main entry point
     * 
     * @param args
     * @throws CompressionException
     * @throws at.cibiv.ngs.tools.util.ParseException
     * @throws Throwable
     * @throws SQLException
     * @throws NoSuchAlgorithmException
     */
    public static void main(String[] args) throws Throwable {

	// create the command line parser
	CommandLineParser parser = new PosixParser();

	// create the Options
	Options options = new Options();
	options.addOption("h", "help", false, "Print this usage information.");

	try {
	    // parse the command line arguments
	    CommandLine line = parser.parse(options, args, true);

	    // validate that block-size has been set
	    if (line.getArgs().length == 0 || line.hasOption("h")) {
		usage(options);
	    }
	    String[] args2 = Arrays.copyOfRange(args, 1, args.length);

	    if (line.getArgs()[0].equalsIgnoreCase(BACKUPCalcreliabilitySignals.CMD)) {
		BACKUPCalcreliabilitySignals.main(args2);
	    } else
		usage(options);
	} catch (Exception exp) {
	    System.out.println("Unexpected exception:" + exp.getMessage());
	}
    }
}
