/*
 * A program that creates structures based on levels of probable pairs.
 * This program is based on the ProbablePair program written by
 * Jessica S. Reuter from the David Mathews lab
 */

#ifndef PROBABLE_PAIRS_RR_H
#define PROBABLE_PAIRS_RR_H

#include "../RNA_class/RNA.h"
#include "../src/ErrorChecker.h"
#include "../src/ParseCommandLine.h"

class ProbablePairRR {
 public:
	// Public constructor and methods.

	/*
	 * Name:        Constructor.
	 * Description: Initializes all private variables.
	 */
	ProbablePairRR();

	/*
	 * Name:        parse
	 * Description: Parses command line arguments to determine what options are required for a particular calculation.
	 * Arguments:
	 *     1.   The number of command line arguments.
	 *     2.   The command line arguments themselves.
	 * Returns:
	 *     True if parsing completed without errors, false if not.
	 */
	bool parse( int argc, char** argv );

	/*
	 * Name:        run
	 * Description: Run calculations.
	 */
	void run();

 private:
	// Private variables.

	// Description of the calculation type.
	string calcType;

	// Input and output file names.
	string input;            // The input plot data file.
	string ctFile;           // The output ct file.

	// Flag signifying if calculation handles RNA (true) or DNA (false).
	bool isRNA;

	// Boolean flag signifying if input is a sequence file (true) or not (false).
	bool isSequence;

	// Float designating what threshold should be used.
	float threshold;
};

#endif /* PROBABLE_PAIR_H */
