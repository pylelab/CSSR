#ifndef STRUCTURE_H
#define STRUCTURE_H

#include <string>
#include <stdlib.h>
#include <vector>
#include "defines.h"
#include "rna_library.h"
#include "CTCommentProvider.h"
#include "DotBracketFormat.h"
#include "../src/phmm/utils/xmath/log/xlog_math.h"

using namespace std;

#ifdef _WINDOWS_GUI
	#include "../Windows_interface_deprecated/platform.h"
#else
	#include "platform.h"
#endif //_WINDOWS

#ifdef EXTENDED_DOUBLE
	#include "extended_double.h" //inlcude code for extended double if needed
#endif//defined EXTENDED_DOUBLE

//This is a deprecated array size that must be removed:
#define maxforce 3000
#define DEFAULT_SHAPE_INTERCEPT -0.6  // default intercept used in SHAPE data pseudoenergy restraint generation (kcal/mol)
#define DEFAULT_SHAPE_SLOPE 1.8       // default slope used in SHAPE data pseudoenergy restraint generation (kcal/mol)

// Enable VERIFY_ARRAY_BOUNDS to force bounds-checking on arrays (for debugging and development only. Should be disabled for releases!)
//#define VERIFY_ARRAY_BOUNDS 

#ifdef VERIFY_ARRAY_BOUNDS 
	#define VERIFY_NUC_INDEX(INDEX) while (INDEX<1||INDEX>GetSequenceLength()) { cerr << "Nucleotide index out of bounds: " << INDEX << "(length " << GetSequenceLength() << ") at " << __FILE__ << ":" << __LINE__ << " base: " << GetBase(i) << endl; break; }
	#define VERIFY_NUC_INDEX_CT(INDEX,CT) while (INDEX<1||INDEX>CT->GetSequenceLength()) { cerr << "Nucleotide index out of bounds: " << INDEX << "(length " << CT->GetSequenceLength() << ") at " << __FILE__ << ":" << __LINE__ << endl; break; }
	#define VERIFY_STRUCTURE_INDEX(INDEX) while (INDEX<1||INDEX>GetNumberofStructures()) { cerr << "Structure index out of bounds: " << INDEX << "(length " << CT->GetNumberofStructures() << ") at " << __FILE__ << ":" << __LINE__ << endl; 
	#define VERIFY_STRUCTURE_INDEX_CT(INDEX,CT) while (INDEX<1||INDEX>CT->GetNumberofStructures()) { cerr << "Structure index out of bounds: " << INDEX << "(length " << CT->GetNumberofStructures() << ") at " << __FILE__ << ":" << __LINE__ << endl; break; }
#else
	// Define all the VERIFY_ARRAY_BOUNDS macros to do nothing.
	#define VERIFY_NUC_INDEX(INDEX)			
	#define VERIFY_STRUCTURE_INDEX(INDEX)
	#define VERIFY_NUC_INDEX_CT(INDEX,CT)
	#define VERIFY_STRUCTURE_INDEX_CT(INDEX,CT)
#endif

//Single structure is a wrapper for the information associated with just a single structure, i.e. pairs, energy, and labels.
struct singlestructure {

	//This function sizes the vectors to the approprate size:
	singlestructure(int sequencelength);

	//keep track of the basepairs
	vector<int> basepr;

	//keep track of the energy of the structure, if available
	int energy;

	//keep a string, from a ct file or sequence file with the sequence desription
	string ctlabel;

};

//! Constants used for the 'modifier' parameter for pseudo-energy-related functions
enum RestraintType { RESTRAINT_SHAPE, RESTRAINT_SHAPE_DIFF, RESTRAINT_SHAPE_AC,	RESTRAINT_SHAPE_GU, RESTRAINT_DMS, RESTRAINT_CMCT, RESTRAINT_DMSNT };

//! structure Class.
/*!
	The structure class provides a class for handling sequences and structures.
*/
//////////////////////////////////////////////////////////////////////
class structure //this structure contains all the info for a structure
{
	public:
		//******************************
		//Constructor:
		//******************************

		//!Constructor.
		//!	\param sructures is an int that specifies how many structures should be anticipated.  This sets up an initial memory allocation, but this can expand as needed.
		structure(int structures = maxstructures+1);

		//!Destructor.
		~structure();

		//*********************************
		//Get and receive sequence and structure information:
		//*********************************


		//! Get the energy for structure number structurenumber.

		//! This function requires that an energy calculation has been performed.  It does not do an eenrgy calculation.
		//! \param structurenumber is an int that gives the structure number being indexed, this is one indexed.
		//! \return an int that gives the energy in kcal/mol*conversionfactor, where conversionfactor is set in /src/defines.h.
		int GetEnergy(int structurenumber) const;
		
		//! Get the number of structures stored.

		//! \return An integer that is the number of structures encoded.
		int GetNumberofStructures() const;

		//! Get the pairing partner for i in structure structurenumber.

		//! \param i is the nucleotide index, which is one-indexed.
		//! \param structurenumber is the structure number, which is one-indexed.
		//! \return The pairing partner, as a nucleotide index.
		int GetPair(int i, int structurenumber=1) const;

		//! Get the base number at position i  (i.e. the index of the base in the alphabet).
		inline int GetBase(int i) const { return numseq[i]; }

		//! Get the label associated with the sequence.

		//! \return The string read from the sequence file.
		string GetSequenceLabel() const;

		//! Get the length of the sequence.

		//! \return An integer that is the sequence length.
		inline int GetSequenceLength() const {
			//Return the value of numofbases:
			return numofbases;
		}

		const char* GetSequence() const;
		

		//! Set the label for a structure, using a string.

		//! \param label is a string that will be strored.
		//! \param structurenumber is the index to which structure will hold the label.  This is one-indexed.
		void SetCtLabel(const string &label, const int structurenumber);
		

		//! Set the label for a structure,using a pointer to cstring.

		//! \param label is a pointer to char that provides the label.
		//! \param structurenumber is the index to which structure will hold the label.
		void SetCtLabel(const char *label, const int structurenumber);

		//! Set the energy for structure numer structurenumber.

		//! \param structurenumber is an int that gives the structure number being indexed, this is one indexed.
		//! \param energy is an int that sets the energy in kcal/mol*conversionfactor, where conversionfactor is set in /src/defines.h.
		void SetEnergy(int structurenumber, int energy);

		
		//Set the pairing partner for i in structure structurenumber.

		//! This function sets nucleotide i paired to j in structure number structurenumber.
		//! \param i is the nucleotide index of the first pairing partner, which is one-indexed.
		//! \param j is the nucleotide index of the second pairing partner.
		//! \param structurenumber is the structure number, which is one-indexed.
		void SetPair(int i, int j, int structurenumber=1);

		//! Set the label from a sequence, using a string.

		//! \param label is a string that will be stored.
		void SetSequenceLabel(const string& label);

		//! Set the sequence for this structure.
		//! This includes the following operations:
		//!   - Allocate space for the bases.
		//!   - Verify each base exists in the alphabet.
		//!   - Setup the arrays nucs, numseq, and hnumber.
		//!   - Check to see if any nucleotide needs to be single-stranded and call AddSingle for them.
		//! The sequence can contain whitespace, which is ignored.
		//! If an error occurs (such as the data-table has not yet been loaded) the return value 
		//! will be an RNA error code (i.e. it corresponds to one of those defined in RNA::GetErrorMessage)
		//! Additionally lastErrorDetails may be set (which can be queried with GetErrorDetails())
		//!
		//! \param sequence The nucleotide sequence for the structure.  This should contain only valid bases. No whitespace is allowed.
		int SetSequence(const string& sequence);


		//********************************
		// Get and Set constraint information
		//*********************************
		

		//! Allocate a bool array with information about whether any given pair is allowed.

		//! This function must be called after reading a sequence and before any template information is read or written.
		//! This mechanism is orthogonal to the functions AddForbiddenPair, GetForbiddenPair5, and GetForbiddenPair3.  It exists for the convenience of coding functiopns that need to forbid a large number of pairs.
		//! The memory use is cleaned up in the destructor.
		void allocatetem();


		//! Add a nucleotide to the list of those that must pair.

		//! \param i is an int that indicates the nucleotide position, one indexed. 
		void AddDouble(int i);

		//! Add a nucleotide to the list of Us in GU pairs.

		//! Note that there is no error checking.  This nucleotide must be a U.
		//! \param i is an int that indicates the nucleotide position, one indexed.
		void AddGUPair(int i);

		//! Add a nucleotide to the list of those accessible to traditional chemical modification.

		//! These nucleotides can only be unpaired, at the end of a helix, in a GU pair, or adjacent to a GU pair.
		//! \param i is an int that indicates the nucleotide position, one indexed.
		void AddModified(int i);

		//! Add a pair of nucleotides to the list of those that must form.

		//! Note that there is no error checking.  This should be an allowed base pair.
		//! \param i is an int that indicates the 5' nucleotide position, one indexed.
		//! \param j is an int that indicates the 3' nucleotide position, one indexed.
		void AddPair(int i, int j);

		//! Add a nucleotide to the list of those not able to pair.

		//! \param i is an int that indicates the nucleotide position, one indexed. 
		void AddSingle(int i);


		//!	Indicate if pairing distance is limited for structrure prediction methods.


		//! \return A bool that is true if the pairing distance has a limit.
		inline bool DistanceLimited() {

			return limitdistance;

		}

		//! Get a nucleotide that must be base paired.

		//! \return An int that gives the nucleotide position, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofDoubles()-1, inclusive.
		int GetDouble(int i);

		//! Get a nucleotide that must not be in a specific pair.

		//! \return An int that gives the nucleotide position for the 5' partner, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofForbiddenPairs()-1, inclusive.
		int GetForbiddenPair5(int i);

		//! Get a nucleotide that must not be in a specific pair.

		//! \return An int that gives the nucleotide position for the 3' partner, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofForbiddenPairs()-1, inclusive.
		int GetForbiddenPair3(int i);

		//! Get a nucleotide that must be a U in a GU pair.

		//! \return An int that gives the nucleotide position, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofGU()-1, inclusive.
		int GetGUpair(int i);

		//! Get a nucleotide that is accessible to chemical modification.

		//! These nucleotides can only be unpaired, at the end of a helix, in a GU pair, or adjacent to a GU pair.
		//! \return An int that gives the nucleotide position, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofMofidied()-1, inclusive.
		int GetModified(int i);

		//! Get the number of nucleotides forced to be double-stranded.

		//! \return An int that is the number of nucleotides constrained to be double stranded.
		int GetNumberofDoubles();

		//! Get the number of pairs that are forbidden.

		//! \return An int that is the number of forbidden pairs.
		int GetNumberofForbiddenPairs();
		
		//! Get the number of Us forced to be in GU pairs.

		//! \return An int that is the number of nucleotides constrained to be in GU pairs.
		int GetNumberofGU();
		
		//! Get the number of nucleotides that are accessible to chemical modification.

		//! \return An int that is the number of nucleotides constrained to be chemically modified.
		int GetNumberofModified();
		
		//! Get the number of nucleotides forced to be single-stranded.

		//! \return An int that is the number of nucleotides constrained to be single stranded.
		int GetNumberofSingles();
		
		//! Get the number of pairs that are constrained to occur.

		//! \return An int that is the number of forced pairs.
		int GetNumberofPairs();

		//! \return An int that is the number of folding domains.
		int GetNumberofDomains();
		
		//! Get a nucleotide that must be in a specific pair.

		//! \return An int that gives the nucleotide position for the 5' partner, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofPairs()-1, inclusive.
		int GetPair5(int i);

		//! Get a nucleotide that must be in a specific pair.

		//! \return An int that gives the nucleotide position for the 3' partner, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofPairs()-1, inclusive.
		int GetPair3(int i);

		//!	Provide the maximum distance between nucleotides that can pair.

		//! \return An int that is the maximum distance in sequence for nucleotides that can pair.
		inline int GetPairingDistanceLimit() {

			return maxdistance;

		}

		//! Get a nucleotide that must be single stranded.

		//! \return An int that gives the nucleotide position, one indexed.
		//! \param i is an int that gives the constraint number, which should be between 0 and GetNumberofSingles()-1, inclusive.
		int GetSingle(int i);
		


		//! Returns a nucleotide that defines the 5' end of a domain
		int GetDomain5(int i);
		
		
		//! Returns a nucleotide that defines the 3' end of a domain
		int GetDomain3(int i);
		

		//! Reset, i.e. remove, all constraints
		void RemoveConstraints();



		//*********************************
		//Functions for disk I/O
		//*********************************

		//! Open a CT File.
		//! This opens a ct file and stores all the information in this instance of structure.
		//! A non-zero return indicates and error in reading the file.
		//! \param ctfile is a pointer to a Null-terminated cstring that gives the filename, including any necessary path information.
		//! \return The return value is 0 on success and non-zero on error. The error code returned corresponds to those defined in RNA::GetErrorMessage.
		int openct(const char *ctfile);

		//! This opens a dot-bracket (DBN) file and stores all the information in this instance of structure.
		//! A non-zero return indicates and error in reading the file.
	    //! Note: This can parse dot-bracket files with or without pseudoknots and automatically distinguishes 
	    //!       between  the following DBN formats (see DotBracketFormat for details): 
		//!       DBN_FMT_SINGLE_TITLE,  DBN_FMT_SIDE_TITLES, or DBN_FMT_MULTI_TITLE. 
		//!       This function cannot parse the format DBN_FMT_MULTI_TITLE_AND_SEQ because a structure class can only contain a single sequence.
		//! \param bracketFile is a pointer to a Null-terminated cstring that gives the filename, including any necessary path information.
		//! \return The return value is 0 on success and non-zero on error. The error code returned corresponds to those defined in RNA::GetErrorMessage.
		int opendbn(const char *bracketFile);

		//! Open a sequence file.
		//! This function works on .seq and FASTA files as well as plain-text sequences.
		//! Importantly, this function differs from the legacy openseq in its return codes. openseq returns 1 on 
		//! success and 0 on error, while this function returns 0 on success and an RNA error code on error.
		//! 
		//! \param seqfile is a const char pointer to a cstring that gives the filename, including any path information.
		//! \return An int that indicates an error state: 0 on no error or a number indicating a more detailed error code. See RNA::GetErrorMessage(int)
		int openseqx (const char *seqfile);

		//*******************************
		//Functions that act on whole structures
		//*******************************

		//! Add another empty structure to the list of singlestructures.
		// DHM: Remember if this is structure 1 to set the label from some sequence label.!
		void AddStructure();

		//! Remove all pairs from a structure, i.e. make it a clean slate for a new set of pairs
		
		//! \param struturenumber is an index to the structure to be cleaned.
		void CleanStructure(int structurenumber);

		

		//! Remove the last structure.
		void RemoveLastStructure();
		

		//********************************
		//Functions for accessing the datatable pointer or opening the datatables
		//********************************

		//! Get the datatable pointer.

		//! This function returns the pointer to the datatable, data.
		//! \return The pointer to the underlying datatable.
		datatable *GetThermodynamicDataTable();

		//! Set the datatable pointer.

		//! This function sets the pointer to the datatable, data.
		//! \param DataTablePointer is the pointer to a datatable.
		void SetThermodynamicDataTable(datatable *DataTablePointer);

		//! Returns true if the data property has been set to a valid datatable and its alphabet has been successfully read.
		bool IsAlphabetLoaded();

		//********************************
		//Additional functions
		//********************************

		//! Sort structures by energy.

		//! This function sorts structures in energy from lowest to highest.
		//! It is important that the structure energies be present by structure prediction or by efn2.
		void sort();



		//! This function reads a SHAPE reactivity datafile and parse the data into single-stranded amd chemical modification constraints.
		//! This function is largely depracated by the pseudo-free energy approach.  It is still available for experimentation.
		//! \return 0 on success or an error code compatible with RNA::GetErrorMessage (for example 1=file not found, 2=could not open file)
		int ReadSHAPE(const char *filename, float SingleStrandThreshold, float ModificationThreshold);//Read SHAPE reactivity data from a file

		//void ReadSHAPE(const char *filename, bool calculate=true);//Read SHAPE reactivity data from a file
		//! This function reads a SHAPE reactivity datafile and saves the data for a linear penalty.
		//! \param calculatePseudoEnergies (default true) indicate whether these data are being read for folding.  
		//!        (false means the raw values need to be stored.)
		//! \param modifier One of the RestraintType enum values to indicate which type of restraint is to be calculated (e.g. SHAPE, DMS, CMCT, etc)
		//! \return 0 on success or an error code compatible with RNA::GetErrorMessage (for example 1=file not found, 2=could not open file)
		int ReadSHAPE(const char *filename, RestraintType modifier=RESTRAINT_SHAPE, bool calculatePseudoEnergies=true);//Read SHAPE reactivity data from a file

		//! This function reads an experimental pair bonus file, similar to SHAPE, but just straightforward
		//! application as kcal bonuses.  As with SHAPE, bonus is applied at 0x, 1x, and 2x for
		//!  single stranded, edge base pairs, and internal base pairs.
		//! \return 0 on success or an error code compatible with RNA::GetErrorMessage (for example 1=file not found, 2=could not open file)		
		int ReadExperimentalPairBonus(const char *filename, double const experimentalOffset = 0.0, double const experimentalScaling = 1.0 );



		//! Deletes the SHAPE, SHAPEss, and SHAPEss_region arrays if the have been allocated.
		//! Also sets the `shaped` field to false and sets the SHAPE, SHAPEss, and SHAPEss_region 
		//! pointers to NULL. This can be called any time to remove SHAPE data.
		void DeleteSHAPE();

		//! Allocates the SHAPE and SHAPEss arrays if they haven't already been created.
		//! Also initializes SHAPEss_region, a triangular 2-D array that stores ss SHAPE energies for loops.
		void AllocateSHAPE();

		//! Return an exact copy of SHAPE data (e.g. for backup-and-restore purposes). 
		//! The return value is a pointer to a dynamically allocated array of doubles.
		//! \param includeSHAPEss Indicates whether SHAPEss data should be included. 
		//!        If true, the array will contain both SHAPE and SHAPEss values and will
		//!          have a size of 4*numbases+2 (with SHAPE data from arr[0] to arr[2N] and 
		//!          SHAPEss from arr[2N+1] to arr[4N+1]).
		//! 	   If false, the array will only contain SHAPE data and will have a size of 
		//!          2*numofbases+1.
		double* CopySHAPE(const bool includeSHAPEss);

		//! Load data from an array of doubles into the SHAPE array (e.g. for backup-and-restore purposes).
		//! \param shapeArray The array of values to load into the SHAPE (and optionally SHAPEss) arrays.
		//!    The data is copied -- this function does NOT take ownership of the passed-in array.
		//! \param includeSHAPEss Indicates whether or not SHAPEss data is included in the passed-in array. 
		//!          If true, the array must contain both SHAPE and SHAPEss values and must
		//!            have a size of at least 4*numbases+2 (with SHAPE data from arr[0] to arr[2N] 
		//!            and SHAPEss from arr[2N+1] to arr[4N+1])
		//!          If false, the array need only contain SHAPE data and must have a size of 
		//!            at least 2*numofbases+1.
		//! If the passed-in array is NULL, this structure's SHAPE data will be deleted.
		void LoadSHAPE(const double* shapeArray, const bool includeSHAPEss);

		double **constant;//constant is used to hold an array of equilibrium constants.  In partition function calculations, 
					//the equilibrium constant is multiplied by constant[j][i] when the i-j pair is formed. 
					//NOTE: The use of constant is NOT orthogonal to using chemical modification data.  They cannot
					//both be used at once.
		bool SHAPEFileRead;


		//!	Returns extended details about the last error. (e.g. error messages produced during file read operations that are otherwise lost.)
		const string& GetErrorDetails();


		//Is a nucleotide index a certain type of nucleotide?:
		//return true if i is type c or false otherwise
		// inline bool IsNuc(int i, char c) {
			
		// 	if (std::find(data->alphabet[numseq[i]].begin(), data->alphabet[numseq[i]].end(), c) != data->alphabet[numseq[i]].end())
		// 		return true;
		// 	else return false;

		// }
		inline bool IsNuc(int i, char c) {
			VERIFY_NUC_INDEX(i);
			return (std::find(data->alphabet[numseq[i]].begin(), data->alphabet[numseq[i]].end(), c) != data->alphabet[numseq[i]].end());
		}
		inline bool IsNuc(int i, char c, datatable *data) {
			VERIFY_NUC_INDEX(i);
			return (std::find(data->alphabet[numseq[i]].begin(), data->alphabet[numseq[i]].end(), c) != data->alphabet[numseq[i]].end());
		}

		//! Sets the default behavior regarding non-critical warnings. 
		//! 1 (ON)  -- Write warnings to stdout.
		//! 2 (ERR) -- Write warnings to STDERR for better error detection in scripts etc..
		//! 0 (OFF) -- Suppress warnings (not recommended unless the user is aware of 
		//!            problematic data and does not need the clutter of warning output.
		//! The default value of this variable is taken from the environment variable
		//! RNA_WARNINGS which should be set to "OFF", "ON", or "ERR" or 
		//! the numeric equivalents "0", "1", or "2". 
		static int ShowWarnings;

		//! SumShapeRepeats affects how multiple datapoints for the same nucleboase are handled 
		//! in SHAPE data. True (the default) causes ReadSHAPE to add multiple datapoints for use
		//! in the "resample with replacement" technique, but setting the environment variable 
		//! AVG_SHAPE_REPEATS to 1 causes multiple values to be averaged instead.
		static bool SumShapeRepeats; // default true.

		//static const double kT;


#ifdef SWIG
		// this hides the following definitions from SWIG when generating proxy-code for Java etc.
		private:
#endif // SWIG

		//! Returns a reference to an output stream appropriate for warnings (according to 
		//! the structure::ShowWarnings setting.)
		ostream& cwarn(); 

		string sequencelabel;//a label that was read from disk along with a sequence

		
		short int *numseq,*hnumber;
		
		int inter[3],allocatedstructures;
		char *nucs;
		bool intermolecular,allocated,templated,stacking;
		bool **tem;//tem stores template information as to whether a pair is allowed

		void allocate(int size = maxbases);
		void allocatestructure(int structures);
		
		
		
		 
		short int min_gu, min_g_or_u;//NMR-derived constraint variables
		short int neighbors[maxforce][maxneighborlength],nneighbors;//also NMR-derived index this from zero in both dimensions
		//regional NMR constraints:
		short int nregion,rmin_gu[maxregions],rmin_g_or_u[maxregions];
		short int rneighbors[maxregions][maxforce][maxneighborlength],rnneighbors[maxregions],start[maxregions],stop[maxregions];
		//microarray type constraints:
		short int nmicroarray,microstart[maxregions],microstop[maxregions],microunpair[maxregions];
		bool *fcedbl;//pointer to a 2-D array used in Dynalign to track nucleotides that must be double-stranded
		
		
		double *SHAPE;//double array to contain SHAPE data -- values less than -500 are ignored
		double **EX;// double array that contains experimental bonuses/penalties
		bool shaped;//keeps track of whether SHAPE data was loaded
		bool experimentalPairBonusExists;//keeps track of whether experimental bonus data was loaded
		bool ssoffset;//keeps track of whether a single stranded offset was read from disk
		double SHAPEslope,SHAPEintercept;//values of slope and intercept for SHAPE data modification of pairing stability
		//SINGLE STRANDED SHAPE ENERGY VARIABLES AND FUNCTIONS
		double *SHAPEss; //short int array that contains SHAPE data for single-stranded segments
		double SHAPEslope_ss, SHAPEintercept_ss; //values of the slope and intercept for SHAPE data modifying single stranded loop stability
		short int **SHAPEss_region;  //2-d short int array containing energy values for hairpin loop combinations
		int SHAPEss_calc(int index_i, int index_j);  //Returns pseudoenergy term for a hairpin loop using single stranded SHAPE data
		short int SHAPEss_give_value(int index);  //Returns the single stranded SHAPE pseudo energy for a given nucleotide
		double CalculatePseudoEnergy(const double data, const RestraintType modifier, const double, const double, const int ntcode, const bool);
		double Gammadist(const double data, const double shape, const double loc, const double scale);
		double Potential(const double data, const std::vector< std::vector<double> > &params, const double kT, const int ntcode = 1);
        
		void ReadProbabilisticPotentialParams();//Read chemical modifier distributions from file
		
        //Parameters for distributions
		std::vector< std::vector<double> > SHAPE_params;
		std::vector< std::vector<double> > DMS_params;
		std::vector< std::vector<double> > DMS_paramsnt;
        std::vector< std::vector<double> > CMCT_params;

		bool distsread;//keep track if the distribution files have been read from disk.

		// energyLabelWriter is a pointer to a function that creates an energy label to be inserted into structure 
		// titles when writing CT and dot-bracket files
		const char* const (*energyLabelWriter)(const int structurenumber);

	private:
		
		//!	Set extended details about the last error. (e.g. error messages produced during file read operations that are otherwise lost.)
		void SetErrorDetails(const string &details); // , bool outputToCErr = true

		//! Fills SHAPEss_region, a 2-d array with pseudo energy terms for loops from i-j using 
		//! single-stranded (ss) SHAPE parameters or offsets.
		//! The SHAPEss and SHAPEss_region arrays must already be allocated (by calling AllocateSHAPE)
		//! and SHAPEss must already contain ss SHAPE pseudo-energies.
		//! This is called from e.g. ReadSHAPE and ReadOffset.
		void FillSHAPEssRegions();

		int numofbases;//number of nucleotides in sequence
		bool limitdistance;//toggle to indicate that there is a limit on the maximum distance between nucs in base pairs
		int maxdistance;//maximum distance between nucs in base pairs
		
		vector<singlestructure> arrayofstructures;//This holds an array of structures, i.e. base pairing information and comments
			
		//variables for holding folding constraints:
		vector<int> doublestranded; //nucleotides that must be double stranded
		vector<int> singlestranded; //nucleotides that must be single stranded
		vector<int> GUpair; //Us in GU pairs
		vector<int>	modified; //nucleotides accessible to tradictional chemical modification agents
		vector<int> pair5; //5' partner in forced pair
		vector<int> pair3; //3' partner in forced pair
		vector<int> forbid5; //5' partner in a forbidden pair
		vector<int> forbid3; //3' partner in a forbidden pair
		vector<int> domains5; //domain definitions
		vector<int> domains3; //domain definitions

		string lastErrorDetails;
		datatable *data;//store the thermodynamic data, and access the alphabet of nucleotides
};


//char *tobase (int i);//convert a numeric value for a base to the familiar
								//character


//void tonum(char *base,structure *ct,int count); //converts base to a numeric

#ifndef SWIG // The following should be excluded from SWIG code generation


integersize ergcoaxflushbases(int i, int j, int ip, int jp, datatable *data);
//this function calculates flush coaxial stacking
//it requires sequence in i,j,ip, and jp
integersize ergcoaxinterbases1(int i, int j, int ip, int jp, int k, int l, datatable *data); 
//this funtion calculates an intervening mismatch coaxial stack
integersize ergcoaxinterbases2(int i, int j, int ip, int jp, int k, int l, datatable *data); 
//this funtion calculates an intervening mismatch coaxial stack
integersize ergcoaxflushbases(int i, int j, int ip, int jp, structure *ct, datatable *data);
integersize ergcoaxinterbases1(int i, int j, int ip, int jp, structure *ct, datatable *data);
integersize ergcoaxinterbases2(int i, int j, int ip, int jp, structure *ct, datatable *data);


integersize erg1(int i,int j,int ip,int jp,structure *ct,datatable *data);
		//calculates energy of stacked base pairs
integersize erg2(int i,int j,int ip,int jp,structure *ct,datatable *data,char a,
	char b);
		//calculates energy of a bulge/internal loop
integersize erg3(int i,int j,structure *ct,datatable *data,char dbl);
		//calculates energy of a hairpin loop
integersize erg4(int i,int j,int ip,int jp,structure *ct,datatable *data,
	bool lfce);
		//calculates energy of a dangling base
//erg4 without parameter usage counting
integersize erg4_nc(int i,int j,int ip,int jp,structure *ct, datatable *data, bool lfce);

inline integersize SHAPEend(int i, structure *ct);//calculate the SHAPE pseudo energy for a single
									//paired nucleotide

//this function calculates whether a terminal pair i,j requires the end penalty
inline integersize penalty(int i,int j,structure* ct, datatable *data) {
	integersize energy;
	VERIFY_NUC_INDEX_CT(i,ct);
	VERIFY_NUC_INDEX_CT(j,ct);
	if (ct->IsNuc(i,'U')||ct->IsNuc(j,'U'))
   	energy= data->auend;
	else energy= 0;//no end penalty
	
	return (energy/*+SHAPEend(i,ct)+SHAPEend(j,ct)*/);

}

inline integersize penalty_nc(int i,int j,structure* ct, datatable *data) {
	integersize energy;

#ifdef COUNTING
	int tmp_count = data->auend.get;
#endif

	energy = penalty(i,j,ct,data);

#ifdef COUNTING	
	data->auend.get = tmp_count;
#endif

	return energy;
}


#endif // Place functions that should be available to client languages (e.g. Java, python etc) BELOW this line.


#endif //STRUCTURE_H
