//Class RNA -- Wrapper for RNAstructure for use in object-oriented applications

//Use the precompiler to make sure the class definition is not included more than once.
#if !defined(RNA_H)
#define RNA_H



//Include all required source here to ease use by end user.  
#include <string>
#include <cstring>
#include <utility>
#include "../src/defines.h"
#include "../src/rna_library.h"
#include "../src/pfunction.h"
#include "thermodynamics.h"
#include "../src/TProgressDialog.h"
#include "../src/phmm/utils/xmath/log/xlog_math.h"


// Allow some types to be used without being fully defined by header inclusion.
// This is useful for minimal builds such as partition-rosetta.
struct coordinates;

//! RNA Class.
/*!
	The RNA class provides an entry point for all the single sequence operations of RNAstructure.
*/

//Note the stylized comments provide facility for automatic documentation via doxygen.

//! The following constants can be supplied to the RNA class constructor RNAInputType parameter 
//! to indicate the type of file represented by the filepathOrSequence parameter.
//!  0 => A literal sequence (NOT a file)
//!  1 => Structure File (CT, DBN/DOT/BRACKET)
//!  2 => Sequence File (SEQ, FASTA, or text)
//!  3 => Partition-Function-Save-File (PFS)
//!  4 => Folding-Save-File (FSV or SAV)
//!  5 => Bracket File (DBN/DOT/BRACKET) (1 will also work, but is slighly less efficient. 5 is more strict)
typedef int RNAInputType; // Just to help new programmers find the right constants (defined below)

//! A string containing the sequence -- NOT a file name.
#define SEQUENCE_STRING 0
//! CT Structure file
#define	FILE_CT 1
//! Sequence file (FASTA, SEQ or plain-text sequence)
#define	FILE_SEQ 2
//! Partition function save file (*.pfs)
#define	FILE_PFS 3
//! Folding Save file (*.sav)
#define	FILE_SAV 4
//! Dot-Bracket Notation Structure file (*.dbn, *.dot, *.bracket)
#define	FILE_DBN 5

class RNA: public Thermodynamics {
	public:
		//******************************
		//Constructors:
		//******************************
		

		//!Constructor - user provides a sequence as a c string.

		//!	Input sequence should contain A,C,G,T,U,a,c,g,t,u,x,X.
		//!	The sequence is case-sensitive. Lower case letters are marked as unpaired.
		//!	T=t=u=U.  If IsRNA is true, the backbone is RNA, so U is assumed.  If IsRNA is false, the backbone is DNA, so T is assumed.
		//!	x=X= nucleotide that neither stacks nor pairs.
		//!	Unknown nucs result in an error.
		//! Note that sequences will subsequently be indexed starting at 1 (like a biologist), so that the 0th position in the sequence array will be nucleotide 1.
		//!	\param sequence is a NULL terminated c string.
		//!	\param IsRNA is a bool that indicates whether this sequence is RNA or DNA.  true=RNA.  false=DNA.  Default is true.
		RNA(const char sequence[], const bool IsRNA=true); 

		//!Constructor - user provides a filename for existing file as a c string.
		//!	The existing file, specified by filename, can either be a ct file, a sequence, or an RNAstructure save file. 
		//!	This constructor generates internal error codes that can be accessed by GetErrorCode() after the constructor is called.  0 = no error.
		//! The errorcode can be resolved to a c string using GetErrorMessage.		
		//!	Note that the contructor needs to be explicitly told, via IsRNA, what the backbone is because files do not store this information.
		//! Note also that save files explicitly store the thermodynamic parameters, therefore changing the backbone type as compared to the original calculation will not change structure predictions.
		//! \param filepathOrSequence is null terminated c string that can be either a path to a CT, SEQ, FASTA, PFS, SAV, or DBN file or a full RNA sequence (whitepsace allowed)
		//! \param fileType is an RNAInputType (e.g FILE_CT, FILE_SEQ etc) that indicates what type of file filepathOrSequence refers to. 
		//!             If filepathOrSequence is an in-place sequence, fileType should be SEQUENCE_STRING.
		//!	\param alphabetName is a c-string that specifies the alphabet name (e.g. "rna",  "dna", etc)
		//!	\param allowUnknownBases boolean that determines the handling of unknown bases (i.e. those not defined in the alphabet)
		//!	       If allowUnknownBases is true, they will be interpreted as bases that neither pair nor stack. 
		//!	       Otherwise they will result in an error.
		//!	\param skipThermoTables If true, only the alphabet specification will be loaded, not the full thermodynamic energy tables.
		//!             The default is false -- i.e. load all tables. This should only be set to true in programs that never make use of energies -- i.e. ct2dot etc.
		RNA(const char filepathOrSequence[], const RNAInputType fileType, const char* const alphabetName, const bool allowUnknownBases=false, const bool skipThermoTables=false);

		//!Constructor - user provides a filename for existing file as a c string.
		//! Construct an RNA from either a sequence or file. Copy the Thermodynamics info and datatable from another Thermodynamics instance.
		//! \param filepathOrSequence is null terminated c string that can be either a path to a CT, SEQ, FASTA, PFS, SAV, or DBN file or a full RNA sequence (whitepsace allowed)
		//! \param fileType is an RNAInputType (e.g FILE_CT, FILE_SEQ etc) that indicates what type of file filepathOrSequence refers to. 
		//!             If filepathOrSequence is an in-place sequence, fileType should be SEQUENCE_STRING.
		//! \param copyThermo A pointer to an existing Thermodynamics from which to copy the datatables.
		RNA(const char filepathOrSequence[], const RNAInputType fileType, const Thermodynamics *copyThermo);

		//!Constructor - user provides a filename for existing file as a c string.

		//!	The existing file, specified by filename, can either be a ct file, a sequence, or an RNAstructure save file. 
		//!	Therefore, the user provides a flag for the file: 
		//!		1 => .ct file,  2 => .seq file,  3 => partition function save (.pfs) file,  4 => folding save file (.sav), 5 => DotBracket (DBN) file
		//!	This constructor generates internal error codes that can be accessed by GetErrorCode() after the constructor is called.  0 = no error.
		//! The errorcode can be resolved to a c string using GetErrorMessage.		
		//!	Note that the contructor needs to be explicitly told, via IsRNA, what the backbone is because files do not store this information.
		//! Note also that save files explicitly store the thermodynamic parameters, therefore changing the backbone type as compaared to the original calculation will not change structure predictions.
		//! \param filename is null terminated c string.
		//! \param fileType is an integer that indicates the file type.
		//!	\param IsRNA is a bool that indicates whether this sequence is RNA or DNA.  true=RNA.  false=DNA.  Default is true.
		RNA(const char filename[], const RNAInputType fileType, const bool IsRNA=true);

		//! Default Constructor - user provides nothing.
		//! This basic constructor is provided for bimolecular folding and should not generally need to be accessed by end users of the RNA class.
		//!	\param IsRNA is a bool that indicates whether this sequence is RNA or DNA.  true=RNA.  false=DNA.  Default is true.
		RNA(const bool IsRNA=true);


		//******************************************************
		//Functions to return error information
		//******************************************************

		//!Return an error code, where a return of zero is no error.

		//!	This function returns and error flag that is generated during construction by RNA(const char &filename, const int type, const bool IsRNA=true) or from CalculateFreeEnergy().
		//!		An error of zero is always no error.  Other codes are errors and a c-string can be fetched for the error with GetErrorMessage().
		//!\return An integer that provides the error code.
		int GetErrorCode() const;

		//!	Return error messages based on code from GetErrorCode and other error codes.		

		//!		0 = no error
		//!		1 = input file not found
		//!		2 = error opening file
		//!		3 = structure number out of range
		//!		4 = nucleotide number out of range
		//!		5 = error reading thermodynamic parameters
		//!		6 = pseudoknot formation
		//!		7 = non-canonical pair
		//!		8 = too many restraints specified
		//!		9 = same nucleotide in conflicting restraint
		//!		10 = no structures to write
		//!		11 = nucleotide not a U (caused by ForceFMNCleavage()
		//!		12 = distance too short
		//!		13 = error reading constraint file
		//!		14 = traceback error
		//!		15 = no partition function data present
		//!		16 = incorrect save file version used
		//!		17 = cannot be performed without having read a save file (.sav)
		//!		18 = threshold is too low to be valid
		//!		19 = drawing coordinates have not been determined
		//!		20 = no sequence has been read
		//!		21 = over 1 probability error on stochastic traceback
		//!		22 = programming error, unrecognized input to constructor
		//!		23 = no structures present
		//!		24 = too few iterations
		//!		25 = index (for drawing) is not a multiple of 10
		//!\param error is the integer error code provided by GetErrorCode() or from other functions that return integer error codes.
		//!\return A pointer to a c string that provides an error message.
		static const char* GetErrorMessage(const int error);

		//!	Returns extended details about the last error. (e.g. error messages produced during file read operations that are otherwise lost.)
		const string GetErrorDetails() const;

		//!	Set extended details about the last error. (e.g. error messages produced during file read operations that are otherwise lost.)
		void SetErrorDetails(const string& details);

		//!	Set the name of the sequence. (This must be called before structure prediction for it to affect the resulting structures.)
		void SetSequenceLabel(const string& label);


		
		//******************************************************************
		//Functions that calculate folding energies for existing structures:
		//******************************************************************

		//!Return the predicted Gibb's free energy change for structure # structurenumber, defaulted to 1.

		//!	Free energies are in kcal/mol.
		//!	The first time this is called, if no other free energy calculation has been performed and the folding temperature has not been specifed,
		//!		thermodynamic parameter files (.dat) files will be read from disk.
		//!	The parameter files should be located in the directory specified by environment 
		//!		variable $DATAPATH, or the pwd. 
		//!	In case of error, the function returns a free energy change of zero.
		//!		Note!: That a free energy change of zero is also a valid folding free energy change.
		//!	Errors will also generate an internal error code, accessible with GetErrorCode().
		//! GetErrorCode() will return 0 when there is no error and other codes can be parsed by GetErrorMessage() or GetErrorMessageString(). 
		//! \param structurenumber is an integer that refers to the index of the structure for which to calculate the folding free energy change.  This defaults to 1.
		//! \param UseSimpleMBLoopRules is a bool that indicates what energy rules to use.  The default, false, uses the complete nearest neighbor model for multibranch loops.  When true is passed, the energy model is instead a simplified model that is the one used by the dynamic programming algorithms.
		//!	\return A double which is the folding free energy change in kcal/mol.
		double CalculateFreeEnergy(const int structurenumber = 1, const bool UseSimpleMBLoopRules = false);

		double ExteriorLoopCorrection(const int structurenumber, const bool UseSimpleMBLoopRules, int min_index, int max_index);

		//!Calculate the folding free energy change for all structures and write the details of the calculation to a file.

		//!	Free energies are in kcal/mol.
		//!	The first time this is called, if no other free energy calculation has been performed and the folding temperature has not been specifed,
		//!		thermodynamic parameter files (.dat) files will be read from disk.
		//!	The parameter files should be located in the directory specified by environment 
		//!		variable $DATAPATH, or the pwd. 
		//!	In case of error, the function returns a non-zero.
		//! \param filename is a c-string that provides the name of the output file to be written.
		//!        If filename==NULL, free energies will be calculated, but no file will be written.
		//! \param UseSimpleMBLoopRules is a bool that indicates what energy rules to use.  The default, false, uses the complete nearest neighbor model for multibranch loops.  When true is passed, the energy model is instead a simplified model that is the one used by the dynamic programming algorithms.
		//!	\return An int that indicates whether an error occurred (0 = no error; 5 = error reading parameter files).
		int WriteThermodynamicDetails(const char filename[], const bool UseSimpleMBLoopRules = false);

		//***********************************************
		//Functions that predict RNA secondary structures
		//***********************************************

		//! Predict the lowest free energy secondary structure and generate suboptimal structures using a heuristic.

		//! This function predicts the lowest free energy structure and suboptimal structures.
		//! If the temperature has not been specified using SetTemperature and no free energies have been calculated, the
		//!		thermodynamic parameters have not been read and therefore they will be read by this function call.  The 
		//!		parameter files should be located in the directory specified by the environment variable $DATAPATH of the pwd.
		//!	In case of error, the function returns a non-zero that can be parsed by GetErrorMessage() or GetErrorMessageString().
		//!	\param percent is the maximum % difference in free energy in suboptimal structures from the lowest free energy structure.  The default is 20.
		//!	\param maximumstructures is the maximum number of suboptimal structures to generate.  The default is 20.
		//!	\param window is a parameter that specifies how different the suboptimal structures should be from each other (0=no restriction and larger integers require structures to be more different).  The defaults is 5, but this should be customized based on sequence length.
		//!	\param savefile is c string containing a file path and name for a savefile (.sav)that can be used to generate energy dot plots and to refold the secondary structure using different suboptimal structure parameters.  The default is "", which results in no save file written.
		//!	\param maxinternalloopsize is the maximum number of unpaired nucleotides in bulge and internal loops.  This is used to accelerate the prediction speed.  The default is 30.
		//! \param mfeonly is a bool that indicates whether only the minimum free energy structure will be generated.  This saves half the calculation time, but no save file can be generated.  Default is false.
		//! \return An int that indicates an error code (0 = no error, 5 = error reading thermodynamic parameter files, 14 = traceback error).
	int FoldSingleStrand(const float percent=20, const int maximumstructures=20, const int window=5, const char savefile[]="", const int maxinternalloopsize = 30, bool mfeonly=false, bool simple_iloops = true, bool disablecoax=false);


		//! Predict the lowest free energy secondary structure and generate all suboptimal structures.

		//! This function predicts the lowest free energy structure and suboptimal structures.  
		//! If the temperature has not been specified using SetTemperature and no free energies have been calculated, the
		//!		thermodynamic parameters have not been read and therefore they will be read by this function call.  The 
		//!		parameter files should be located in the directory specified by the environment variable $DATAPATH of the pwd.
		//!	In case of error, the function returns a non-zero that can be parsed by GetErrorMessage() or GetErrorMessageString().
		//! Two controls are available for limiting the number of structures, the maximum % difference in energy (percent) and the maximum absolute change in energy (deltaG).  The smaller of the two will be used as the limit.
		//!	\param percent is the maximum % difference in free energy in suboptimal structures from the lowest free energy structure.  The default is 5.
		//!	\param deltaG is the maximum difference in free energy change above the lowest free energy structure (in kcal/mol).  The defaults is 0.6 kcal/mol.
		//! \return An int that indicates an error code (0 = no error, non-zero = error).
		int GenerateAllSuboptimalStructures(const float percent=5, const double deltaG=0.6);


		//! Predict the structure with maximum expected accuracy and suboptimal structures.

		//! This function predicts structures composed of probable base pairs and single-srtranded nucleotide, weighted by gamma.
		//! The score for a structure is = gamma * 2 * (sum of pairing probabilities for pairs) + (sum of unpairing probabilities for single stranded nucleotides).
		//! This function requires partition function data from either a previous partition function calculations or
		//!		from having read a partition function save file during construction of the class.
		//!	In case of error, the function returns a non-zero that can be parsed by GetErrorMessage() or GetErrorMessageString().
		//!	\param maxPercent is the maximum percent difference is score in generating suboptimal structures.  The default is 20.
		//!	\param maxStructures is the maximum number of suboptimal structures to generate.  The default is 20.
		//! \param window is the window parameter, where a higher value generates suboptimal structures that are more different from each other.  The default is 1.
		//! \param gamma is the weight given to base pairs
		//! \return An int that indicates an error code (0 = no error, non-zero = error). 
		int MaximizeExpectedAccuracy(const double maxPercent=20, const int maxStructures=20, const int window=1, const double gamma=1.0);

		//!Predict the partition function for a sequence.

		//!This function must be called to predict base pair probabilities, perform stochastic traceback, or for maximizing expected accuracy.
		//! If the temperature has not been specified using SetTemperature and no free energies have been calculated, the
		//!		thermodynamic parameters have not been read and therefore they will be read by this function call.  The 
		//!		parameter files should be located in the directory specified by the environment variable $DATAPATH of the pwd.
		//!	In case of error, the function returns a non-zero that can be parsed by GetErrorMessage() or GetErrorMessageString().
		//! Note that the parameter temperature is used when calculating equilibrium constants, but does not change the temperature
		//!		at which the free energies are determined.  SetTemperature, from the underlying base class Thermodynamics, should be used to
		//!		change the temperature for most calculations.  This parameter should generally not be used.  The default is -10.0 and values below zero cause
		//!		this parameter to be ignored (the correct default behavior).  Note also that if SetTemperature is not used, the temperature defaults to
		//!		310.15 K (37 deg. C), which is the desired behavior for most purposes.
		//! \param savefile is a c string that contains the path and filename for creating a save file.  This defaults to "", which indicates no file is to be written.
		//! \param temperature is a double that indicates a pseudo-temperature for calculating equilibrium constants from  free energies at fixed temperature previously specified.
		//! \param restoreSHAPE Whether SHAPE data should be restored after this call. The default is true. 
		//!                     PartitionFunction converts SHAPE pseudo-energies into equillibrium constants, overwriting the 
		//!                     GetStructure()->SHAPE array. If restoreSHAPE is true, the original values will be restored before the
		//!                     function returns. Otherwise the equillibrium constants will remain in the GetStructure()->SHAPE array.
		//! \return An int that indicates an error code (0 = no error, 5 = error reading thermodynamic parameter files).
		int PartitionFunction(const char savefile[]="",double temperature=-10.0, bool disablecoax=false, bool restoreSHAPE=true);
		
		//! Predict structures containing highly probable pairs.

		//! This function predicts structures composed of probable base pairs.
		//! This function requires partition function data from either a previous partition function calculations or
		//!		from having read a partition function save file during construction of the class.
		//!	In case of error, the function returns a non-zero that can be parsed by GetErrorMessage() or GetErrorMessageString().
		//!	\param probability is the pairing probability threshold, where pairs will be predicted if they have a higher probability.  Note that a value of less than 0.5 (50%), will cause an error.  The default value of zero will trigger the creation of 8 structures, with thresholds of >=0.99, >=0.97, >=0.95, >=0.90, >=0.80, >=0.70, >=0.60, >0.50.
		//! \return An int that indicates an error code (0 = no error, non-zero = error).
		int PredictProbablePairsRR(vector<pair<float,pair<int,int> > >&RR_list,
            const float probability=1e-5);

		//! Predict maximum expected accuracy structures that contain pseudoknots from either a sequence or a partition function save file.

		//! This function uses base pair probabilities to predict structures that contains pseudoknots.
		//! This function requires partition function data from either a previous partition function calculations or
		//!		from having read a partition function save file during construction of the class.
		//!	In case of error, the function returns a non-zero that can be parsed by GetErrorMessage() or GetErrorMessageString().
		//!	\param iterations is the number of iterations of pair selection that are performed.  The default and recommended value is 1.
		//!	\param MinHelixLength is the shortest helix that is allowed.  If this is set >1, a post-processing step is performed to remove short helices.  Default = 1, i.e. no post-processing.
		//! \param threshold is the required minimum probability for a pair.  Only include pairs with a higher probability.  Default = 0, i.e. include all pairs.
		//! \return An int that indicates an error code (0 = no error, non-zero = error). 
		int ProbKnot(int iterations=1, int MinHelixLength = 1, double threshold = 0);

		//! Predict maximum expected accuracy structures that contain pseudoknots from a file containing ensemble of structures.

		//! This function uses base pair probabilities to predict structures that contains pseudoknots.
        //! This function requires a file with ensemble of structures. This function processes the file to calculate pair probabilities.
		//!	In case of error, the function returns a non-zero that can be parsed by GetErrorMessage() or GetErrorMessageString().
		//!	\param iterations is the number of iterations of pair selection that are performed.  The default and recommended value is 1.
		//!	\param MinHelixLength is the shortest helix that is allowed.  If this is set >1, a post-processing step is performed to remove short helices.  Default = 1, i.e. no post-processing.
		//! \param threshold is the required minimum probability for a pair.  Only include pairs with a higher probability.  Default = 0, i.e. include all pairs.
		//! \return An int that indicates an error code (0 = no error, non-zero = error). 
		int ProbKnotFromSample(int iterations=1, int MinHelixLength = 1, double threshold = 0 );

		//! Re-predict the lowest free energy secondary structure and generate suboptimal structures using a heuristic.

		//! This function predicts the lowest free energy structure and suboptimal structure after a save file (.sav) was specified to the constructor.
		//! The step of predicting structures from the save file is rapid, so this is laregely a method to quickly generate a different set of suboptimal structures.
		//! Refolding can only be performed if the RNA constructor was called with a save file name.  (That is, you cannot call this after calling fold single strand, without loading the data from disk with a new instance of RNA.  This is for historical reasons.)
		//!	In case of error, the function returns a non-zero that can be parsed by GetErrorMessage() or GetErrorMessageString().
		//!	\param percent is the maximum % difference in free energy in suboptimal structures from the lowest free energy structure.  The default is 20.
		//!	\param maximumstructures is the maximum number of suboptimal structures to generate.  The default is 20.
		//!	\param window is a parameter that specifies how different the suboptimal structures should be from each other (0=no restriction and larger integers require structures to be more different).  The default is 5.
		//! \return An int that indicates an error code (0 = no error, 5 = error reading thermodynamic parameter files, 14 = traceback error).
		int ReFoldSingleStrand(const float percent=20, const int maximumstructures=20, const int window=5);

		//! Sample structures from the Boltzman ensemable.

		//! This function requires partition function data from either a previous partition function calculations or
		//!		from having read a partition function save file during construction of the class.
		//!	In case of error, the function returns a non-zero that can be parsed by GetErrorMessage() or GetErrorMessageString().
		//!	\param structures is the number of structures to be sampled.  The default is 1000.
		//!	\param seed is an integer that seeds the random number generator that is required for sampling, which defaults to 1.
		//! \return An int that indicates an error code (0 = no error, non-zero = error).
		int Stochastic(const int structures=1000, const int seed=1);
		

		//********************************************************************************
		//Functions that specify or report constraints on folding.
		//Also, functions that read or write constraints from disk.
		//These constraints only affect subsequent calls to structure prediction routines.
		//********************************************************************************


		//!Read SHAPE data from disk.
		
		//!The SHAPE data is used to constrain structure prediction on subsequent structure predictions.
		//!The function returns 0 with no error and a non-zero otherwise that can be parsed by GetErrorMessage() or GetErrorMessageString().
		//!Pseudo folding free energy change parameters should be in units of kcal/mol.
		//!\param filename is a c string that indicates a file that contains SHAPE data.
		//!\param IsPseudoEnergy indicates whether this is the pseudo folding free energy constraint (the preferred method).  This defaults to true.
		//!\param slope is the slope when IsPseudoEnergy=true and is a threshold above which nucleotides are forced single stranded otherwise.
		//!\param intercept is the intercept when IsPseudoEnergy=true and is a threshold above which a nucleotide is considered chemically modified otherwise.
		//!\param modifier is the type of chemical modification probe that was used (currently accepted values are SHAPE, diffSHAPE, DMS, and CMCT). Defaults to SHAPE.
		//!\return An integer that indicates an error code (0 = no error, 1 = input file not found).
		int ReadSHAPE(const char filename[], const double slope, const double intercept, RestraintType modifier=RESTRAINT_SHAPE, const bool IsPseudoEnergy=true);

		//!Read SHAPE data from disk including single-stranded SHAPE pseudo free energys.
		
		//!The SHAPE data is used to constrain structure prediction on subsequent structure predictions.
		//!This version of the overloaded function includes a single-stranded pseudo free energy change.
		//!The function returns 0 with no error and a non-zero otherwise that can be parsed by GetErrorMessage() or GetErrorMessageString().
		//!Pseudo folding free energy change parameters should be in units of kcal/mol.
		//!\param filename is a c string that indicates a file that contains SHAPE data.
		//!\param dsSlope is the double-stranded slope.
		//!\param dsIntercept is the double-stranded intercept.
		//!\param modifier is the type of chemical modification probe that was used (currently accepted values are SHAPE, DMS, and CMCT). Defaults to SHAPE.
		//!\param ssSlope is the single-stranded slope.
		//!\param ssIntercept is the single-stranded intercept.
		//!\return An integer that indicates an error code (0 = no error, 1 = input file not found).
		int ReadSHAPE(const char filename[], const double dsSlope, const double dsIntercept, const double ssSlope, const double ssIntercept, RestraintType modifier=RESTRAINT_SHAPE);

		//! Read experimental pair bonuses from disk.

		//! This is a quantity that results in a bonus added to a specific pair, once per stack, so that pairs in the middle of a helix get the bonus twice
		//! and those at the end of a helix get the bonus once.
		//! The bonus is in the form of experimentalScaling*value + experimentalOffset.
		//! The data is formatted using a simple square matrix of values and no headers.  The format requires that there be N^2 entries for a sequence of N nucleotides. 
		//!\param filename is a c string that indicates a file that contains data.
		//!\param experimentalOffset is a double that is added to each value. 
		//!\param experimentalScaling is a double by which each value is multiplied.
		//!\return An integer that indicates an error code (0 = no error, 1 = input file not found).
		int ReadExperimentalPairBonus(const char filename[], double const experimentalOffset, double const experimentalScaling );

		//!Remove all folding constraints.

		//!This function strips all previously assigned folding constraints.
		//!Note that this function does not delete SHAPE constraints or pseudo free energies.
		void RemoveConstraints();
        


		//****************************************
		//Functions that write output information:
		//****************************************


		//*******************************************************
		//Functions that return information about structures:
		//*******************************************************

		//!Break any pseudoknots that might be in a structure.

		//! This function uses the method of Smit et al. to break pseudoknots by running a dynamic programming algorithm that cannot predict pseudoknots while only allowing the pairs that already exist in the structure.
		//! When minimum_energy = true (the default), this function predicts the lowest free energy structure that has no pseudoknots.  
		//! Note that when minimum_energy = true, this function might additionally break pairs that are not pseudoknotted if the pairs increase the folding free energy change or are forbidden (as in an isolated pair).
		//! Also note that when minum_energy=true, this function uses the GenerateAllSuboptimalStructures methodology behind the scenes, so large internal loops in the input would lead to a loss of pairs.
		//! When minumum_energy is set to false, this function maximizes the number of base pairs in the pseudoknot free structure.	
		//!	Return 0 if no error and non-zero errors can be parsed by GetErrorMessage() or GetErrorMessageString().
		//!\param minimum_energy is a bool thgat indicates where the structure should be minimum in free energy (true) or maximize pairs (false).
		//!\param structurenumber is an int that indicates a specific structure for which to break pseudoknots (indexed from 1).  The default value, 0, indicates that all structures should have pseudoknots broken.
		//!\param useFastMethod If true (default) a fast method that is energy and sequence agnostic will be used to find pseudoknots. This method should be identical in results to MEAFill, but much faster. 
		//!                     However, if minimum_energy is also true, AllTrace will be used instead of the fast method (i.e. this parameter will be ignored).
		//!\return an int that provides an error code.  0 = no error.
		int BreakPseudoknot(const bool minimum_energy=true, const int structurenumber = 0, const bool useFastMethod = true);


		//! Get the nucleotide to which the specified nucleotide is paired.

		//! Returns the pairing partner of the ith nucleotide in structure number structurenumber.
		//! Zero means the nucleotide is unpaired.
		//! This function generates internal error codes that can be accessed by GetErrorCode() after the constructor is called: 0 = no error, nonzero = error.
		//! The errorcode can be resolved to a c string using GetErrorMessage.
		//!\param i is an int that indicates the nucleotide to which the pairing partner is being queried.
		//!\param structurenumber is an int that indicates the structure number, where the default is 1.
		//!\return An int that indicates the other nucleotide in pair, where 0 is no paired.
		int GetPair(const int i, const int structurenumber=1);
		
		//!Get the total number of specified or predicted structures.

		//!\return An integer specify the total number of structures.
		int GetStructureNumber() const;


		//*****************************************************
		//Functions that return information about the sequence:
		//*****************************************************

		//!Get the total length of the sequence.

		//!\return An integer that specifies the total length of the sequence.
		int GetSequenceLength() const;

		//!\return The full nucleotide sequence.
		const char* GetSequence() const;

		//!\return A segment of the nucleotide sequence.
		std::string GetSequence(size_t start, size_t length=std::string::npos) const;

		//******************************************************************
		//Functions that access or populate underlying classes:
		//******************************************************************
		
		//!Access the underlying structure class.
		//!This is provided for use with two sequence methods.
		//!Generally, there is no need for end users to use this function because the RNA class provides an convenient wrapper for accessing the information in an RNA class.
		//!\return A pointer to structure.
		structure *GetStructure();

		//******************************************************************
		//Functions that provide a connection to TProgressDialog, for following calculation progress:
		//******************************************************************
		
		//!Return the current pointer to an instance of a class derived from ProgressHandler 
		//!(such as TProgressDialog for text-mode programs or a simple ProgressHandler for GUI progams or scripts)
		//!This is used during inheritance to provide access to the underlying TProgressDialog.
		//\return A pointer to the TProgressDialog class.
		ProgressHandler* GetProgress();

		//*****************************
		//Destructor:
		//*****************************

		//!Destructor.

		//! The destructor automatically cleans up all allocated memory for predicted or specified structures.
		~RNA();


		//!Copy thermodynamic parameters from an instance of an RNA class.

		//!This is generally not needed because functions automatically populate
		//!the parameters from disk.  It is helpful, however, in constructors of derived classes.
		//! Normally the Thermodyanamic(Thermodynamic* copy) constructor should be used instead.
		//!    Reason: Since the extended alphabet, the datatable is loaded in the RNA class constructor, 
		//!    so copying the datatable must also happen in the constructor to avoid reading a 
		//!    (possibly different) datatable first.
		void CopyThermo(Thermodynamics& copy); // override the method in Thermodynamics

	protected:


		//Integer to keep track of error codes.
		//These errors result on file i/o problems during construction.
		//The errors can be accessed using GetErrorCode().
		int ErrorCode;


		//The following are needed to provide calculation progress
		ProgressHandler *progress;

		//Read Files
		int FileReader(const char filename[], const RNAInputType fileType);
		
		//The following set of variables are needed for partition function calculations.
		//they are in protected instead of private so derived classes can access them for probability calculations.
		PFPRECISION *w5,*w3,**wca;
		pfdatatable *pfdata;
		DynProgArray<PFPRECISION> *w,*v,*wmb,*wl,*wmbl,*wcoax, *wlc;
		PFPRECISION Q;

		// Common constuctor logic:
		//   1. set member fields to defaults.
		//   2. If Thermodynamics->currentAlphabetName is not NULL, call ReadThermodynamic to read the datatables
		//   3. If sequenceOrFileName is not NULL:
		//        - If type is SEQUENCE_STRING, load the sequence contained in sequenceOrFileName directly.
		//        - Otherewise load the file specified by sequenceOrFileName.
		void init(const char * sequenceOrFileName, const RNAInputType fileType, const bool allowUnknownBases=false, const bool skipThermoTables=false);

	private:
		//The primitive class for storing sequence and structure data.
		//Private inheritance to force the user to use the interface provided by RNA.
		//call GetStructure() if you need to access the data
		structure *ct;

		//The following bool is used to indicate whether the partion function arrays have been allocated and therefore need to be deleted.
		bool partitionfunctionallocated;

		

		//The following bool is used to indicate whether the folding free energy arrays are allocated and therefore need to be deleted.
		bool energyallocated;
		

		//The following set of variables are used for restoring folding save files (.sav) for refolding and for energy dot plots.
		DynProgArray<integersize> *ew2,*ewmb2;
		integersize *ew5,*ew3;
		int vmin;
		DynProgArray<integersize> *ev,*ew,*ewmb;
		
		//The following variables are used by both partition function calculations and save file restoration.
		bool *lfce,*mod;
		forceclass *fce;

		// Holds error messages that occur when reading files, so additional information can be shown to the users (e.g. in GUI programs where console output is not seen).
		string lastErrorDetails;

};

#endif//RNA_H defined
