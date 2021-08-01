
//Class RNA -- Wrapper for RNAstructure for use in object-oriented applications
// #define ENABLE_DEBUG_LOGS
#include "../src/debug_logging.h"

#include "RNA.h"

#include "../src/rna_library.h"
#include "../src/structure.h"
#include "../src/DynProgArray.h"
#include "../src/forceclass.h"
#include "../src/pfunction.h"
#include "../src/boltzmann.h"
#define PARTITION_CORE

#include <iostream>
#include <iomanip>

const float epsilon = 1e-6; // a small number for a tolerance in comparing floats

// Common constuctor logic:
//   1. set member fields to defaults.
//   2. If alphabetName is not NULL, call ReadThermodynamic to read the datatables
//   3. If sequenceOrFileName is not NULL:
//        - If type is SEQUENCE_STRING, load the sequence contained in sequenceOrFileName directly.
//        - Otherewise load the file specified by sequenceOrFileName.
void RNA::init(const char * sequenceOrFileName, const RNAInputType fileType, const bool allowUnknownBases /* default: false */, const bool skipThermoTables /* default: false */) {
    //set error status to zero
    ErrorCode=0;
    lastErrorDetails = "";

    //allocate ct
    ct = new structure();

    //Indicate that the partition function calculation has not been performed.
    partitionfunctionallocated = false;

    //Indicate that the energy data is not read.
    energyallocated = false;

    //Do not report progress by default:
    progress=NULL;

    // Do not call ReadThermodynamic if any of the following are true:
    //  - The table is already loaded (which would also be the case if Thermodynamics was copied.)
    //  - No alphabet has been specified. (i.e. with the RNA(IsRNA) constructor)
    //  - We are loading a sav file (PFS or FSV/SAV)
    bool readThermo = !(IsAlphabetRead() || GetAlphabetName().empty() || 
                          fileType==FILE_PFS || fileType==FILE_SAV);
    if (readThermo) {
        this->skipThermoTables = skipThermoTables; // set this property in the Thermodynamics base class. If true, only the alphabet will be loaded..Not the energy files.
        //Now read the thermodynamic parameters immediately, so that the alphabet is known
        // Note that this uses the directory, alphabetName, and temperature already set in the inherited Thermodynamics object.
        if ((ErrorCode=ReadThermodynamic())!=0) return; //The thermodynamic parameter files were not found, so report an error.
        // Thermodynamics->data should now be loaded
        data->allowUnknownBases = allowUnknownBases;
    }
    //Now also add the pointer to the data_tables to the structure class:
    if (data!=NULL)ct->SetThermodynamicDataTable(data);

    if (sequenceOrFileName == NULL) return; // do NOT load sequence information if sequenceOrFileName is NULL

    // Load the sequence or file
    if (fileType==SEQUENCE_STRING) { // type==SEQUENCE_STRING indicates that sequenceOrFileName represents a sequence and should be loaded directly.
        // SetSequence returns 0 if successful. Otherwise it returns an error code (and may also set lastErrorDetails)
        ErrorCode = ct->SetSequence(sequenceOrFileName);
    } else {
        // sequenceOrFileName represents the path to a file
        ErrorCode = FileReader(sequenceOrFileName, fileType);
    }
}

//constructor where user provides a string with the sequence
RNA::RNA(const char sequence[], const bool IsRNA):Thermodynamics(IsRNA, IsRNA?DT_RNA:DT_DNA) {
    init(sequence, SEQUENCE_STRING);  // init sets fields to default values, loads the thermodynamic parameters, and then loads the sequence
}

//RNA::RNA(const char sequence[], const char* const alphabetName, const bool allowUnknownBases) {
//  init(sequence, SEQUENCE_STRING, alphabetName, allowUnknownBases);  // init sets fields to default values, loads the thermodynamic parameters, and then loads the sequence
//}

//  This constructor is deprecated. Instead use RNA(const char filename[], const RNAFileType type, const bool IsRNA)
//  For type: 1=>CT_FILEtype, 2 => FILE_SEQ, 3 => FILE_PFS, 4 => FILE_SAV, and 5 => FILE_DBN
RNA::RNA(const char filename[], const RNAInputType type, const bool IsRNA ):Thermodynamics(IsRNA, IsRNA?DT_RNA:DT_DNA) {
    init(filename, (RNAInputType)type);  // init sets fields to default values, loads the thermodynamic parameters, and then loads the sequence
}

RNA::RNA(const char filepathOrSequence[], const RNAInputType type, const char* const alphabetName, const bool allowUnknownBases, const bool skipThermoTables) 
    : Thermodynamics(isAlphabetRNA(alphabetName), alphabetName) {
    init(filepathOrSequence, type, allowUnknownBases, skipThermoTables);  // init sets fields to default values, loads the thermodynamic parameters, and then loads the sequence
}

// Copy Thermodynamics info and datatable from another Thermodynamics instance.
RNA::RNA(const char filepathOrSequence[], const RNAInputType type, const Thermodynamics *copyThermo) : Thermodynamics(*copyThermo) {
    init(filepathOrSequence, type);  // loads the sequence, but does not call ReadThermo if Thermodynamics is already loaded.
}

// Default constructor.
RNA::RNA(const bool IsRNA):Thermodynamics(IsRNA, NULL) {
    // Allocate the underlying structure class and nothing more.
    
    // The Thermodynamics alphabet name was set to NULL, so the thermodynamic parameter tables will NOT be loaded in init.
    // So the programmer is required to call ReadThermodynamic and set the sequence and structural information explicitly.
    // Setting the sequence to NULL in init indicates that no sequence should be loaded.
    init(NULL, SEQUENCE_STRING); 
}

// Override Thermodynamics::CopyThermo so we can copy the datatable to our CT.
void RNA::CopyThermo(Thermodynamics& copy) {
    Thermodynamics::CopyThermo(copy);
    ct->SetThermodynamicDataTable(copy.GetDatatable());
}

//Return the value of ErrorCode
int RNA::GetErrorCode() const {
    return ErrorCode;
}

//Return a c string that describes errors from GetErrorCode and other errors.
const char* RNA::GetErrorMessage(const int error) {
    switch(error) {
        case 0: return "No Error.\n";
        /* Standard File Errors */
        case 1: return "Input file not found.\n";
        case 2: return "Error opening file.\n"; // This should be returned to indicate a file IO error, NOT an error with the content of the file, which should be #28 or #29.
        /* Misc Errors */
        case 3: return "Structure number out of range.\n";
        case 4: return "Nucleotide number out of range.\n";
        case 5: return "Error reading thermodynamic parameters.\n";
        case 6: return "This would form a pseudoknot and is not allowed.\n";
        case 7: return "This pair is non-canonical and is therefore not allowed.\n";
        case 8: return "Too many restraints specified.\n";
        case 9: return "This nucleotide already under a conflicting constraint.\n";
        case 10: return "There are no structures to write to file.\n";
        case 11: return "Nucleotide is not a U.\n";
        case 12: return "Maximum pairing distance is too short.\n";
        case 13: return "Error reading constraint file.\n";
        case 14: return "A traceback error occurred.\n";
        case 15: return "No partition function data is available.\n";
        case 16: return "Wrong save file version used or file format not recognized.\n";
        case 17: return "This function cannot be performed unless a save file (.sav) was correctly loaded by the RNA constructor.\n";
        case 18: return "This threshold is too low to generate valid secondary structures.\n";
        case 20: return "No sequence has been read.\n";
        case 21: return "Probabilities summed to greater than 1 in stochastic traceback.\n";
        case 22: return "Programming error.  Incorrect file type passed to constructor.\n";
        case 23: return "There are no structures present.\n";
        case 24: return "Too few iterations.  There must be at least one iteration.\n";
        case 25: return "Index is not a multiple of 10.\n";
        case 26: return "k, the equilibrium constant, needs to be greater than or equal to 0.\n";
        case 27: return "Lyngso O(N^3) internal loop search is not compatible with a parallel calculation.\n";  
        case 28: return "Error reading sequence.\n";
        case 29: return "Invalid file format.\n";
        case 30: return "Programming error: The thermodynamic parameters have not been read.\n";
        case 31: return "Length mismatch between sequence and annotation file.\n"; // Tried to load probability annotations from a PFS file obtained from folding a shorter sequence.
        case 32: return "Array size mismatch.\n"; // Caller passed in an array that is too small to hold requested information.
        case 33: return "Error opening pseudoknot penalty constants file.\n"; // Caller passed in an array that is too small to hold requested information.
        case 34: return "Error opening output file for writing.\n"; 
		case 35: return "Error writing output file.\n"; 
		case 36: return "Pairs must have probability greater than zero.  Therefore, the probknot threshold must be >= 0.";
        /* SHAPE, Experimental Pair Bonus, etc Restraint data */
        case 201: return "Restraint File Not Found (SHAPE or other experimental data).\n";
        case 202: return "Could not open or read Restraint file (SHAPE or other experimental data).\n";
        case 203: return "Wrong number of restraints in file (SHAPE or other experimental data).\n";
        case 204: return "Invalid Nucleotide Number in Restraint file (SHAPE or other experimental data).\n";
        case 215: return "This function is incompatible with other restraints (SHAPE or other experimental data).\n"; // for Rsample.
        /* General Errors */
        case 99: return "The calculation was canceled.\n";
        default: return "Unknown Error\n";
    }
}

const string RNA::GetErrorDetails() const {
    return lastErrorDetails.empty() ? ct->GetErrorDetails() : lastErrorDetails;
}
void RNA::SetErrorDetails(const string& details) {
    lastErrorDetails = details;
}

void RNA::SetSequenceLabel(const string& label) {
    GetStructure()->SetSequenceLabel(label);
}


int RNA::PredictProbablePairsRR(vector<pair<float,pair<int,int> > >&RR_list,
    const float probability)
{
    int i,j;
    float prob;
    
    if (probability > 1 || probability < 0) return 18;
    if (!partitionfunctionallocated) return 15; //no partition function data

    //Get one clean structure
    int NumberofStructures=ct->GetNumberofStructures();
    if (NumberofStructures>0)
    {
        ct->CleanStructure(1);
        for (i=ct->GetNumberofStructures();i>1;--i)
            ct->RemoveLastStructure();
    }
    else ct->AddStructure();

    //cout<<"probability="<<probability<<endl;
    for (i=1;i<ct->GetSequenceLength();i++)
    {
        for (j=i+1;j<=ct->GetSequenceLength();j++)
        {
            prob=calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce);
            if (prob>probability) RR_list.push_back(
                pair<float,pair<int,int> >(prob,pair<int,int>(i,j)));
        }
    }
    return 0;
}


//Calculate the partition function for the current sequence.
int RNA::PartitionFunction(const char savefile[], double temperature, bool disablecoax, bool restoreSHAPE) {
    int i,j;
    char *savefilename;
    //check to make sure that a sequence has been read
    if (ct->GetSequenceLength()==0) return 20;
    
    if (!VerifyThermodynamic()) return 5; //The thermodynamic data tables have not been read 


    //savefile will be passed to the function dynamic for structure prediction.
    //Dynamic expects a null pointer if no file is to be created, so set savefilename to null if savefile is an empty string.
    if (is_blank(savefile)) 
	    savefilename=NULL;
    else {
        savefilename=new char[((int) strlen(savefile))+1];
        strcpy(savefilename,savefile);
    }

    if (partitionfunctionallocated) {
        delete v;
        delete w;
        delete wmb;
        delete wl;
        delete wlc;
        delete wmbl;
        delete wcoax;
        delete fce;
        delete[] lfce;
        delete[] mod;
        delete[] w3;
        delete[] w5;
        delete pfdata;
    }
    //Allocate the memory needed (only if this is the first call to pfunction):
    //indicate that the memory has been allocated so that the destructor will delete it.
    partitionfunctionallocated = true;

    //allocate space for the v and w arrays:
    w = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
    v = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
    wmb = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
    wl = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
    wlc = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
    wmbl = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
    wcoax = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
    fce = new forceclass(ct->GetSequenceLength());
    lfce = new bool [2*ct->GetSequenceLength()+1];
    mod = new bool [2*ct->GetSequenceLength()+1];

    for (i=0;i<=2*ct->GetSequenceLength();i++) {
        lfce[i] = false;
        mod[i] = false;
    }

    for (i=0;i<ct->GetNumberofModified();i++) {

        if (ct->GetModified(i)!=1&&ct->GetModified(i)!=ct->GetSequenceLength()) {
            mod[ct->GetModified(i)]=true;
            mod[ct->GetModified(i)+ct->GetSequenceLength()]=true;
        }
    }

    w5 = new PFPRECISION [ct->GetSequenceLength()+1];
    w3 = new PFPRECISION [ct->GetSequenceLength()+2];

    if (ct->intermolecular) {
        //take advantage of templating to prevent intramolecular base pairs

        ct->allocatetem();//This allocates a bool array that indicates what pairs are allowed.  It is initilaized to true, i.e. all nucs can pair with each other.
        for (i=1;i<ct->inter[0];i++) {
            for (j=i+1;j<=ct->inter[2];j++) {

                //Set intermolecular pairs to false, i.e. not allowed.  Note the indexing with the high index and then the low index number.
                ct->tem[j][i]=false;

            }
        }
        for (i=ct->inter[2]+1;i<ct->GetSequenceLength();i++) {
            for (j=i+1;j<=ct->GetSequenceLength();j++) {

                //Another set of intermolecular pairs to forbid.
                ct->tem[j][i]=false;

            }
        }
    }

    //Initialize the partition function datatable:
        //Ignore the setting of parameter temperature if it is less than zero.
        //Generally, this parameter should be left at the default.
    
    pfdata=new pfdatatable(data,TO_PFP(scalingdefinition),(temperature<0?GetTemperature():temperature));

    //This code converts the SHAPE array of data to equilibrium constants, which is
    //required for the partition function. If the restoreSHAPE parameter is true 
    // (which it is by default), then a backup copy of the SHAPE data is made and is
    // restored after partition, so that it can be used for Fold, etc.
    double *SHAPE_backup = NULL;
    if (ct->shaped) {
      if (restoreSHAPE) SHAPE_backup = ct->CopySHAPE(false/*includeSHAPEss*/);
      for (i=1;i<=2*ct->GetSequenceLength();i++) ct->SHAPE[i]=(double) boltzman( ct->SHAPE[i], pfdata->pftemp);
    }

    // Converts Experimental Pair Bonus array to equilibrium constants
    if (ct->experimentalPairBonusExists){
      for (i = 1; i <= 2*ct->GetSequenceLength(); i++){
        for(j = i; j <= 2*ct->GetSequenceLength(); j++){
          double avg = 0.5*(ct->EX[i][j] + ct->EX[j][i]);
          ct->EX[i][j] = boltzman(avg, pfdata->pftemp);
          ct->EX[j][i] = ct->EX[i][j]; //symmetrize as well...
        }
      }
    }

    //The next section handles the case where base pairs are not
            //not allowed to form between nucs more distant
            //than ct->GetPairingDistanceLimit()
    if (ct->DistanceLimited()) {
        //This allocates a bool array that indicates what pairs are allowed.  It is initilaized to true, i.e. all nucs can pair with each other.
        if (!ct->templated) ct->allocatetem();
        for (j=minloop+2;j<=ct->GetSequenceLength();j++) {
            for (i=1;i<j;i++) {
                //Set distant pairs to false, i.e. not allowed.  Note the indexing with the high index and then the low index number.
                if (j-i>=ct->GetPairingDistanceLimit()) ct->tem[j][i]=false;
            }
        }
    }

#ifndef _CUDA_CALC_
    //cout << "Performing CPU Partition Function Code" << endl;
    //default behavior: calculate the partition function on the CPU
    calculatepfunction(ct,pfdata,progress,savefilename,false,&Q,w,v,wmb,wl,wlc,wmbl,wcoax,fce,w5,w3,mod,lfce,disablecoax);
#else //ifdef _CUDA_CALC_
    //if cuda flag is set, calculate on GPU
    //this requires compilation with nvcc
    //it will work for pair probabilities but it will break stochastic traceback

//read nearest neighbor parameters at desired temperature
    const char *path = getDataPath();
    if (!path)
        die("%s: need to set environment variable $DATAPATH", "Could not find nearest neighbor parameters");
    else
        cout << endl << "Read in data tables" << endl;
    struct param par;
    //param *par = new param();
    //cout << "Size (cpp) : " << sizeof(*par) << endl;
    //cout << "Size (cpp int22) : " << sizeof(par->int22) << endl;
    

    param_read_from_text(path, &par, 0,0);
//    param_read_alternative_temperature(path, &par, pfdata->pftemp, !isrna);
//copy the sequence, correcting for 1-indexing
    
    char* buf = new char[ct->GetSequenceLength()+1];
    for(int i=0;i<GetSequenceLength();i++){
        buf[i] = ct->nucs[i+1];
    }
    buf[ct->GetSequenceLength()] = '\0'; //make it null-terminated
//calculate the partition function on the GPU
    int *bcp;
    bcp = ct->generate_constraint_matrix();

    prna_t p = prna_new(buf, &par, 1, bcp);
//transfer over the partition function from the GPU data structure
//because the GPU partition function is calculated in log space, we do not directly transfer the partition function
//instead, we put the pair probs, calculated in log space in V, 1 everywhere else
//so RNA::GetPairProbability(i,j) will return V(i,j)*1/1, which will be the pairing probability
    for(int i=1;i<=ct->GetSequenceLength();i++){
        for(int j=i+1;j<=ct->GetSequenceLength();j++){
            //cout << "i:\t" << i << "\tj:\t" << j <<endl;
            // Get V(i,j)
            v->f(i,j) = XLOG_TO_PF(-1*get_v_array(p, i-1, j-1)); //V(i,j) = P(i,j)
            
            //Get V'(i,j)
            //In the partition-cuda code, V'(i,j) is stored in V(j,i)
            v->f(j,i+ct->GetSequenceLength()) = XLOG_TO_PF(-1*get_v_array(p, j-1, i-1)); //V'(i,j) = 1
        }
        //cout << "i:\t" << i << endl;
        w5[i] = PFSCALE(XLOG_TO_PF(-1*get_w5_array(p, i-1)), pfdata->scaling, -2); //have to divide by scaling*2
        w3[i] = PFSCALE(XLOG_TO_PF(-1*get_w3_array(p, i-1)), pfdata->scaling, -2); //because RNA::calculateprobability expects scaling values
    }
    delete[] buf;
    prna_delete(p);
    delete[] bcp;
    //delete &par;
#endif
	if (savefilename!=NULL) {
        if (!progress||!progress->canceled()) writepfsave(savefilename,ct,w5,w3,v,w,wmb,wl,wlc,wmbl,wcoax,fce,mod,lfce,pfdata);

        //clean up some memory use:
        delete[] savefilename;
    }
    if (SHAPE_backup!=NULL) { ct->LoadSHAPE(SHAPE_backup, false/*includeSHAPEss*/); delete[] SHAPE_backup; }
    return progress&&progress->canceled() ? 99 : 0;
}


//Read SHAPE data to constrain structure prediction on subsequent structure predictions.
//filename is a c string that indicates a file that contains SHAPE data.
//IsPseudoEnergy indicates whether this is the pseudo folding free energy constraint (the preferred method).  This defaults to true.
//slope is the slope when IsPseudoEnergy=true and is a threshold above which nucleotides are forced single stranded otherwise.
//intercept is the intercept when IsPseudoEnergy=true and is a threshold above which a nucleotide is considered chemically modified otherwise.
//modifier is the type of chemical modification probe that was used (currently accepted values are SHAPE, diffSHAPE, DMS, and CMCT). Defaults to SHAPE.
//Returns an integer that indicates an error code (0 = no error, 1 = input file not found).
int RNA::ReadSHAPE(const char filename[], const double slope, const double intercept, RestraintType modifier, const bool IsPseudoEnergy) {
    // ct->ReadSHAPE will verify that the SHAPE input file exists
    int code;
    if (IsPseudoEnergy) {
        //This is the pseudo energy version
        ct->SHAPEslope=slope*conversionfactor;//register the slope in tenths of kcal/mol
        ct->SHAPEintercept=intercept*conversionfactor;//register the intercept in tenths of a kcal/mol
        code = ct->ReadSHAPE(filename, modifier);//call ReadSHAPE() to read the file and determine pseudo energies
    }
    else
        code = ct->ReadSHAPE(filename,(float) slope,(float) intercept);//call ReadSHAPE() with parameters to parse thresholds
    if (ErrorCode==0) ErrorCode=code; // set the error code, but ONLY if it hasn't already been set (we don't want to hide a previous error)
    return code;
}

//Read SHAPE data to constrain structure prediction on subsequent structure predictions - overloaded version for including single-stranded SHAPE.
//filename is a c string that indicates a file that contains SHAPE data.
//dsSlope is the double-stranded slope.
//dsIntercept is the double-stranded intercept.
//modifier is the type of chemical modification probe that was used (currently accepted values are SHAPE, DMS, and CMCT). Defaults to SHAPE.
//ssSlope is the single-stranded slope.
//ssIntercept in the single-stranded intercept.
//Returns an integer that indicates an error code (0 = no error, 1 = input file not found).
int RNA::ReadSHAPE(const char filename[], const double dsSlope, const double dsIntercept, const double ssSlope, const double ssIntercept, RestraintType modifier) {
    // ct->ReadSHAPE will verify that the SHAPE input file exists

    ct->SHAPEslope=dsSlope*conversionfactor;//register the slope in tenths of kcal/mol
    ct->SHAPEintercept=dsIntercept*conversionfactor;//register the intercept in tenths of a kcal/mol
    ct->SHAPEslope_ss=ssSlope*conversionfactor;//register the slope in tenths of kcal/mol
    ct->SHAPEintercept_ss=ssIntercept*conversionfactor;//register the intercept in tenths of a kcal/mol
    int code = ct->ReadSHAPE(filename, modifier);//call ReadSHAPE() to read the file and determine pseudo energies
    if (ErrorCode==0) ErrorCode=code; // set the error code, but ONLY if it hasn't already been set (we don't want to hide a previous error)
    return code; 
}

//Read SHAPE data to constrain structure prediction on subsequent structure predictions.
//filename is a c string that indicates a file that contains SHAPE data.
//parameter1 is the slope.
//parameter2 is the intercept.
//Returns an integer that indicates an error code (0 = no error, 1 = input file not found).
int RNA::ReadExperimentalPairBonus(const char filename[], double const experimentalOffset, double const experimentalScaling ) {
    // ct->ReadExperimentalPairBonus will verify that the input file exists
    int code = ct->ReadExperimentalPairBonus(filename, experimentalOffset, experimentalScaling );
    if (ErrorCode==0) ErrorCode=code; // set the error code, but ONLY if it hasn't already been set (we don't want to hide a previous error)
    return code; 
}


//Remove all previously defined constraints.
void RNA::RemoveConstraints() {
    ct->RemoveConstraints();
    ct->min_gu=0;
    ct->min_g_or_u=0;
    ct->nneighbors=0;
    ct->nregion=0;
    ct->nmicroarray=0;
}


// Get the nucleotide to which the specified nucleotide is paired.
int RNA::GetPair(const int i, const int structurenumber) {

    //check to make sure i is a valid nucleotide
    if (i<1||i>ct->GetSequenceLength()) {
        ErrorCode=4;
        return 0;
    }
    //make sure there is structure information for this RNA
    else if (ct->GetNumberofStructures()==0){
        ErrorCode = 23;
        return 0;
    }
    //now make sure structurenumber is a valid structure
    else if (structurenumber<1||structurenumber>ct->GetNumberofStructures()) {
        ErrorCode=3;
        return 0;
    }
    else {
        //error trapping complete
        return ct->GetPair(i,structurenumber);
    }
}



//Get the total length of the sequence
int RNA::GetSequenceLength() const {
    return ct->GetSequenceLength();
}

const char* RNA::GetSequence() const {
    return ct->GetSequence();
}

std::string RNA::GetSequence(size_t start, size_t length) const {
    if (start<1) start=1;
    if (start>GetSequenceLength()) return "";
    if (length==string::npos)
        length=GetSequenceLength()-start;
    length = std::min(length, GetSequenceLength()-start);
    return string(ct->nucs+start, length);
}

//Access the underlying structure class.
structure *RNA::GetStructure() {
    return ct;
}


RNA::~RNA() {





    if (partitionfunctionallocated) {
        //The partition function calculation was performed, so the memory allocated for the partition function needs to be deleted.
        delete[] lfce;
        delete[] mod;
        delete[] w5;
        delete[] w3;
        delete v;
        delete w;
        delete wmb;
        delete wl;
        delete wmbl;
        delete wcoax;
        delete fce;
        delete pfdata;

    }

    if (energyallocated) {
        //A folding save file was opened, so clean up the memory use.

        delete[] lfce;
        delete[] mod;

        delete[] ew5;
        delete[] ew3;

        if (ct->intermolecular) {
            delete ew2;
            delete ewmb2;
        }

        delete ev;
        delete ew;
        delete ewmb;
        delete fce;
    }

    delete ct;//delete the structure
}

// string getFileNameWithoutExt(const char*)


//This is a protected function for handling file input.
int RNA::FileReader(const char filename[], const RNAInputType type) {

    if (!isStdIoFile(filename) && !fileExists(filename)) {
        SetErrorDetails(sfmt("The path '%s' is invalid or does not exist.", filename));
        return 1; // file not found.
    }

    if (type==FILE_CT||type==FILE_SEQ||type==FILE_DBN)
        if (!IsAlphabetRead()) return 30;

    // RMW 2015-03-12: try/catch to prevent any unexpected errors from crashing the program. At least show an error message.
    // (Previous to this, passing the wrong type of file to openct caused an unhandled memory allocation exception.)
    try {
        //open the file based on type:
        switch(type) {
            case FILE_CT: // ct file
                return ct->openct(filename);  // openct returns 0 on success and uses the same non-zero error codes as RNA::GetErrorMessage. So the result can be returned directly.
            case FILE_DBN: // dot-bracket file
                return ct->opendbn(filename); // opendbn returns 0 on success and uses the same non-zero error codes as RNA::GetErrorMessage. So the result can be returned directly.
            case FILE_SEQ:
                //type indicates a .seq file
                return ct->openseqx(filename);
            case FILE_PFS: { //partition function save file
                LOG_INFO("Reading PFS file")
                short vers;
                //allocate the ct file by reading the save file to get the sequence length:
                std::ifstream sav(filename,std::ios::binary);

                read(&sav,&(vers));//read the version of the save file

                if (vers!=pfsaveversion) {
                    //Wrong version!
                    sav.close();
                    return 16;
                }

                int sequencelength;
                //read the length of the sequence
                read(&sav,&(sequencelength));
                sav.close();

                LOG_INFO("Allocate Arrays")
                //allocate everything
                ct->allocate(sequencelength);

                w = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
                v = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
                wmb = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
                wmbl = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
                wcoax = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
                wl = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
                wlc = new DynProgArray<PFPRECISION>(ct->GetSequenceLength());
                fce = new forceclass(ct->GetSequenceLength());

                w5 = new PFPRECISION [ct->GetSequenceLength()+1];
                w3 = new PFPRECISION [ct->GetSequenceLength()+2];

                lfce = new bool [2*ct->GetSequenceLength()+1];
                mod = new bool [2*ct->GetSequenceLength()+1];

                pfdata = new pfdatatable();
				data = new datatable();

                //indicate that the memory has been allocated so that the destructor will delete it.
                partitionfunctionallocated = true;
                LOG_INFO("Reading save data")

                //load all the data from the pfsavefile:
                readpfsave(filename, ct, w5, w3,v, w, wmb,wl, wlc,wmbl, wcoax, fce,&pfdata->scaling,mod,lfce,pfdata,data);
                return 0;
            }  // FILE_PFS
            default:
                return 22; // error - invalid file type
        } // SWITCH type
    } catch (std::exception* ex) {
        SetErrorDetails(ex->what());
        return 2;
    }
}
