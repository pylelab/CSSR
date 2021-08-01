


//#if !defined(STRUCTURE_CPP)
//#define STRUCTURE_CPP

//#include <stdlib.h>
#include "structure.h"
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "defines.h"
#include "pfunction_math.h"

using namespace std;

#define DBN_BRACKET_SYMBOLS "()<>{}[]AaBbCcDd" // List of characters that represent basepairs in dot-bracket notation. The symbols are listed in pairs: <open1><close1><open2><close2>....
#define DBN_UNPAIRED_SYMBOLS ".-," // List of characters accepted as unpaired bases in dot-bracket notation.

#if defined(WIN32) && !defined(__GNUC__)
	//This is compiling on Windows, where there is no definition for tgamma in cmath
	//Instead, use mathimf.h by including a separate file
	//Note: If compiling on Windows with GCC (e.g. MinGW) we do NOT need to include gamma.h

	#include "gamma.h"
	#define tgamma gammafunction
	//#define isnan _isnan
#endif

#ifdef COUNTING
	#define NO_COUNT(...) __VA_ARGS__._value
#else
	#define NO_COUNT(...) __VA_ARGS__
#endif //COUNTING

//Size the contained vectors in singlestructure for the length of the sequence:
singlestructure::singlestructure(int sequencelength) :
	basepr(sequencelength+1) //Add one to the sequence length so that the arrays can be one-indexed:
{
	//Initialize the energy of a new structure to 0.
	energy = 0;

}

//Establish an operator for comparing singlestructures bases on folding free energy change:
inline bool operator<(const singlestructure &a, const singlestructure &b) {

	return (a.energy < b.energy);
}


//Constructor
structure::structure(int structures)
{

	
	int i;

	allocatestructure(structures);


	allocated = false;

	
	
	intermolecular = false;
	
	templated = false;
	
	min_gu = 0;
	min_g_or_u = 0;
	nneighbors = 0;
	nregion = 0;
	nmicroarray=0;
	

	//initialize values for SHAPE slope and intercept as zero for both double and single stranded restraints.
	SHAPEslope_ss = 0;
	SHAPEintercept_ss = 0;

	SHAPEslope = 0;
	SHAPEintercept = 0;

	data = NULL;  //This is the pointer to the thermodynamic data, stored in a datatable.


	for (i=0;i<maxregions;i++) rnneighbors[i]=0;
	//for (i = 0; i < 100;i++) {
	//	for ( j = 0; j < 25; j++) {
	//		neighbors[i][j]=0;
	//	}
	//}

	stacking = false;
	limitdistance=false;//toogle to true to limit base pairing distance
	maxdistance=600;//default maximum distance between paired nucs
	shaped = false;//by default, a structure does not have SHAPE data associated with it
	SHAPE=SHAPEss=NULL;
	SHAPEss_region = NULL;//by default, the SHAPEss_region does not need to be allocated
	ssoffset = false;//by default, a structure does not have a single stranded offset read from disk
	experimentalPairBonusExists = false;//by default, no pairwise bonuses provided
	EX=NULL;

	constant = NULL;//by default, indicate that the array constant is not allocated


	numofbases = 0;//To start, set the length of the sequence to zero.  This indicates that no sequence has been read.
	distsread = false;

	sequencelabel="\n";//define a default sequence label

	lastErrorDetails = "";

}


//Get an energy associated with a structure
int structure::GetEnergy(int structurenumber) const {

	return arrayofstructures[structurenumber-1].energy;


}

// Get the number of structures stored.
int structure::GetNumberofStructures() const {

	return (int) arrayofstructures.size();

}

//Get the pairing partner of i in structure structurenumber
int structure::GetPair(int i, int structurenumber) const {
	VERIFY_STRUCTURE_INDEX(structurenumber)
	VERIFY_NUC_INDEX(i)
	//Use structurenumber-1 so that the underlying zero indexing works with the 1 indexing used in the rest of the software.
	return arrayofstructures[structurenumber-1].basepr[i];

}

//Get the label associated with the sequence.
string structure::GetSequenceLabel() const {

	return sequencelabel;

}

const char* structure::GetSequence() const {
    if (GetSequenceLength()==0)
        return "";
    return nucs+1; // add one char position, so the first (unused) char is not included.
}



// A stack useful for pseudoknot calculations. It stores intervals [i, j].
// An interval is pushed onto the stack with push(i', j')  (which does NOT affect the members i or j)
// An interval is popped off the stack with pop(), which sets the members i and j.
struct IntervalStack {
	vector<short> stack; unsigned short i, j; unsigned int pos;  /** i and j are the Current values */
	IntervalStack(int capacity = 8) :  stack(capacity), pos(0) {} 
};

//Set the label for structure structurenumber
void structure::SetCtLabel(const string &label, const int structurenumber) {

	arrayofstructures[structurenumber-1].ctlabel = label;

}

//A version of lable setting that uses char pointer:
void structure::SetCtLabel(const char *label, const int structurenumber) {

	arrayofstructures[structurenumber-1].ctlabel = label;

}

//Assign the energy for a structure:
void structure::SetEnergy(int structurenumber, int energy) {
	
	arrayofstructures[structurenumber-1].energy = energy;

}


//Specify a basepair in a specific structure:
void structure::SetPair(int i, int j, int structurenumber) {

	arrayofstructures[structurenumber-1].basepr[i] = j;
	arrayofstructures[structurenumber-1].basepr[j] = i;

}


//Set a sequence label.
void structure::SetSequenceLabel(const string& label) {

	sequencelabel = label;

}


// Set the sequence for this structure.
// This includes the following operations:
//   - Allocate space for the bases.
//   - Verify each base exists in the alphabet.
//   - Setup the arrays nucs, numseq, and hnumber.
//   - Check to see if any nucleotide needs to be single-stranded and call AddSingle for them.
// The sequence can contain whitespace, which is ignored.
// If an error occurs (such as the data-table has not yet been loaded) the return value 
// will be an RNA error code (i.e. it corresponds to one of those defined in RNA::GetErrorMessage)
// Additionally lastErrorDetails may be set (which can be queried with GetErrorDetails())
int structure::SetSequence(const string& sequence) {
	if (!IsAlphabetLoaded()) return 30;
	int count = 0; string::const_iterator it, ite=sequence.end();
	//iterate through the string and count the number of non-whitespace characters
	for(it=sequence.begin();it!=ite;it++) if (!::isspace(*it)) count++; // skip whitespace

	allocate(count);
	count=0;
	nucs[0]='\0';
	hnumber[0]=0;
	// Setqup the CT sequence (i.e. arrays nucs, numseq, numofbases)
	for (unsigned int i=0;i<sequence.size();i++) {
		// record each nucleotide
		char base=sequence[i]; // sequence indices are 0-based while the CT arrays (and pos) are 1-based, so subtract 1 from pos.
		if (::isspace(base)) continue;
		nucs[++count]=base;
		int num = data->basetonum(base); 
		if (num==-1) {
			SetErrorDetails(sfmt("Invalid nucleobase %c at position %i.", base, i+1));
			return 28; // error reading sequence.
		}
		numseq[count]=num;
		hnumber[count]=count;
		//Check to see if any nucleotide needs to be single-stranded:
		//scan through the notpairing list for nucleotides that cannot pair
		for (unsigned int j=0;j<data->not_pairing.size();++j) {
			if (base==data->not_pairing[j]) {
				AddSingle(count); 
				break;
			}
		}
	}
	nucs[numofbases+1]='\0'; // add terminating null character.
	return 0; //success
}

//This allocates space in an array that is used for folding with phylogenetic data.
//	tem == template for allowed and disallowed pairs
void structure::allocatetem()

{
	int i,j;
	//Size = size;//save the size of the array so that the destructor can
   				//deallocate the space

   tem = new bool *[numofbases+1];
   for (i=0;i<=numofbases;i++) {
    	tem[i] = new bool [i+1];
   }
   templated = true;

   //initialize all positions to true:
	for (i=0;i<=numofbases;i++) {
		for (j=i;j<=numofbases;j++) {
    		tem[j][i] = true;
		}
   }

}

// Add a nucleotide to the list of those that must pair.
void structure::AddDouble(int i) {

	doublestranded.push_back(i);
}


//Add a nucleotide to the list of Us in GU pairs.
void structure::AddGUPair(int i) {

	GUpair.push_back(i);

}

//Add a nucleotide to the list of those accessible to traditional chemical modification.
void structure::AddModified(int i) {

	modified.push_back(i);

}

//Add a pair of nucleotides to the list of those that must form.
void structure::AddPair(int i, int j) {

	pair5.push_back(i);
	pair3.push_back(j);

}

//Add a nucleotide to the list of those not able to pair. 
void structure::AddSingle(int i) {

	singlestranded.push_back(i);

}

//! Get a nucleotide that must be base paired.
int structure::GetDouble(int i) {

	return doublestranded[i];

}

//Get a nucleotide that must not be in a specific pair.
int structure::GetForbiddenPair5(int i) {

	return forbid5[i];

}


//Get a nucleotide that must not be in a specific pair.
int structure::GetForbiddenPair3(int i) {


	return forbid3[i];

}

//Get a nucleotide that must be a U in a GU pair.
int structure::GetGUpair(int i) {

	return GUpair[i];

}

//Get a nucleotide that is accessible to chemical modification.
int structure::GetModified(int i) {

	return modified[i];

}

//Get the number of nucleotides forced to be double-stranded.
int structure::GetNumberofDoubles() {

	return doublestranded.size();

}

//Get the number of pairs that are forbidden.
int structure::GetNumberofForbiddenPairs() {

	return forbid5.size();

}
		
//Get the number of Us forced to be in GU pairs.
int structure::GetNumberofGU() {

	return GUpair.size();

}
		
//Get the number of nucleotides that are accessible to chemical modification.
int structure::GetNumberofModified() {

	return modified.size();

}
		
//Get the number of nucleotides forced to be single-stranded.
int structure::GetNumberofSingles() {

	return singlestranded.size();

}
		
//Get the number of pairs that are constrained to occur.
int structure::GetNumberofPairs() {

	return pair5.size();

}

//Get the number of folding domains.
int structure::GetNumberofDomains() {

	return domains5.size();

}

//Get a nucleotide that must be in a specific pair.
int structure::GetPair5(int i) {

	return pair5[i];

}

//Get a nucleotide that must be in a specific pair.
int structure::GetPair3(int i) {

	return pair3[i];

}


//Get a nucleotide that defines a domain.
int structure::GetDomain5(int i) {

	return domains5[i];

}

//Get a nucleotide that defines a domain.
int structure::GetDomain3(int i) {

	return domains3[i];

}


//Get a nucleotide that must be single stranded.
int structure::GetSingle(int i) {

	return singlestranded[i];

}



//Reset, i.e. remove, all constraints
void structure::RemoveConstraints() {

	doublestranded.clear(); 
	singlestranded.clear();
	GUpair.clear();
	modified.clear();
	pair5.clear();
	pair3.clear();
	forbid5.clear();
	forbid3.clear();

}



//Open a Dot-Bracket (DBN) Structure File from disk and store all the information.
//Note that the thermodynamic data must be read from disk and SetThermodynamicDataTable used to set the datatable pointer
//before this function is called.
// This can parse dot-bracket files in any of the following formats (see DotBracketFormat): 
//  DBN_FMT_SINGLE_TITLE,  DBN_FMT_SIDE_TITLES, or DBN_FMT_MULTI_TITLE
//  Note: DBN_FMT_MULTI_TITLE_AND_SEQ is NOT supported because a structure class can only contain a single sequence.
int structure::opendbn(const char *bracketFile) {
	if (!IsAlphabetLoaded()) return 30; //error -- parameters have not been read.
	// Create a variable that handles errors.
	long linenumber = 0;
	string line;
	string sequence;
	size_t start = string::npos;

	std::istream input(std::cin.rdbuf()); // if the filename is "-", use STDIN for input.
	std::ifstream file_in;
	if (!isStdIoFile(bracketFile)) { 
		if (!fileExists(bracketFile)) {
			SetErrorDetails(sfmt("The path '%s' is invalid or does not exist.", bracketFile));
			return 1; // file not found.
		}
		//Open the file for reading:
		file_in.open(bracketFile);
		//Make sure the file opened.  If not, return an error indicator.
		if (!file_in.is_open()) {
			SetErrorDetails(sfmt("Failed to open the file '%s'. Please verify the file and its permissions.", bracketFile));
			return 2; // Error opening file
		}
		input.rdbuf(file_in.rdbuf());
	}

	//scan for first relevant line. Ignore lines that are empty, contain only whitespace, or start with a semi-colon (;).
	#define WHITESPACE " \t\r"
	while (getline(input,line)) {
		++linenumber;
		if (line.length() == 0 || line[0] == '\r') continue; // line is empty
		if (line[0] == ';')
			continue;
		if (string::npos != (start = line.find_first_not_of(WHITESPACE))) // if the line has any non-whitespace character, exit the loop
			break;
	}
	if (input.bad()) {
		//Some kind of read error occurred.
		SetErrorDetails(sfmt("An error occured while reading the file '%s'.", bracketFile));
		return 2;
	}

	// Now past the comments. See if a valid line was found.
	if (start == string::npos || line[start] != '>') {
		//The end of the file was reached before finding the end of the comments, return with error
		SetErrorDetails("The dot-bracket file did not contain a sequence label starting with '>'.");
		return 29;
	}

	line.erase(0, start+1); // remove all characters up to and including start. (they were either whitespace or '>')
	SetSequenceLabel(line);

	// Read the next line, which is the sequence.
	if (getline(input,line).fail()) {
		SetErrorDetails("The dot-bracket file did not contain a sequence.");
		return 29;
	}
	linenumber++;

	//Process each character in the sequence
	for (unsigned int i=0;i<line.length();++i) {
		if (line[i] < 33) continue; // skip whitepace.
		// If the character is valid, add it. Otherwise return with an error.
		if (-1==data->basetonum(line[i])) {
			//This nucleotide wasn't recognized
			// report this error via the lastErrorDetails mechanism
			SetErrorDetails(sfmt("Invalid nucleobase '%c' at line %li column %i.", line[i], linenumber, i+1));
			return 28;
		} else {
			//Put the nucleotide in the sequence
			sequence+=line[i];
		}
	}//end of iteration over sequence characters

	if (sequence.size()==0) {
		SetErrorDetails("The file did not contain any nucleotides.");
		return 29;
	}

	SetSequence(sequence);
	
	// Read subsequent lines. Each is a new structure. e.g.: 
	//   >SequenceLabel
	//   AGUGACU...
	//   ((..(((...  (Str1 comment)
	//   ..(((((...  (Str2 comment)
	//   ((((..((..  (Str3 comment)
	int count=0; // number of structures
	string label("");
	while (!getline(input,line).fail()) {
		++linenumber;
		if (line.find_first_not_of(WHITESPACE)==string::npos) continue; // line is empty or whitespace

		if (line[0] == ';')
			continue; // a comment

		if (line[0] == '>') { // A structure label.
			label = line.substr(1);
			continue;
		}

		//Add a structure to which pairs can be placed.
		count++;
		AddStructure();
		SetCtLabel(label.empty()?GetSequenceLabel():label, count);
		label.clear();

		//Now parse the brackets.
		// Multiple bracket types can be used. Pseudoknots can be encoded by using different types of brackets: e.g.  ..<<<..(...>>>..)..
		struct Bracket { unsigned int id; unsigned int pos; };
		#define INVALID_BRACKET_ID (unsigned int)(-1)  // get the max unsigned int.
		int bracketCount = 0;
		Bracket *brackets= new Bracket[line.size()];
		
		//DEBUG: cout << "Starting bracket parsing." << endl;
		//DEBUG: cout << "Seq:" << sequence << endl;
		//DEBUG: cout << "Brk:" << line << endl;

		for(unsigned int i=0; i<line.size();i++) {
			char c = line[i];
			if (c=='\r'||findchr(DBN_UNPAIRED_SYMBOLS, c)!=string::npos) continue; // skip if it is a CR or a dot-symbol (. or -).
			if (c==' '||c=='\t') {
				//whitespace ends the structure. the remainder is a sequence label
				label=line.substr(i+1);
				cout << "Found side label: '" << label << "'" << endl;
				trim(label);
				if (!label.empty())
					SetCtLabel(label, count);
				break; // quit parsing the structure since the side-label has been found.
			}

			 size_t id = findchr(DBN_BRACKET_SYMBOLS, c);
			if (id==string::npos) {
				SetErrorDetails(sfmt("Invalid character in dot-bracket file: '%c' at line %i column %i.", c, linenumber, i+1));
				return 29;
			}
			if ((id&1)==0) {
				// it is a left-side (opening) bracket. Just add it to the list.
				//DEBUG: cout << "LeftBracket " << c << endl;
				brackets[bracketCount].id=id+1; // +1 because we store the closing bracket id, which is one greater than the opening bracket.
				brackets[bracketCount].pos=i;
				bracketCount++;
			} else {
				// it is a right-side (closing) bracket.
				// find its match and add the pair.
				//DEBUG: cout << "RightBracket " << c << endl;
				int found = -1;
				for(int j=bracketCount-1;j>-1;j--) { // look backwards through the pairs vector, searching for the matching bracket.
					if (brackets[j].id==id) {
						found = j; break;
					}
				}
				if (found==-1) {
					SetErrorDetails(sfmt("Unmatched bracket in dot-bracket file: '%c' at line %li column %i.", c, linenumber, i+1));
					return 29;
				}
				//DEBUG: cout << "SetPair " << (brackets[found].pos+1) << "," << (i+1) << endl;
				// Register the pair that was found.
				SetPair(brackets[found].pos+1,i+1, count);  // +1 because pos is the 1-based nucleotide position, while i is the 0-based string index.

				// If we are at the end of the array, decrement the count
				if (found==bracketCount-1) {
					bracketCount--;
					while(bracketCount>0 && brackets[bracketCount].id==INVALID_BRACKET_ID) // this loop removes all invalidated brackets that might have been introduced by a pseudoknot (i.e. when pairs of different-typed-brackets cross each other).
						bracketCount--;
				} else
					// this bracket is NOT at the end. It must be part of a pseudoknot. Just invalidate it for now.
					brackets[found].id = INVALID_BRACKET_ID;
			}
		} // end of parsing brackets
		delete[] brackets;
	} // end of getLine loop
	
	if (GetNumberofStructures()==0) {
		SetErrorDetails("The dot-bracket file did not contain any structures.");
		return 29;
	}
	return 0; // success
}

// Open a CT File from disk and store all the information.
// Note that the thermodynamic data must be read from disk and SetThermodynamicDataTable 
// used to set the datatable pointer before this function is called.
// int structure::openct(const char *ctfile) {
// 	// if ctfile is "-", read from stdin.
// 	if (isStdIoFile(ctfile))
// 		return openct(std::cin, "STDIN");
// 	else {
// 		if (!fileExists(ctfile)) return 1;
// 		std::ifstream fin(ctfile);
// 		//Open the file for reading:
// 		fin.open(ctfile);
// 		if (!fin.is_open()) return 2;
// 		return openct(fin, ctfile);
// 		//Make sure the file opened.  If not, return an error indicator.
// 	}
// }

//Open a CT File from disk and store all the information.
//Note that the thermodynamic data must be read from disk and SetThermodynamicDataTable used to set the datatable pointer
//before this function is called.
int structure::openct(const char *ctfile) {
	if (!IsAlphabetLoaded()) return 30; //error -- parameters have not been read.
	const int MAX_CT_BASES = 20000; // arbitrary number used mostly to prevent huge memory allocations if the CT file has text/garabage at the top that looks like a large number.
	int count, basepair;
	char base[2];
	string sequenceLabel;
	
	std::istream in(std::cin.rdbuf()); // if the filename is "-", use STDIN for input.
	std::ifstream file_in;
	if (!isStdIoFile(ctfile)) { 
		if (!fileExists(ctfile)) return 1;
		//Open the file for reading:
		file_in.open(ctfile);
		//Make sure the file opened.  If not, return an error indicator.
		if(!file_in.is_open()) return 2;
		in.rdbuf(file_in.rdbuf());
	}

	// detect whether this is a dot-bracket or a CT file.
	if (in.peek()=='>') {
		// if the first character is '>' this is a dot-bracket file
		file_in.close(); // close the input stream (but only if it was a file -- not stdin) before calling opendbn. Ideally we could just pass our ifstream to an overload of opendbn
		return opendbn(ctfile);
	}

	//First read the first item in the file.  If this is -100, then the file was flagged as a CCT-formatted file.  (A more compact format.)
	in >> count;
	if (!in) { // test the istream to make sure a valid number was read (and not just text/garbage etc)
		SetErrorDetails("Invalid character data at line 1 column 0. Expected a valid number of bases.");
		return 29;
	}

	if (count == -100) { 
		//this is a CCT formatted file:
		in >> numofbases;
		in >> count;
		getline(in, sequenceLabel);
		SetSequenceLabel(sequenceLabel);
		allocate(numofbases);
		for (int i=1;i<=numofbases;i++) {
   			in >> numseq[i];
			nucs[i]=GetThermodynamicDataTable()->numtobase(numseq[i]);
			// Todo: check for invalid base (-1)
			hnumber[i] = i;
		}
		for (int numofstructures=1;numofstructures<=count;numofstructures++) {
			AddStructure();
			SetCtLabel(sequenceLabel,numofstructures);
			int energy;
			in >> energy;
			SetEnergy(numofstructures,energy);
    		for (int i=1;i<=numofbases;i++) {
				in >> basepair;
				if (basepair>i) SetPair(i,basepair,numofstructures);
			}
		}
		nucs[numofbases+1]='\0'; // add terminating null character.
		return 0; //success (CCT format)
	} else if (count < 0 || count > MAX_CT_BASES) {
	    // Count would only be negative (other than -100) if the file is corrupt or it is not a true CT file. 
		// Report this error via the SetErrorDetails mechanism
		SetErrorDetails(count < 0 ? "Negative number of bases." : "Total number of bases exceeds maximum.");
		return 29; // invalid file format
	} else { // this is a ct file:
		
		//FOR NOW, DISABLE READING CT FILES WITH STACKING INFO: THIS SEEMS TO BE A PROBLEM WITH THE WEBSERVER BECAUSE THERE ARE MANY FORMATS IN CURRENT USE
		
		//first decide if it contains nucleotide stacking data at the far right
		//determine this based on the length of the line:
		//in.getline(temp,1000);
		//in.getline(temp,1000);

		//i = (int) strlen(temp);
		//if (i==35) {
		//	ct->stacking = true;
		//}
		//by default, ct->stacking is false


		//Allocate memory in the structure, the first thing read was the sequence length
		int position,linker=0; // linker is the number of intermolecular linkers found (i.e. linker+1 is the number of RNA molecules/strands)
		allocate(count);
		numofbases = count;
		long linenumber=1;//keep track of linenumber to report errors

		for (int numofstructures = 1;!in.eof();++(numofstructures))	{
			AddStructure();
			getline(in, sequenceLabel);
			SetCtLabel(sequenceLabel,numofstructures);
			//Also set the sequencelabel if this is structure #1
			if (numofstructures==1) SetSequenceLabel(sequenceLabel);
			
			for (count=1;count<=((numofbases));count++)	{
				++linenumber;//reading from a new line, so increment the line number
				in >> position;//base number

				//Check that the bases are incrementing correctly
				if (position!=count) {
					//There was a problem
					SetErrorDetails(sfmt("Invalid nucleobase index %i (expected %i) in structure %i at line %li.", position, count, numofstructures, linenumber));
					return 29;
				}

				in >> base;//read the base

				nucs[count]=base[0];

				numseq[count]=GetThermodynamicDataTable()->basetonum(base[0]);
				if (numseq[count]==-1) {
					//This means the bases was not recognized, and an error should be reported
					SetErrorDetails(sfmt("Invalid nucleobase '%c' in structure %i at line %li.", base[0], numofstructures, linenumber));
					return 29;
				}
				
				if (data->isLinker(numseq[count])&&numofstructures==1) {
      				intermolecular = true;
       				inter[linker++] = count;
				}
				in >> position;//read the previous connected nucleotide
				//Check the numbering
				if (position!=(count-1)&&position!=0) {
					//There was a problem
					SetErrorDetails(sfmt("Unexpected backbone connection between nucleotides %i and %i in structure %i at line %li.", count, position, numofstructures, linenumber));
					return 29;
				}
				in >> position;//read the next connected nucleotide
				//Check the numbering
				if (position!=(count+1)&&position!=0) {
					//There was a problem
					SetErrorDetails(sfmt("Unexpected backbone connection between nucleotides %i and %i in structure %i at line %li.", count, position, numofstructures, linenumber));
					return 29;
				}

				// basepair is the index that the current base is paired with, according to the CURRENT line in the file.
				// (basepair can be either before or after this one)
				in >> basepair;

				// Verify that the base is not paired with itself.
				if (basepair==count) {
					SetErrorDetails(sfmt("Base %i is paired with itself in structure %i at line %li.", count, numofstructures, linenumber));
					return 29;
				}

				// Get the index that THIS pair is ALREADY assigned to pair with.
				// This will be 0 unless the CURRENT nucleotide is the 3' base of a pair that was 
				// first mentioned on a PREVIOUS line in the file.
				int assigned = GetPair(count,numofstructures);

				// For base-pairing consistency one of the following must be true:
				//  1).  basepair == 0    AND assigned == 0         (This base is not paired to any others)
				//  2).  basepair > count AND assigned == 0         (This is the 5' base -- paired to one that comes later in the sequence)
				//  3).  basepair < count AND assigned == basepair  (This is the 3' base -- paired to one that came earlier in the sequence)
				//       Note that #1 is a specific case of #3, so if #3 is checked, #1 can be ignored.
				// As a corollary, here are some examples of inconsistencies:
				//  1)   basepair == 0      AND  assigned != 0         (A previous base said it was paired to this one, but this one says it is unpaired.)
				//  2)   basepair > count   AND  assigned != 0         (This base says it is paired to a later base, but the later base is already paired (which it shouldn't be until after THIS line is processed.)
				//  3)   basepair < count   AND  assigned != basepair  (A previous base said it was paired to this one, but this one says it is paired to a DIFFERENT base (not the previous one).

				if ((basepair<count&&assigned!=basepair) || (basepair>count&&assigned!=0)) {
					// For error reporting, we might have to swap some numbers. count is always THIS base. But either assigned or basepair could be 0. If so, swap them.
					if (assigned==0) assigned = basepair;
					SetErrorDetails(sfmt("Inconsistent base pairing information between nucleotides %i and %i in structure %i at line %li.", count, assigned, numofstructures, linenumber));
					return 29;
				}

				// One scenario that isn't covered above is that two bases are both paired to the same one (which comes AFTER the other two)
				// e.g.:  Pairs: 1->5  3->5  5->3   Steps: 1. SetPair(1,5);   2. Verify: GetPair(3)==0 and basepair>count; SetPair(3,5)  3. Verify: GetPair(5)==3==basepair.  (no errors detected so far)
				if (basepair>count) {
					assigned = GetPair(basepair,numofstructures); //here assigned is the base to which the TARGET base is assigned. 
					//assigned should be 0. If not, it means a previous base already listed it as being paired to it.
					if (assigned!=0&&assigned!=basepair) {
						SetErrorDetails(sfmt("Inconsistent base pairing information. Bases %i and %i are both paired to base %i in structure %i at line %li.", assigned, count, basepair, numofstructures, linenumber));
						return 29;
					}
				}

				if (basepair>count) SetPair(count,basepair,numofstructures);//set base pairing info
				in >> hnumber[count];//read historical numbering
				//if (stacking) in >> basepr[numofstructures][count+numofbases];
			}

			++linenumber;//read another line
			in >> count; //start on next structure and see whether the end of file is reached
		}
	}

	nucs[numofbases+1]='\0'; // add terminating null character.
	return 0; //success
}

//# same as openseq except that the return value is 0 on success or an error code corresponding to RNA::GetErrorMessage(int)
int structure::openseqx (const char *seqfile) {
	if (!IsAlphabetLoaded()) return 30; //error -- parameters have not been read.
	string line;
	string sequence;
	size_t start = string::npos;
	long linenumber = 0;
	enum SeqFileType { TEXT, SEQ, FASTA } fileType = TEXT;  // TEXT means plain text (AUGC...) with no comments or other delimiters.

	std::istream input(std::cin.rdbuf()); // if the filename is "-", use STDIN for input.
	std::ifstream file_in;
	if (!isStdIoFile(seqfile)) { 
		if (!fileExists(seqfile)) return 1; // file not found.
		//Open the file for reading:
		file_in.open(seqfile);
		//Make sure the file opened.  If not, return an error indicator.
		if(!file_in.is_open()) return 2;
		input.rdbuf(file_in.rdbuf());
	}

	//Now identify the file type.
    //  Usually, a starting > means FASTA while a starting ; means .seq 
	//  However, the FASTA/Pearson format allows for comment lines starting with ; See: https://en.wikipedia.org/wiki/FASTA_format

	//scan for first relevant line. Ignore lines that are empty, contain only whitespace, or start with a semi-colon (;).
	while (!getline(input,line).fail()) {
		++linenumber;
		if (line.length() == 0 || line[0] == '\r') continue; // line is empty
		if (line[0] == ';') {
			fileType = SEQ; // assume it is a SEQ, but this will be change to FASTA if the first character of the sequence label is '>'
			continue;
		}
		if (string::npos != (start = line.find_first_not_of(" \t\r"))) // if the line has any non-whitespace character, exit the loop
			break;
	}

	// Now past the comments. See if a valid line was found.
	if (start == string::npos) {
		//The end of the file was reached before finding the end of the comments, return with error
		SetErrorDetails("The file did not contain any nucleotides.");
		return 28; //Error reading sequence
	}

	if (line[start] == '>') { fileType = FASTA; start++; } // FASTA files must start with a '>' (following any ;-comments)

	if (start != 0) line.erase(0, start); // remove all characters up to, but not including start. (they were either whitespace or '>')
	
	bool reuseLine = false; // used for plain-text sequences

	//Set the current line as the sequence label (if FASTA or SEQ)
	if (fileType == SEQ || fileType == FASTA)
		SetSequenceLabel(line);
	else {
		// The filetype is plain text. 
		// *** Uncomment the following to disable parsing of plain-text sequences: ***
		//	SetErrorDetails("The format of the sequence file was invalid. It did not appear to be a SEQ or FASTA file.");
		//	return 29;
		SetSequenceLabel(getFileName(seqfile, true)); // Use the file name as a sequence label.
		reuseLine = true; // I'd prefer to use seekg, but there seems to be a bug in it on Windows.
	}

	bool foundEndChar=false; // indicates that a character was found that indiates the termination of a sequence. (for SEQ this is a '1')

	//now read each line and process it:
	while ((reuseLine||getline(input,line))&&!foundEndChar) {
		if (reuseLine) 
			reuseLine = false;
		else
			++linenumber;

		//Process each character in the current line
		for (int i=0;i<line.length();++i) {
			if (line[i] < 33) continue; // skip whitepace.

			if ((fileType==SEQ && line[i]=='1') || (fileType==FASTA && line[i]=='>')) {
				//Done reading sequence
				foundEndChar=true;
				break;
			} 

			// If the character is valid, add it. Otherwise return with an error.
			if (-1==data->basetonum(line[i])) {
				//This nucleotide wasn't recognized
				// report this error via the lastErrorDetails mechanism
				SetErrorDetails(sfmt("Invalid nucleobase '%c' at line %li column %i.", line[i], linenumber, i+1));
				return 28; //Error reading sequence
			} else {
				//Put the nucleotide in the sequence
				sequence+=line[i];
			}
		}//end of iteration over characters in the line
	} //end of while getLine loop

	if (fileType == SEQ && !foundEndChar) { 
		//The end of the file was reached before finding the 1 that indicates the sequence end
		// report this error via the SetErrorDetails mechanism
		SetErrorDetails("The file was missing the required '1' (one) that indicates the end of a sequence in a SEQ file.");
		return 29; //Invalid file format
	}
	
	if (sequence.size()==0) {
		SetErrorDetails("The file did not contain any nucleotides.");
		return 28; //Error reading sequence.
	}

	//The sequence string is populated with the sequence, now allocate the needed space in structure, and record the sequence
	SetSequence(sequence);
	return 0; //success (no error code)
}

//Add a new structure:
void structure::AddStructure() {


	arrayofstructures.push_back(singlestructure(numofbases));

	//If this is the first structure, go ahead and copy the sequence label to the structure label:
	if(arrayofstructures.size()==1) {
		arrayofstructures[0].ctlabel = sequencelabel;

	}

	


}

//Remove all pairs from a structure:
void structure::CleanStructure(int structurenumber) {
	int i;

	for (i=1;i<=numofbases;++i) arrayofstructures[structurenumber-1].basepr[i]=0;



}





//Remove the last structure:
void structure::RemoveLastStructure() {

	arrayofstructures.pop_back();

}

//This function sizes the vectors ctlable and energy to the right size for the maximum expected number of structures.
//This function is called multiple times if the number of structures is going the exceed the allocatedstructures.
void structure::allocatestructure(int structures) {
	
	
	//Set the maximum number of structures, to aid the vector in finding memory once:
	arrayofstructures.reserve(structures+1);

	

}



//Get the pointer to the underlying thermodynamic data
datatable *structure::GetThermodynamicDataTable() {
	return data;
}

//Set the pointer to the underlying thermodynamic data 
void structure::SetThermodynamicDataTable(datatable *DataTablePointer) {

	data=DataTablePointer;
}

// Returns true if the data property has been set to a valid datatable and its alphabet has been successfully read.
bool structure::IsAlphabetLoaded() {
	return data!=NULL&&data->loadedAlphabet;
}

//sort the structures from lowest to highest free energy
void structure::sort() {
	
	
	//Use sort, and refer to the energy sort function: 
		
	std::sort(arrayofstructures.begin(),arrayofstructures.end());

	

}


//The destructor:
structure::~structure()
{
	int i;



	
	if (allocated) {
		delete[] numseq;

		delete[] hnumber;
		delete[] nucs;
	}
	if (templated) {
		for (i=0;i<=numofbases;i++) {
    		delete[] tem[i];
   		}

   		delete[] tem;
	}
	DeleteSHAPE();

	if ( experimentalPairBonusExists ) {
		delete[] EX;
	}
	if (constant!=NULL) {
		//delete the equilibrium constants
		for (i=0;i<=numofbases;i++) {
    		delete[] constant[i];
   		}

   		delete[] constant;

	}
	
}



void structure::allocate(int size)

{
	
	numofbases = size;

	//Try not to include the following, unless found necessary in the future:
	//if (allocated) {
	//		//Delete memory if already allocated
	//	delete[] numseq;
	//	delete[] hnumber;
	//	delete[] nucs;

	//}
	numseq = new short int [2*size+1];
	hnumber = new short int [size+1];
	nucs = new char [size+2];
   
	allocated = true;



}



// Returns the textual name for the type of restraints (e.g. SHAPE)
// Used in warning messages etc.
const char* restraintTypeName(RestraintType modifier) {
	switch(modifier) {
		case RESTRAINT_SHAPE:return "SHAPE";
		case RESTRAINT_SHAPE_DIFF:return "diffSHAPE";
		case RESTRAINT_SHAPE_AC: return "SHAPE_AC";
		case RESTRAINT_SHAPE_GU: return "SHAPE_GU";
		case RESTRAINT_DMS: return "DMS";
		case RESTRAINT_CMCT:return "CMCT"; 
		// For unknown types, just return the simple description "restraint"
		default: return "restraint"; // e.g. "Warning: Error in >>restraint<< file ..."
	}
}

double structure::Gammadist(const double data, const double shape, const double loc, const double scale){
	return (1/scale)*pow((data - loc)*(1/scale), (shape - 1))*exp(-(1/scale)*(data - loc))/tgamma(shape);
}


/*
double structure::Potential(const double data, const std::vector< std::vector<double> > &params, const double kT){
	// params[0] is for paired, params[0] for unpaired...params[][j], j=0,1,2 for shape, loc scale of component 1
	// j=3,4,5 for shape, loc, scale of component 2 and j=6,7 for weights of components 1 and 2 respectively.
	double pairedprob = params[0][6]*Gammadist(data, params[0][0], params[0][1], params[0][2]) + 
	                   params[0][7]*Gammadist(data, params[0][3], params[0][4], params[0][5]); 
	double unpairedprob = params[1][6]*Gammadist(data, params[1][0], params[1][1], params[1][2]) + 
	                     params[1][7]*Gammadist(data, params[1][3], params[1][4], params[1][5]);
	return -kT*log(pairedprob/unpairedprob);
}
*/

double structure::Potential(const double data, const std::vector< std::vector<double> > &params, const double kT, const int ntcode) {
    // params[i] is for paired, params[i+1] for unpaired...
    // params[][j], j=0,1,2 for shape, loc scale of component 1
	// j=3,4,5 for shape, loc, scale of component 2 
    // j=6,7 for weights of components 1 and 2 respectively.
    
    int idx;
    switch(ntcode) {
        case 1: idx=0;  //A (or for default case without nuc-specific potentials)
                break;
        case 2: idx=2;  //C
                break;
        case 3: idx=4;  //G
                break;
        case 4: idx=6;  //U
                break;
        case 0: return 0.0; // nuc is X,N
        case 5: return 0.0; // nuc is I 
    }

	double pairedprob = params[idx][6]*Gammadist(data, params[idx][0], params[idx][1], params[idx][2]) + 
	                    params[idx][7]*Gammadist(data, params[idx][3], params[idx][4], params[idx][5]); 
	double unpairedprob = params[idx+1][6]*Gammadist(data, params[idx+1][0], params[idx+1][1], params[idx+1][2]) + 
	                      params[idx+1][7]*Gammadist(data, params[idx+1][3], params[idx+1][4], params[idx+1][5]);
    
    return -kT*log(pairedprob/unpairedprob);        
}



void structure::ReadProbabilisticPotentialParams() {
	string filedir(getDataPath()); // getDataPath will return DATAPATH if it exists or "." otherwise.
	filedir += "/dists/";

	string line;
	int start, end, i, j;
	int nparams = 8;
	// Initialize parameters
	std::vector<double> shapecol1;
	std::vector<double> shapecol2;
	for(i=0; i<nparams; i++){
		shapecol1.push_back(0.0);
		shapecol2.push_back(0.0);
	}
	SHAPE_params.push_back(shapecol1);
	SHAPE_params.push_back(shapecol2);
    
    std::vector<double> dmscol1;
	std::vector<double> dmscol2;
	for(i=0; i<nparams; i++){
		dmscol1.push_back(0.0);
		dmscol2.push_back(0.0);
	}
	DMS_params.push_back(dmscol1);
	DMS_params.push_back(dmscol2);

    for(i=0; i<8; i++) {
        std::vector<double> dmscol;
	    for(j=0; j<nparams; j++){
		    dmscol.push_back(0.0);
	    }
	    DMS_paramsnt.push_back(dmscol);
    }
    
	std::vector<double> cmctcol1;
	std::vector<double> cmctcol2;
	for(i=0; i<nparams; i++){
		cmctcol1.push_back(0.0);
		cmctcol2.push_back(0.0);
	}
	CMCT_params.push_back(cmctcol1);
	CMCT_params.push_back(cmctcol2);

	// Start reading distribution parameters....note that this code is ugly, can be collapsed into a single parameter array
	// SHAPE
	string tmp = filedir + "SHAPEdist.txt";
	char *filename = (char*)tmp.c_str();
	ifstream SHAPEfile;
	SHAPEfile.open(filename);
	if (SHAPEfile.good()){
		getline(SHAPEfile, line);
		for(i=0;i<2;i++) {
			start = 0;
			end = 0;
			getline(SHAPEfile, line);
			for(j=0;j<nparams;j++) {
				end = line.find(" ", start);
				SHAPE_params[i][j] = atof(line.substr(start, end).c_str());
				start = end + 1;
			}
		}
		SHAPEfile.close();
	}
	else {
		cout << "Cannot open file " + tmp + "\n";
	}

	// DMS
	tmp = filedir + "DMSdist.txt";
	filename = (char*)tmp.c_str();
	ifstream DMSfile;
	DMSfile.open(filename);
	if (DMSfile.good()){
		getline(DMSfile, line); // pop off the header
		for(i=0;i<2;i++) {
			start = 0;
			end = 0;
			getline(DMSfile, line);
			for(j=0;j<nparams;j++) {
				end = line.find(" ", start);
				DMS_params[i][j] = atof(line.substr(start, end).c_str());
				start = end + 1;
			}
		}
		DMSfile.close();
	}
	else {
		cout << "Cannot open file " + tmp + "\n";
	}
    
    
    // DMS NT
	tmp = filedir + "DMSdist_nt.txt";
	filename = (char*)tmp.c_str();
	ifstream DMSNTfile;
	DMSNTfile.open(filename);
	if (DMSNTfile.good()){
		getline(DMSNTfile, line);
		for(i=0;i<8;i++) {
			start = 0;
			end = 0;
			getline(DMSNTfile, line);
			for(j=0;j<nparams;j++) {
				end = line.find(" ", start);
				DMS_paramsnt[i][j] = atof(line.substr(start, end).c_str());
				start = end + 1;
			}
		}
		DMSNTfile.close();
	}
	else {
		cout << "Cannot open file " + tmp + "\n";
	}
    
	// CMCT
	tmp = filedir + "CMCTdist.txt";
	filename = (char*)tmp.c_str();
	ifstream CMCTfile;
	CMCTfile.open(filename);
	if (CMCTfile.good()){
		getline(CMCTfile, line);
		for(i=0;i<2;i++) {
			start = 0;
			end = 0;
			getline(CMCTfile, line);
			for(j=0;j<nparams;j++) {
				end = line.find(" ", start);
				CMCT_params[i][j] = atof(line.substr(start, end).c_str());
				start = end + 1;
			}
		}
		CMCTfile.close();
	}
	else {
		cout << "Cannot open file " + tmp + "\n";
	}
}

// This function calculates the pseudoenergy for a given reactivity data. It changes the calculation
// depending on the modifier specified, giving either the log-likelihood-ratio of the unpaired/paired
// probabilities given a reactivity distribution per modifier, or the "classic" Deigan et al. bonus
// term when no modifier or an unrecognized modifier is provided. 
double structure::CalculatePseudoEnergy(const double data, const RestraintType modifier, const double slope, const double intercept, const int ntcode, const bool use_params_if_available ) {
	

    std::vector< std::vector<double> > *params; // use a pointer to avoid making a copy
    int ntcode2;
    double data2 = data;
    double kT = 5.904976983149999;  // kT at 24 C in 10ths of kcal/mol. Value used by Cordero et al.
    

    // slope==0 && intercept==0 & !use_params_if_available is triggered by SHAPEss calls to CalculatePseudoEnergy
    if( data <= -500 || (slope == 0 && intercept == 0 && !use_params_if_available )) {
		return 0;
    }
    
    switch(modifier) {
		case RESTRAINT_SHAPE_AC:
		case RESTRAINT_SHAPE_GU:
			// This is only applied if SHAPE_AC or SHAPE_GU is specified
			// By default, use the Deigan "default" method when the modifier is "SHAPE"
			params = &SHAPE_params;
            ntcode2 = 1;
			break;
		case RESTRAINT_DMS:
			params = &DMS_params;
            ntcode2 = 1;
			break;
        case RESTRAINT_DMSNT:
            params = &DMS_paramsnt;
            ntcode2 = ntcode;
            kT = 3.0816567; // 0.5*kT at 37 C in 10ths of kcal/mol; 0.5 adjust is for double counting of internal pairs
            break;
		case RESTRAINT_CMCT:
			params = &CMCT_params; 
            ntcode2 = 1;
			break;
		default: 
			// Calculate energies for SHAPE and diffSHAPE.
			// For negative values of data, just return the intercept (baseline).
			// Otherwise, return the calculated pseudo-energy.
			return data>0 ? log(data+1.0)*slope+intercept : intercept;
	}
	// If we get here, we are using DMS, DMSNT, CMCT, SHAPE_AC or SHAPE_GU.
    
    if( params->empty() ) 
        return 0;
    
    // SHAPE default returns `intercept` when data<0
    // For Gamma potentials, data<=0 is undefined, and is poorly behaved for data->0
    // For DMSNT params, 0.01 is well-behaved valuea, so evalute here as `intercept`
    // Behavior for earlier (Cordero et al 2012) parameters has not been benchmarked, 
    // so preserve `return 0` for these restraints. However, I think it would make sense
    // to consider also evaluting these potentials at data~0.01 to return the equivalent 
    // of an intercept.
    if (modifier == RESTRAINT_DMSNT && data < 0.005) {
        data2 = 0.005;
    } else if (data < 0) {
        return 0;
    }

	const double val  = Potential(data2, *params, kT, ntcode2);
    


 	// RMW: why not return `intercept` instead of 0?  If this function returns `intercept` for SHAPE data <= 0.
	return std::isnan(val) ? 0 : val;

}

// Allocate the SHAPE and SHAPEss arrays if they haven't already been created.
// Also initializes SHAPEss_region, a triangular 2-D array that stores ss SHAPE energies for loops.
void structure::AllocateSHAPE() {
	if (shaped) return; // already created.
	//initializes an array to hold SHAPE data for double stranded segments
	SHAPE = new double [2*numofbases+1];
	//initializes an array to hold SHAPE data for single stranded segments
	SHAPEss = new double [2*numofbases+1];
	shaped = true;
	// Initialize all points to zero.
	for (int position=0; position<=2*numofbases; position++) {
		SHAPE[position] = 0;
		SHAPEss[position] = 0;
	}
	// initialize triangular 2-d array that stores ss SHAPE energies for loops. 1st index is ending location, 2nd index is starting location
	SHAPEss_region = new short int*[numofbases+1];
	for (int i=1; i<=numofbases; i++) SHAPEss_region[i] = new short int[i];
}

// Deletes the SHAPE, SHAPEss, and SHAPEss_region arrays if the have been allocated.
// Also sets the `shaped` field to false and sets the SHAPE, SHAPEss, and SHAPEss_region
// pointers to NULL.
void structure::DeleteSHAPE() {
	if (!shaped) return; // already created.
	if (SHAPE!=NULL)    delete[] SHAPE;
	if (SHAPEss!=NULL)  delete[] SHAPEss;
	if (SHAPEss_region!=NULL) {
		for (int i=1; i<=numofbases; i++)
			delete[] SHAPEss_region[i];
		delete[] SHAPEss_region;
	}
	shaped=false;
	SHAPE=SHAPEss=NULL;
	SHAPEss_region=NULL;
}

// This function reads a SHAPE reactivity datafile and saves the data for a linear penalty.
// calculatePseudoEnergies (default true) indicate whether these data are being read for folding.  
//  (false means the raw values need to be stored.)
int structure::ReadSHAPE(const char *filename, RestraintType modifier, bool calculatePseudoEnergies) {
	int position;
	double data;

	//Only read the dists if they need to be read.  These are not used for SHAPE or diffSHAPE.
	if( !distsread && !(modifier==RESTRAINT_SHAPE||modifier==RESTRAINT_SHAPE_DIFF) ){
		ReadProbabilisticPotentialParams();
		distsread=true;
	}

	// Allocate the SHAPE and SHAPEss arrays if they haven't already been created.
	// Also initializes SHAPEss_region, a triangular 2-D array that stores ss SHAPE energies for loops.
	AllocateSHAPE();

	// temporary array to hold reactivity data for double-stranded segments
	double *SHAPEnew = new double [2*numofbases+1];
	// temporary array to hold reactivity data for single stranded segments
	double *SHAPEssnew = new double [2*numofbases+1];
	// useful for bootstrapping -- keep count of how many times the data for a given position is specified.
	int *num_data_points =  new int [numofbases+1];

	// the array pointers above will be automatically deleted if return is called before the end of the function.
	auto_delete<double,true> d1(SHAPEnew), d2(SHAPEssnew); auto_delete<int,true> d3(num_data_points);

	for (position=0; position <= 2*numofbases; position++) {
		SHAPEnew[position] = 0.0;
		SHAPEssnew[position] = 0.0;
	}
	for (position=0; position <= numofbases; position++) {
		num_data_points[position] = 0;
	}

	if (!fileExists(filename)) return 201; //file not found
	ifstream in(filename);
	if (!in.good()) return 202; // error opening file


	vector<int> invalidPositions;
	bool hasRepeats = false; //track whether any values have been repeated.
	while (!(in >> position >> data).fail()) { // this will ensure the last data point is read, even if there is no newline at the end of the file. (EOF can be set after a valid datapoint is read if there is no whitespace following it.)
		//read and parse all data
			//required format is rows with sequence position followed by reactivity
		if (position<1||position>numofbases) {
			invalidPositions.push_back(position);
			continue;
		}
			if (calculatePseudoEnergies) { 
				//store calculated pseudo-energy values.
				SHAPEnew[position]   += CalculatePseudoEnergy(data, modifier, SHAPEslope, SHAPEintercept, numseq[position], true /*use_params_if_available*/ );
				SHAPEssnew[position] += CalculatePseudoEnergy(data, modifier, SHAPEslope_ss, SHAPEintercept_ss, numseq[position],  false /*use_params_if_available*/ );
                
                //cout << position << " " << nucs[position] << " " << data << " " << SHAPEnew[position]/10 << "\n";

			} else { 
				//store the raw SHAPE values
				SHAPE[position] = data;
				SHAPEss[position] = data;
			}
			if (num_data_points[ position ]) hasRepeats = true;
			num_data_points[ position ]++;
	}
	in.close();
    

	if (!invalidPositions.empty())
		cwarn() << "Warning: Invalid nucleobase positions in " << restraintTypeName(modifier) << " file " << filename << ": " << invalidPositions << ". (Sequence length is " << numofbases << ".)" << endl;

	if (calculatePseudoEnergies) {
		for (position=1;position<=numofbases;position++) {
			if(num_data_points[position] > 0){
				// SumShapeRepeats is true by default, but can be turned off by setting the environment variable AVG_SHAPE_REPEATS to 1.
				if (structure::SumShapeRepeats) { 
					// SHAPE/DMS/CMCT files can be ‘resampled with replacement’ (bootstrapped) to
					//  mock simulate having extra or no data at some positions — we need
					//  to *accumulate* the pseudo energies there. 
					SHAPE[position]+=SHAPEnew[position];
					SHAPEss[position]+=SHAPEssnew[position];
				} else {
					SHAPE[position]+=SHAPEnew[position]/num_data_points[position];
					SHAPEss[position]+=SHAPEssnew[position]/num_data_points[position];
				}
			}
		}

		// SHAPE values in ‘second set’ must exactly match SHAPE values in first set.
		for (position=1;position<=numofbases;position++) {
			SHAPE[position+numofbases]=SHAPE[position];
			SHAPEss[position+numofbases]=SHAPEss[position];
		}
	}

	if (hasRepeats && ShowWarnings!=0 && SumShapeRepeats) {
		// Show a warning about duplicate SHAPE values
		// unless warnings are suppressed or the user has expressly 
		// permitted multiple values by setting the AVG_SHAPE_REPEATS environment variable.
		ostream& warn = cwarn(); // (cwarn() gives the user some control over how warnings are shown.)
		warn << "Warning: The following nucleobase positions were repeated in " << restraintTypeName(modifier) << " file " << filename << ":";
		for (position=1;position<=numofbases;position++)
			if (num_data_points[position])
				warn << " " << position;
		warn << endl << "(This may be OK if you are bootstrapping -- resampling with replacement.)" << endl;
	}

	//fills SHAPEss_region, a 2-d array with pseudo energy terms for loops from i-j using ss SHAPE parameters or offsets
	FillSHAPEssRegions();

	// (deletion of temporary array is handled by the auto_delete objects defined above)
	return 0;
}

// Return an exact copy of SHAPE data.
// The return value is a pointer to a dynamically allocated array of doubles.
// If includeSHAPEss is true, the array will contain both SHAPE and SHAPEss values and will
// have a size of 4*numbases+2 (with SHAPE data from arr[0] to arr[2N] and SHAPEss from arr[2N+1] to arr[4N+1])
// If includeSHAPEss is false, the array will only contain SHAPE data and will have a size of 
// 2*numofbases+1.
double* structure::CopySHAPE(const bool includeSHAPEss) {
	if (!shaped) return NULL;
	double* arr = new double[(includeSHAPEss?2:1)*(2*numofbases+1)]; // 2N+1 or 4N+2
	// Fill arr[0] to arr[2N] with SHAPE   (note that arr[0] is unused since nucleotide indexing is 1-based)
	for (int i=0; i<=2*numofbases; i++) 
		arr[i]=SHAPE[i];
	if (includeSHAPEss)
		// Fill arr[2N+1] to arr[4N+1] with SHAPEss  (note that arr[2N+1] is unused since nucleotide indexing is 1-based)
		for (int i=0; i<=2*numofbases; i++) 
			arr[2*numofbases+1+i]=SHAPEss[i];
	return arr;
}

// Load data from an array of doubles into the SHAPE array.
// This function copies the data (it does NOT take ownership of the passed-in array).
// If includeSHAPEss is true, the array must contain both SHAPE and SHAPEss values and must
// have a size of at least 4*numbases+2 (with SHAPE data from arr[0] to arr[2N] and SHAPEss from arr[2N+1] to arr[4N+1])
// If includeSHAPEss is false, the array need only contain SHAPE data and must have a size of 
// at least 2*numofbases+1.
// If the passed-in array is NULL, this structure's SHAPE data will be deleted.
void structure::LoadSHAPE(const double* shapeArray, const bool includeSHAPEss) {
	if (shapeArray==NULL)
		DeleteSHAPE();
	else {
		AllocateSHAPE(); // create arrays if not already created.
		// Fill arr[0] to arr[2N] with SHAPE   (note that arr[0] is unused since nucleotide indexing is 1-based)
		for (int i=0; i<=2*numofbases; i++) 
			SHAPE[i]=shapeArray[i];
		if (includeSHAPEss)
			// Fill arr[2N+1] to arr[4N+1] with SHAPEss  (note that arr[2N+1] is unused since nucleotide indexing is 1-based)
			for (int i=0; i<=2*numofbases; i++) 
				SHAPEss[i]=shapeArray[2*numofbases+1+i];
	}
}

// Fills SHAPEss_region, a 2-d array with pseudo energy terms for loops from i-j using
//   single-stranded (ss) SHAPE parameters or offsets.
// The SHAPEss and SHAPEss_region arrays must already be allocated (by calling AllocateSHAPE)
//   and SHAPEss must already contain ss SHAPE pseudo-energies.
// This is called from e.g. ReadSHAPE and ReadOffset.
void structure::FillSHAPEssRegions() {
	//fills 2-d array with pseudo energy terms for loops from i-j using ss SHAPE parameters or offsets
	//If SHAPE was previously read, this is a redo of the action
	for (int j = 2; j <= numofbases; j++) {
		SHAPEss_region[j][j - 1] = (short int)(SHAPEss[j] + SHAPEss[j-1]); //sets energy for "zero sized loop".  Acts as starting value to add onto below
		for (int i = j - 2; i >= 1; i--) {
			SHAPEss_region[j][i] = SHAPEss_region[j][i + 1] + (short int)(SHAPEss[i]); //adds energy for additional loop element
		}
	}
}

//this function returns a psuedo energy term for a single stranded nucleotide based on SHAPE data
short int structure::SHAPEss_give_value(int index) {
	if (shaped) {
		if (index > numofbases) return (short int)SHAPEss[index - numofbases];
		else return (short int)SHAPEss[index];
	}
	else return 0;
}

//this is the function that will return a pseudo energy term for hairpin loops based off of SHAPE data
int structure::SHAPEss_calc(int index_i, int index_j) {
	if (shaped) {
		//accounts for the SHAPEss_region array being only NxN, not 2Nx2N
		if (index_i > numofbases) index_i -= numofbases;
		if (index_j > numofbases) index_j -= numofbases;
		if (index_i > index_j) {
			int temp_index = index_i;
			index_i = index_j;
			index_j = temp_index;
		}
		return SHAPEss_region[index_j][index_i];
	} else return 0;  //if no shaped data is being used, return zero
}


// This function reads a SHAPE reactivity datafile and parse the data into single-stranded amd chemical modification constraints.
// This function is largely deprecated by the pseudo-free energy approach.  It is still available for experimentation.
// Returns 0 on success or an error code compatible with RNA::GetErrorMessage (for example 1=file not found, 2=could not open file)
int structure::ReadSHAPE(const char *filename, float SingleStrandThreshold, float ModificationThreshold) {
	int position;
	float data;

	if (!fileExists(filename)) return 201;  // file not found
	ifstream in(filename);
	if (!in.good()) return 202; // error opening file.

	vector<int> invalidPositions;
	while (!(in >> position >> data).fail()) {
		//read and parse all data
			//required format is rows with sequence position followed by reactivity
		if (position<1||position>numofbases) {
			invalidPositions.push_back(position);
			continue;
		}
		if (data>=SingleStrandThreshold) {
			AddSingle(position);
		}
		else if (data>=ModificationThreshold) {
			AddModified(position);
		}
	}
	in.close();
	if (!invalidPositions.empty())
		cwarn() << "Warning: Invalid nucleobase positions in SHAPE file " << filename << ": " << invalidPositions << ". (Sequence length is " << numofbases << ".)" << endl;
	return 0;
}

// This function reads an experimental pair bonus file, similar to SHAPE, but just straightforward
// application as kcal bonuses.  As with SHAPE, bonus is applied at 0x, 1x, and 2x for
//  single stranded, edge base pairs, and internal base pairs.
// Returns 0 on success or an error code compatible with RNA::GetErrorMessage (for example 1=file not found, 2=could not open file)
int structure::ReadExperimentalPairBonus(const char *filename, double const experimentalOffset, double const experimentalScaling ) {
	//int position;
	//float data;

	//initializing 2-d array that stores bonuses for base pairs.
	// actually this is already done in RNA.cpp. Thanks Pablo!
	//delete[] EX;
	EX = new double *[ 2*numofbases+1 ];
	for (int i = 0; i < 2*numofbases+1; i++) EX[i] = new double [ 2*numofbases+1 ];
	for (int i = 0; i < 2*numofbases+1; i++)
	  for (int j = 0; j < 2*numofbases+1; j++)
	    EX[i][j] = 0.0;

	for (int i = 1; i <= 2*numofbases; i++)
	  for (int j = 1; j <= 2*numofbases; j++)
	    EX[i][j] = experimentalOffset * conversionfactor;

	int i( 1 ), j( 1 ), count( 0 );
	double val;
    char header;
    
    int numofbases2 = numofbases*numofbases;


    // check that file exists 
    if (is_blank(filename) || !fileExists(filename)) return 201; // file doesn't exist
    
    ifstream in(filename);
    if (!in.good()) return 202; //error opening file
    
    // determine file format
    in.get(header);
    if (header == ';') {
        // format is columnar; expecting three columns: i(int) j(int) val(float)

        in.ignore(1000,'\n'); // pop off the complete header line
        
        while (count < numofbases2) {
            
            in >> i;
            in >> j;
            in >> val;
            
            if (!in.good()) break; // either end of file or bad value couldn't be coerced into i/j/val
            
            EX[i           ][j           ] += val * PFPRECISION( conversionfactor ) * experimentalScaling;
	        EX[i+numofbases][j           ] = EX[i][j];
	        EX[i           ][j+numofbases] = EX[i][j];
	        EX[i+numofbases][j+numofbases] = EX[i][j];
            
            EX[j           ][i           ] += val * PFPRECISION( conversionfactor ) * experimentalScaling;
	        EX[j+numofbases][i           ] = EX[j][i];
	        EX[j           ][i+numofbases] = EX[j][i];
	        EX[j+numofbases][i+numofbases] = EX[j][i];
       
            count++;
        }
        

        if (in.eof()) {
            cout << count << " columnar pairing restraints read...";
        } else {
            SetErrorDetails(sfmt("Experimental bonus file '%s' intrepreted as columnar format contains improper value or is incorrectly formatted", filename));
            return 203; // see RNA::GetErrormessage
        }


    // no header, so matrix format
    } else {
        in.unget(); //rewind file
    
	    while (!in.eof() && j <= numofbases) {
	        //read and parse all data
	        //required format is bonuses in square matrix
	        in >> val;

	        EX[i           ][j           ] += val * PFPRECISION( conversionfactor ) * experimentalScaling;
	        EX[i+numofbases][j           ] = EX[i][j];
	        EX[i           ][j+numofbases] = EX[i][j];
	        EX[i+numofbases][j+numofbases] = EX[i][j];
	        count++;

	        i++;
	        if ( i > numofbases ){
	            i = 1;
	            j++;
	        }
	    }
	
        if ( count != numofbases2 ){
		    SetErrorDetails(sfmt("Experimental bonus file '%s' intrepreted as matrix format but did not have expected number of values\nFound %i but expected %i.\nIf columnar format, first line needs to start with ';'", filename, count, numofbases*numofbases));
	        return 203; // see RNA::GetErrormessage

	    }
	}
 

    in.close();
	experimentalPairBonusExists = true;
	return 0;
}




const string& structure::GetErrorDetails() {
	return lastErrorDetails;
}
void structure::SetErrorDetails(const string& details) {
	lastErrorDetails = details;
}

// returns 0, 1 or 2 based on the value of a c-string (ON=1, OFF=0, ERR=2)
int parse_OnOffErrFlag(const char*const cstr) {
	string warn=as_str(cstr); toUpper(warn);
	return (warn=="OFF"||warn=="0")?0:(warn=="ERR"||warn=="2")?2:1;
}
// Set default behavior regarding warnings. 1=ON (stdout) 1=ERR (stderr) 0=OFF
int structure::ShowWarnings(parse_OnOffErrFlag(getenv("RNA_WARNINGS")));

// affects how multiple datapoints for the same nucleboase are handled in SHAPE data.
// the default is to add them for the "resample with replacement" technique 
// but setting the environment variable AVG_SHAPE_REPEATS to 1 causes
// multiple values to be averaged instead.
bool structure::SumShapeRepeats(is_blank(getenv("AVG_SHAPE_REPEATS")));


// Return an appropriate output stream for warnings --
// one of cout, cerr or nullstream, depending on the value of ShowWarnings.
// TODO: Later this can be changed to allow a third option -- a stringstream to capture warnings for library clients (e.g. Java)
ostream& structure::cwarn() {
	return ShowWarnings==0?NullStream::Default:ShowWarnings==2?cerr:cout;
}

// This is the default implementation of a function to return a comment/label like 
// "ENERGY = ..." which is inserted before the existing structure label when writing 
// CT or dot-bracket files.
// If GetEnergy() is zero, this function returns "" to indicate that no comment
// should be inserted.
string CTComments::EnergyCommentProvider::getComment(const structure* ct, const int structurenumber) {
	const int energy = ct->GetEnergy(structurenumber);
	if (energy==0) return "";
	stringstream comment("ENERGY = "); comment.seekp(0, std::ios::end); // move to the end of the string so we append the rest.
	comment << std::fixed << std::setprecision(conversionprecision) << (energy/(float)conversionfactor);
	return comment.str();
}
// Create memory storage location for static members of CTComments.
CTComments::EnergyCommentProvider CTComments::Energy;
CTComments::NoCommentProvider     CTComments::None;


//read a ct file with sequence and structural information
#define linelength 20






//takes a nucleotide input and stores a numeric representation
//This is depracated because the alphabet needs to come from structure, which reads it from disk
/*void tonum(char *base,structure *ct,int count)	{
if (!strcmp(base,"A")) (ct->numseq[count] = 1);
else if(!strcmp(base,"B")) {
	(ct->numseq[count] = 1);

}
else if(!strcmp(base,"a")) {
	ct->numseq[count]=1;
	ct->AddSingle(count);

}
else if(!strcmp(base,"C")) (ct->numseq[count] = 2);
else if(!strcmp(base,"Z")) {
	(ct->numseq[count] = 2);

}
else if(!strcmp(base,"c")) {
	ct->numseq[count] = 2;
	ct->AddSingle(count);
}
else if(!strcmp(base,"G")) (ct->numseq[count] = 3);
else if(!strcmp(base,"H")) {
	(ct->numseq[count] = 3);

}
else if(!strcmp(base,"g")) {
	ct->numseq[count] = 3;
	ct->AddSingle(count);
}

else if(!strcmp(base,"U")||!strcmp(base,"T")) (ct->numseq[count] = 4);
else if(!strcmp(base,"V")||!strcmp(base,"W")) {
	(ct->numseq[count] = 4);

}
else if(!strcmp(base,"u")||!strcmp(base,"t")) {
	ct->numseq[count] = 4;
	ct->AddSingle(count);
}

else if(!strcmp(base,"I")) {
	ct->numseq[count] = 5;
	ct->intermolecular= true;
}

else (ct->numseq[count]=0);  //this is for others, like X
return;
}*/






//Depracated because of new extended alphabet
/*char *tobase (int i)

{  //function is a look up table to convert the base
	// 	integer represention to a familiar character

	if (i==1) return "A";
	 else if (i==2) return "C";
	 else if (i==3) return "G";
	 else if (i==4) return "U";
	 else if (i==0) return "X";
    else if (i==5) return "I";
	 else return "?";   //to track down mistakes

}*/








//#endif


//add the factor from SHAPE calculation
//This pseudo-energy was calculated when the file was loaded (see structure.cpp).
//The pseudo energy is applied twice for each nuc in interior pair and once for each nuc in terminal pair.
inline integersize SHAPEend(int i, structure *ct) {
	if (ct->shaped) return (integersize) ct->SHAPE[i];
	return 0;
}



//calculate the energy of stacked base pairs

//oriented by:
//5' i ip 3'
//   |  |
//   j jp

integersize erg1(int i,int j,int ip,int jp,structure *ct, datatable *data)
{

		integersize energy;

		 if ((i==(ct->GetSequenceLength()))||(j==((ct->GetSequenceLength())+1))) {
      		//this is not allowed because n and n+1 are not cavalently attached
			energy = INFINITE_ENERGY;
		}
		else {
      		energy = data->stack[(ct->numseq[i])][(ct->numseq[j])]
				[(ct->numseq[ip])][(ct->numseq[jp])]+data->eparam[1];
			//if (ct->shaped) {
			energy+=SHAPEend(i,ct);
			energy+=SHAPEend(j,ct);
			energy+=SHAPEend(ip,ct);
			energy+=SHAPEend(jp,ct);

		if (ct->experimentalPairBonusExists ) {
		  energy += 0.5 * ( ct->EX[i][j] + ct->EX[j][i] )
		    + 0.5 * ( ct->EX[ip][jp] + ct->EX[jp][ip] );
		}

				//}

		}
		return energy;
}


//calculate the energy of a bulge/internal loop
//where i is paired to j; ip is paired to jp; ip > i; j > jp
integersize erg2(int i,int j,int ip,int jp,structure *ct, datatable *data,
	char a, char b)
{

	integersize energy;
	int size,size1,size2,loginc, lopsid, energy2,count,k, energy_option;
	/* size,size1,size2 = size of a loop
		energy = energy calculated
		loginc = the value of a log used in large hairpin loops
	*/
   	if (((i<=(ct->GetSequenceLength()))&&(ip>(ct->GetSequenceLength())))||((
      	jp<=(ct->GetSequenceLength()))&&(j>(ct->GetSequenceLength())))) {
         //A loop cannot contain the ends of the sequence

         return INFINITE_ENERGY;
      }



      size1 = ip-i-1;
		size2 = j - jp - 1;



      if ((a>0)||(b>0)) {
      	if ((a&DUBLE)||(b&DUBLE)) return INFINITE_ENERGY;//the loop contains a nuc that
      		//should be double stranded
      	else if ((a&INTER)) {
         	//the loop is actually between two strands (ie: intermolecular)

             	if (size2>1) {//free energy is that of two terminal mismatches
               	//and the intermolecular initiation
                  energy = data->init + data->tstack[ct->numseq[i]][ct->numseq[j]]
                  	[ct->numseq[i+1]][ct->numseq[j-1]] +
                     data->tstack[ct->numseq[jp]][ct->numseq[ip]]
                  	[ct->numseq[jp+1]][ct->numseq[ip-1]];
               }
               else if (size2==1) {//find the best terminal mismatch and terminal
               	//stack free energies combination

                  energy = NO_COUNT(data->init) + NO_COUNT(data->tstack[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]]) +
                     	erg4_nc (jp,ip,ip-1,2,ct,data,false)+penalty_nc(jp,ip,ct,data);
				  
				  energy_option = 1;						 
                  
				  energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstack[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]]) +
                     	erg4_nc (i,j,i+1,1,ct,data,false)+penalty_nc(i,j,ct,data);

#ifdef COUNTING				 
				  if (energy2 < energy){
					energy_option=2;
					energy = energy2;
				  }		 
#else				  
				  energy = min (energy,energy2);
#endif //COUNTING

                  //if ((ct->numseq[i+1]!=5)&&(ct->numseq[ip-1]!=5)) {
                     //now consider if coaxial stacking is better:
                     energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstackcoax[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]])
                        + NO_COUNT(data->coaxstack[ct->numseq[jp+1]][ct->numseq[ip-1]][ct->numseq[j]][ct->numseq[i]])+penalty_nc(i,j,ct,data)+penalty_nc(jp,ip,ct,data);

#ifdef COUNTING				 
				  if (energy2 < energy){
					energy_option=3;
					energy = energy2;
				  }		 
#else				  
				  energy = min (energy,energy2);
#endif //COUNTING

                     energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstackcoax[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[j-1]][ct->numseq[ip-1]])
                        + NO_COUNT(data->coaxstack[ct->numseq[j-1]][ct->numseq[ip-1]][ct->numseq[j]][ct->numseq[i]])+penalty_nc(i,j,ct,data)+penalty_nc(jp,ip,ct,data);

#ifdef COUNTING				 
				  if (energy2 < energy){
					energy_option=4;
					energy = energy2;
				  }	

				  switch (energy_option){
						case 1:
							energy = data->init + data->tstack[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]]
								+ erg4(jp,ip,ip-1,2,ct,data,false)+penalty(jp,ip,ct,data);
							break;

						case 2:
							energy = data->init + data->tstack[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] 
								+ erg4(i,j,i+1,1,ct,data,false)+penalty(i,j,ct,data);
							break;

						case 3:
							energy = data->init + data->tstackcoax[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]]
								+ data->coaxstack[ct->numseq[jp+1]][ct->numseq[ip-1]][ct->numseq[j]][ct->numseq[i]]+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
							break;

						case 4:
							energy = data->init + data->tstackcoax[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[j-1]][ct->numseq[ip-1]] 
								+ data->coaxstack[ct->numseq[j-1]][ct->numseq[ip-1]][ct->numseq[j]][ct->numseq[i]]+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
							break;
					}		 
#else				  
				  energy = min (energy,energy2);
#endif //COUNTING

               }
               else if (size2==0) {//just have dangling ends or flush stacking
				energy = NO_COUNT(data->init) + erg4_nc (jp,ip,ip-1,2,ct,data,false) +
					erg4_nc (i,j,i+1,1,ct,data,false)+penalty_nc(i,j,ct,data)+penalty_nc(jp,ip,ct,data);
				
				energy2 = NO_COUNT(data->init) + NO_COUNT(data->coax[ct->numseq[ip]][ct->numseq[jp]][ct->numseq[j]][ct->numseq[i]])+penalty_nc(i,j,ct,data)+penalty_nc(jp,ip,ct,data);

#ifdef COUNTING                  
				if (energy< energy2){
					energy = data->init + erg4 (jp,ip,ip-1,2,ct,data,false) +
						erg4 (i,j,i+1,1,ct,data,false)+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
				}
				else{
					energy = data->init + data->coax[ct->numseq[ip]][ct->numseq[jp]][ct->numseq[j]][ct->numseq[i]]+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
				}
#else
                  energy = min(energy,energy2);
#endif //COUNTING				  
               }


         		return energy;
      	}
         else if (b&INTER) {
                  	//the loop is actually between two strands (ie: intermolecular)

             	if (size1>1) {//free energy is that of two terminal mismatches
               	//and the intermolecular initiation
                  energy = data->init + data->tstack[ct->numseq[i]][ct->numseq[j]]
                  	[ct->numseq[i+1]][ct->numseq[j-1]] +
                     data->tstack[ct->numseq[jp]][ct->numseq[ip]]
                  	[ct->numseq[jp+1]][ct->numseq[ip-1]];
               }
               else if (size1==1) {//find the best terminal mismatch and terminal
               	//stack free energies combination
				  energy_option=1;
 
				energy = NO_COUNT(data->init) + NO_COUNT(data->tstack[ct->numseq[i]][ct->numseq[j]]
					[ct->numseq[i+1]][ct->numseq[j-1]]) +
					erg4_nc (ip,jp,jp+1,1,ct,data,false)+penalty_nc(ip,jp,ct,data);

				energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstack[ct->numseq[jp]][ct->numseq[ip]]
					[ct->numseq[jp+1]][ct->numseq[ip-1]]) +
					erg4_nc (i,j,j-1,2,ct,data,false)+penalty_nc(i,j,ct,data);

#ifdef COUNTING				 
				  if (energy2 < energy){
					energy_option=2;
					energy = energy2;
				  }		 
#else				  
				  energy = min (energy,energy2);
#endif //COUNTING
				//if ((ct->numseq[i+1]!=5)&&(ct->numseq[ip-1]!=5)) {
				//now consider if coaxial stacking is better:
				energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstackcoax[ct->numseq[i]]
					[ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]])
					+ NO_COUNT(data->coaxstack[ct->numseq[i+1]][ct->numseq[j-1]]
					[ct->numseq[ip]][ct->numseq[jp]])+penalty_nc(i,j,ct,data)+penalty_nc(jp,ip,ct,data);
 		
#ifdef COUNTING				 
				  if (energy2 < energy){
					energy_option=3;
					energy = energy2;
				  }		 
#else				  
				  energy = min (energy,energy2);
#endif //COUNTING

				energy2 = NO_COUNT(data->init) + NO_COUNT(data->tstackcoax[ct->numseq[i]]
					[ct->numseq[j]][ct->numseq[ip-1]][ct->numseq[j-1]])
					+ NO_COUNT(data->coaxstack[ct->numseq[ip-1]][ct->numseq[j-1]]
					[ct->numseq[ip]][ct->numseq[jp]])+penalty_nc(i,j,ct,data)+penalty_nc(jp,ip,ct,data);
  
#ifdef COUNTING				 
				  if (energy2 < energy){
					energy_option=4;
					energy = energy2;
				  }		 

				  switch (energy_option){
					case 1:
						energy = data->init + data->tstack[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] +
							erg4 (ip,jp,jp+1,1,ct,data,false)+penalty(ip,jp,ct,data);
						break;

					case 2:
						energy = data->init + data->tstack[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] +
							erg4 (i,j,j-1,2,ct,data,false)+penalty(i,j,ct,data);
						break;

					case 3:
						energy = data->init + data->tstackcoax[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]]
							+ data->coaxstack[ct->numseq[i+1]][ct->numseq[j-1]][ct->numseq[ip]][ct->numseq[jp]]+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
						break;

					case 4:
						energy = data->init + data->tstackcoax[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip-1]][ct->numseq[j-1]]
							+ data->coaxstack[ct->numseq[ip-1]][ct->numseq[j-1]][ct->numseq[ip]][ct->numseq[jp]]+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
						break;
				}
#else				  
				  energy = min (energy,energy2);
#endif //COUNTING

				
               }
               else if (size1==0) {//just have dangling ends or flush stacking
				energy = NO_COUNT(data->init) + erg4_nc(jp,ip,jp+1,1,ct,data,false) +
					erg4_nc(i,j,j-1,2,ct,data,false)+penalty_nc(i,j,ct,data)+penalty_nc(jp,ip,ct,data);
				
				energy2 = NO_COUNT(data->init) + NO_COUNT(data->coax[ct->numseq[j]][ct->numseq[i]]
					[ct->numseq[ip]][ct->numseq[jp]])+penalty_nc(i,j,ct,data)+penalty_nc(jp,ip,ct,data);

#ifdef COUNTING				
				if (energy < energy2) {
					energy = data->init + erg4 (jp,ip,jp+1,1,ct,data,false) +
						erg4 (i,j,j-1,2,ct,data,false)+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
				}
				else {
					energy = data->init + data->coax[ct->numseq[j]][ct->numseq[i]]
						[ct->numseq[ip]][ct->numseq[jp]]+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
				}
#else
                  energy = min(energy,energy2);
#endif //COUNTING				  
               }


         		return energy;

         }
      }


      //a typical internal or bulge loop:
		//size1 = ip-i-1;
		//size2 = j - jp - 1;

		//Introduces single stranded pseudoenergies from SHAPE data for interior/bulge loops
		int SHAPEss_energy = 0;

		if (size1 == 1) SHAPEss_energy += ct->SHAPEss_give_value(i + 1);
		else if (size1 != 0) SHAPEss_energy += ct->SHAPEss_calc(i + 1, ip - 1);

		if (size2 == 1) SHAPEss_energy += ct->SHAPEss_give_value(j - 1);
		else if (size2 != 0) SHAPEss_energy += ct->SHAPEss_calc(jp + 1, j - 1);


		if (size1==0||size2==0) {//bulge loop


			size = size1+size2;
			if (size==1) {
				count = 1;
				energy = data->stack[ct->numseq[i]][ct->numseq[j]]
						[ct->numseq[ip]][ct->numseq[jp]]
						+ data->bulge[size] + data->eparam[2];
						

				if (size1==1)  {

					//count the number of alternative bulges that exist:

					k = i;
					while (ct->numseq[k]==ct->numseq[i+1]) {
						count++;
						k--;

						//During suboptimal structure prediction (where sequence is doubled), make sure this doesn't walk from one end of the sequence to the other.
						if (k==ct->GetSequenceLength()||k==0) break;
					}

					k=ip;
					while(ct->numseq[k]==ct->numseq[i+1]) {
						count++;
						k++;

						//During suboptimal structure prediction (where sequence is doubled), make sure this doesn't walk from one end of the sequence to the other.
						if (k==ct->GetSequenceLength()+1||k>2*ct->GetSequenceLength()) break;
					}
					//give bonus to C adjacent to single C bulge
					if ((ct->IsNuc(i+1,'C')||ct->IsNuc(i+1,'c'))&&count>1) energy+= data->singlecbulge;

				}

				else {
					//size2 == 1

					//count the number of alternative bulges that exist:

					k = jp;
					while (ct->numseq[k]==ct->numseq[jp+1]) {
						count++;
						k--;

						//During suboptimal structure prediction (where sequence is doubled), make sure this doesn't walk from one end of the sequence to the other.
						if (k==ct->GetSequenceLength()||k==0) break;
					}
					k=j;
					while(ct->numseq[k]==ct->numseq[jp+1]) {
						count++;
						k++;

						//During suboptimal structure prediction (where sequence is doubled), make sure this doesn't walk from one end of the sequence to the other.
						if (k==ct->GetSequenceLength()+1||k>2*ct->GetSequenceLength()) break;
					}
					//give bonus to C adjacent to single C bulge
					if ((ct->IsNuc(j-1,'C')||ct->IsNuc(j-1,'c'))&&count>1) energy+= data->singlecbulge;

				}
				//apply a correction for the number of equivalent states because
					//the bulge can move to adjacent sites
				energy-= (int) (data->RT*conversionfactor* log ((double) count));
			}
			else if (size>30) {

				loginc = int((data->prelog)*log(double ((size)/30.0)));
				energy = data->bulge[30] + loginc + data->eparam[2];
				energy = energy + penalty(i,j,ct,data) + penalty(jp,ip,ct,data);

			}
			else {
         		energy = data->bulge[size] + data->eparam[2];
				energy = energy + penalty(i,j,ct,data) + penalty(jp,ip,ct,data);
			}
		}
		else {//internal loop
			size = size1 + size2;
			lopsid = abs(size1-size2);

			if (size>30) {

				loginc = int((data->prelog)*log((double ((size))/30.0)));
				if (size1==1||size2==1) {
            		energy = data->tstki1n[ct->numseq[i]][ct->numseq[j]]
						[ct->numseq[i+1]][ct->numseq[j-1]] +
						data->tstki1n[ct->numseq[jp]][ct->numseq[ip]]
						[ct->numseq[jp+1]][ct->numseq[ip-1]] +
						data->inter[30] + loginc + data->eparam[3];
#ifdef COUNTING			
					if (NO_COUNT(data->maxpen) < lopsid*NO_COUNT(data->poppen[min(2,min(size1,size2))])) energy += data->maxpen;
					else{
						energy += lopsid*data->poppen[min(2,min(size1,size2))];
						data->poppen[min(2,min(size1,size2))].get+=lopsid-1;						
					}
#else
					energy += min(data->maxpen,lopsid*data->poppen[min(2,min(size1,size2))]);
#endif					
				}

				else {
					energy = data->tstki[ct->numseq[i]][ct->numseq[j]]
						[ct->numseq[i+1]][ct->numseq[j-1]] +
						data->tstki[ct->numseq[jp]][ct->numseq[ip]]
						[ct->numseq[jp+1]][ct->numseq[ip-1]] +
						data->inter[30] + loginc + data->eparam[3];
#ifdef COUNTING						
						if (NO_COUNT(data->maxpen) < lopsid*NO_COUNT(data->poppen[min(2,min(size1,size2))])) energy += data->maxpen;
						else{
							energy += lopsid*data->poppen[min(2,min(size1,size2))];
							data->poppen[min(2,min(size1,size2))].get+=lopsid-1;
						}
#else
						energy += min(data->maxpen,lopsid*data->poppen[min(2,min(size1,size2))]);
#endif //COUNTING						
						
				}
			}
			else if ((size1==2)&&(size2==2))//2x2 internal loop
			    energy = data->iloop22[ct->numseq[i]][ct->numseq[ip]]
					[ct->numseq[j]][ct->numseq[jp]]
					[ct->numseq[i+1]][ct->numseq[i+2]]
					[ct->numseq[j-1]][ct->numseq[j-2]];


			else if ((size1==1)&&(size2==2)) {//2x1 internal loop
				energy = data->iloop21[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]]
					[ct->numseq[j-1]][ct->numseq[jp+1]][ct->numseq[ip]][ct->numseq[jp]];


			}
			else if ((size1==2)&&(size2==1)) {//1x2 internal loop
				energy = data->iloop21[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]]
					[ct->numseq[ip-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]];

			}

			else if (size==2) //a single mismatch

				energy = data->iloop11[ct->numseq[i]][ct->numseq[i+1]][ct->numseq[ip]]
					[ct->numseq[j]][ct->numseq[j-1]][ct->numseq[jp]];
			else if (size1==1||size2==1) { //this loop is lopsided
         	//this is a 1xn loop:
				energy = data->tstki1n[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] +
						data->tstki1n[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] +
						data->inter[size] + data->eparam[3];
#ifdef COUNTING					
				if (NO_COUNT(data->maxpen) < lopsid*NO_COUNT(data->poppen[min(2,min(size1,size2))])) energy += data->maxpen;
				else{
					energy += lopsid*data->poppen[min(2,min(size1,size2))];
					
					data->poppen[min(2,min(size1,size2))].get+=lopsid-1;
				}  
#else	
				energy += min(data->maxpen,(lopsid*data->poppen[min(2,min(size1,size2))]));
#endif	
			}


			else if ((size1==2&&size2==3)||(size1==3&&size2==2)) {
			//this is a 2x3 loop
				energy = data->tstki23[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] +
					data->tstki23[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] +
					data->inter[size] + data->eparam[3];

#ifdef COUNTING					
				if (NO_COUNT(data->maxpen) < lopsid*NO_COUNT(data->poppen[min(2,min(size1,size2))])) energy += data->maxpen;
				else{
					energy += lopsid*data->poppen[min(2,min(size1,size2))];
					
					data->poppen[min(2,min(size1,size2))].get+=lopsid-1;
				}  
#else	
				energy += min(data->maxpen,(lopsid*data->poppen[min(2,min(size1,size2))]));
#endif					


			}
			else {



         		energy = data->tstki[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] +
					data->tstki[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] +
					data->inter[size] + data->eparam[3];

#ifdef COUNTING						
				if (NO_COUNT(data->maxpen) < lopsid*NO_COUNT(data->poppen[min(2,min(size1,size2))])) energy += data->maxpen;
				else{
					energy += lopsid*data->poppen[min(2,min(size1,size2))];
				
					data->poppen[min(2,min(size1,size2))].get+=lopsid-1;
					
				}
#else
				energy += min(data->maxpen,(lopsid*data->poppen[min(2,min(size1,size2))]));
#endif //COUNTING				
			}
			//energy+=(SHAPEend(i,ct)+SHAPEend(j,ct)+SHAPEend(ip,ct)+SHAPEend(jp,ct));

		}

		//adds SHAPE pseudoenergy to energy value

		energy += SHAPEss_energy;

		return energy;
}

//calculate the energy of a hairpin loop:
integersize erg3(int i,int j,structure *ct, datatable *data,char dbl)
{
	integersize energy;
	int size,loginc,count,key,k;
	const int alphabetSize = data->alphabet.size();


	if ((i<=(ct->GetSequenceLength()))&&(j>(ct->GetSequenceLength()))) {
      	//A hairpin cannot contain the ends of the sequence
         return INFINITE_ENERGY;
    }

	if (dbl&DUBLE) return INFINITE_ENERGY;//the loop contains a base that should be
      										//double stranded

    else if (dbl&INTER) {//intermolecular interaction
      	//intermolecular "hairpin" free energy is that of intermolecular
         //	initiation plus the stacked mismatch
		energy = data->init;

#ifdef COUNTING
         if (NO_COUNT(data->tstack[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]]) < erg4_nc(i,j,i+1,1,ct,data,false)) {
         	energy += data->tstack[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]];
         }
         else{
         	energy += erg4(i,j,i+1,1,ct,data,false);
         }
		 energy += penalty(i,j,ct,data);
#else
         energy += min(data->tstack[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]],erg4(i,j,i+1,1,ct,data,false))+penalty(i,j,ct,data);
#endif
         return energy;
    }




		size = j-i-1;



		if (size>30) {

			loginc = int((data->prelog)*log((double ((size))/30.0)));



			energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]]
				+ data->hairpin[30]+loginc+data->eparam[4];
		}
		else if (size<3) {
      		energy = data->hairpin[size] + data->eparam[4];
			energy+=penalty(i,j,ct,data);
		}
		else if (size==4) {

			//key = (ct->numseq[j])*3125 + (ct->numseq[i+4])*625 +
			//	(ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*alphabetSize+(ct->numseq[i]);
			int digit_factor = 1;
			key = 0;
			for (count = 0; count < 6; ++count) {
				key += ct->numseq[i+count] * digit_factor;
				digit_factor = digit_factor * alphabetSize;
			}
			for (count=0;count<data->numoftloops;count++) {
				
					if (key==data->tloop[count][0]) {
						return data->tloop[count][1];
					}
			}
			
			energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]]
				+ data->hairpin[size] + data->eparam[4];
		}
		else if (size==3) {

			//key = (ct->numseq[j])*625 +
			//	(ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*5+(ct->numseq[i]);
			int digit_factor = 1;
			key = 0;
			for (count = 0; count < 5; ++count) {
				key += ct->numseq[i+count] * digit_factor;
				digit_factor = digit_factor * alphabetSize;
			}
			for (count=0;count<data->numoftriloops;count++) {
				if (key==data->triloop[count][0]) return data->triloop[count][1];
			}

			energy =	data->hairpin[size] + data->eparam[4]
         	+penalty(i,j,ct,data);
		}
		else if (size==6) {
			//key = (ct->numseq[j])*78125 + (ct->numseq[i+6])*15625 + (ct->numseq[i+5])*3125
			//	+ (ct->numseq[i+4])*625 +
			//	(ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*5+(ct->numseq[i]);
			int digit_factor = 1;
			key = 0;
			for (count = 0; count < 8; ++count) {
				key += ct->numseq[i+count] * digit_factor;
				digit_factor = digit_factor * alphabetSize;
			}
			for (count=0;count<data->numofhexaloops;count++) {
				if (key==data->hexaloop[count][0]) {
					return data->hexaloop[count][1];
				}
			}

			energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]]
				+ data->hairpin[size] + data->eparam[4];
		}

		else {
			energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]]
				+ data->hairpin[size] + data->eparam[4];
		}




		//check for GU closeure preceded by GG

		if (ct->IsNuc(i,'G')||ct->IsNuc(i,'g')) {
			if (ct->IsNuc(j,'U')||ct->IsNuc(j,'u')) {

		 
      			if ((i>2&&i<ct->GetSequenceLength())||(i>ct->GetSequenceLength()+2)) {
				
					if (ct->IsNuc(i-1,'G')||ct->IsNuc(i-1,'g')) {
						if (ct->IsNuc(i-2,'G')||ct->IsNuc(i-2,'g')) {
				

         					energy = energy + data->gubonus;


						}	
					}//end if (ct->IsNuc(i-1,'G')||ct->IsNuc(i-1,'g'))
				}//end if ((i>2&&i<ct->GetSequenceLength())||(i>ct->GetSequenceLength()+2))
			}//end if (ct->IsNuc(j,'U')||ct->IsNuc(j,'u')) 
		}//end if (ct->IsNuc(i,'G')||ct->IsNuc(i,'g'))

	  //adds pseudo energy term for a hairpin loop based off of SHAPE data
	  if (ct->shaped) energy += ct->SHAPEss_calc(i + 1, j - 1);

      //check for an oligo-c loop

      for (k=1;(k<=size);k++) {
       	if (ct->numseq[i+k] != 2) return energy;//this is not an oligo-C loop
      }
      //this is a poly c loop so penalize
      if (size==3) return (energy + data->c3);
      else return (energy + data->cint + size*data->cslope);


}



//calculate the energy of a dangling end:
integersize erg4(int i,int j,int ip,int jp,structure *ct, datatable *data, bool lfce)
{
	integersize energy;


	//dangling base
		// jp = 1 => 3' dangle
		// jp = 2 => 5' dangle



      if (lfce) return INFINITE_ENERGY;//stacked nuc should be double stranded

	  //commented out 11/8/99
      //if (ip==5) return 0;//dangling nuc is an intermolecular linker

		energy = data->dangle[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][jp];

		//also add SAHPE single stranded contribution:
		energy+=ct->SHAPEss_give_value(ip);

		return energy;

}


integersize erg4_nc(int i,int j,int ip,int jp,structure *ct, datatable *data, bool lfce)
{
	integersize energy;

#ifdef	COUNTING
	int tmp_count = data->dangle[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][jp].get;
#endif

	energy = erg4(i,j,ip,jp,ct,data,lfce);

#ifdef	COUNTING
	data->dangle[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][jp].get = tmp_count;
#endif

	return energy;
}



//these functions calculate the free energy for coaxial stacking of pair i-j onto ip-jp
//	require that the backbone continues directly between j and ip without a nucleotide in
//	a canonical pair
//k (determined by the function called) indicates the intervening nuc in intervening mismatch that is not sandwiched between
//	j and ip
integersize ergcoaxflushbases(int i, int j, int ip, int jp, structure *ct, datatable *data) {
	//flush stacking, k==0

		//Remapped 10/11/2013 to match the order of the helical stack table
		return data->coax[ct->numseq[j]][ct->numseq[i]][ct->numseq[ip]][ct->numseq[jp]];

}
integersize ergcoaxinterbases1(int i, int j, int ip, int jp, structure *ct, datatable *data) {
	//k==i-1
	return data->tstackcoax[ct->numseq[j]][ct->numseq[i]][ct->numseq[j+1]][ct->numseq[i-1]] +
				data->coaxstack[ct->numseq[j+1]][ct->numseq[i-1]][ct->numseq[ip]][ct->numseq[jp]] +
				ct->SHAPEss_give_value(j+1) + ct->SHAPEss_give_value(i-1);

}

integersize ergcoaxinterbases2(int i, int j, int ip, int jp, structure *ct, datatable *data) {
	//k==jp+1
	return data->tstackcoax[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] +
				data->coaxstack[ct->numseq[j]][ct->numseq[i]][ct->numseq[j+1]][ct->numseq[jp+1]] +
				ct->SHAPEss_give_value(jp+1) + ct->SHAPEss_give_value(ip-1);

}


//this function calculates the free energy for coaxial stacking of pair i-j onto ip-jp
//	require that the backbone continues directly between j and ip without a nucleotide in
//	a canonical pair


//Note that this form of the function takes the sequences, not the number of the nuc
//	in i,j,ip,and jp.
integersize ergcoaxflushbases(int i, int j, int ip, int jp, datatable *data) {


		

		//Look up the energy in the coax array.  
		//This was remapped 10/10/2013 to be consistent in order with the stack array.
		return data->coax[j][i][ip][jp];

}



//this function calculates the free energy for coaxial stacking of pair i-j onto ip-jp
//	require that the backbone continues directly between j and ip without a nucleotide in
//	a canonical pair
//k indicates the intervening nuc in intervening mismatch that is not sandwiched between
//	j and ip
//l is the nuc sandwiched between j and ip


//this requires a backbone like:
// j-l-ip
// | . |
// i-k jp
//i.e. a discontinuity between k and jp


//Note that this form of the function takes the sequences, not the number of the nuc
//	in i,j,ip,jp,k and l.
integersize ergcoaxinterbases1(int i, int j, int ip, int jp, int k, int l, datatable *data) {

		//coaxial stacking with an intervening mismatch

			return data->tstackcoax[j][i][l][k] +
				data->coaxstack[l][k][ip][jp];



}


//this function calculates the free energy for coaxial stacking of pair i-j onto ip-jp
//	require that the backbone continues directly between j and ip without a nucleotide in
//	a canonical pair
//k indicates the intervening nuc in intervening mismatch that is not sandwiched between
//	j and ip
//l is the nuc sandwiched between j and ip


//this requires a backbone like:
// j-l-ip
// | . |
// i k-jp
//i.e. a discontinuity between i and jp


//Note that this form of the function takes the sequences, not the number of the nuc
//	in i,j,ip,jp,k and l.
integersize ergcoaxinterbases2(int i, int j, int ip, int jp, int k, int l, datatable *data) {

		//coaxial stacking with an intervening mismatch

			return data->tstackcoax[jp][ip][k][l] +
				data->coaxstack[j][i][l][k];



}


//This function will calculate the free energy of a multiloop starting at nuc i for
//	structure #st
//This uses a recursive algorithm

#define js true //this switch can turn on and off the logarithmic dep. on
				//unpaired nucleotides.  This is helpful for debugging

//simplemb, when true, indicates that the logarithmic dependence should be off.
//This sets the energy function equal to that used by the dynamic programming
//algorithms.

//This function will calculate the free energy of an exterior loop
//  for structure #st
//This uses a recursive algorithm


