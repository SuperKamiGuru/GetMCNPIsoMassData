#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <dirent.h>
#include <iomanip>
#include "include/ElementNames.hh"

#define nMass 1.008664916


using namespace std;

int ProcessFile(stringstream &stream, std::vector<double> &isotopeMass, std::vector<int> &elemNumIso, std::vector<int> &elemBaseA,
                std::vector<int> &elemIsoIndex, char lib);
int GetIsoMass(stringstream &stream, int &Z, int &A, double &mass);
void GetDataStream( string, std::stringstream&);
string CreateMacroName(string geoFileName, string outDirName);
void SetDataStream( string, std::stringstream&);

int main(int argc, char **argv)
{

    string word;
    char lib, version='7';
    string inFileName, outDirName, natAbunFile, fileName , test;
    int result=0;

    std::vector<double> isotopeMass(119, 0);
    std::vector<int> elemNumIso(119, 0);
    std::vector<int> elemBaseA(119, 0);
    std::vector<int> elemIsoIndex;
    elemIsoIndex.reserve(119);
    stringstream numConv;

    ElementNames elemNames;
    elemNames.SetElementNames();

    for(int i=0; i<119; i++)
    {
        elemIsoIndex.push_back(i);
    }

    elemNumIso[0]=1;

    stringstream stream;

    if(argc==5)
    {
        stream << argv[1] << ' ' << argv[2] << ' ' << argv[3] << ' ' << argv[4] ;
        stream >> inFileName >> outDirName >> natAbunFile >> version;
    }
    else if(argc==4)
    {
        stream << argv[1] << ' ' << argv[2] << ' ' << argv[3];
        stream >> inFileName >> outDirName >> natAbunFile;
    }
    else
    {
        cout << "Incorrect number of inputs; give the MCNP files for the data to be extracted from, the output directory name, and the file containing the natural abundances" << endl;
        return 1;
    }
    numConv << "endf" << version;
    numConv >> test;
    numConv.clear();
    numConv.str("");

    stream.clear();
    stream.str("");

    if(inFileName.back()=='/')
    {
        lib=version;
        DIR *dir;
        struct dirent *ent;
        if ((dir = opendir (inFileName.c_str())) != NULL)
        {
            while ((ent = readdir (dir)) != NULL)
            {
                fileName = string(ent->d_name).substr(0,5);

                if(fileName==test)
                {
                    GetDataStream(inFileName+'/'+ent->d_name, stream);
                    result+=ProcessFile(stream, isotopeMass, elemNumIso, elemBaseA, elemIsoIndex, lib);
                    stream.str("");
                    stream.clear();
                }
            }
            closedir(dir);
        }
    }
    else
    {
        GetDataStream(inFileName, stream);
        lib = inFileName[int(inFileName.length()-3)];
        result+=ProcessFile(stream, isotopeMass, elemNumIso, elemBaseA, elemIsoIndex, lib);
        stream.str("");
        stream.clear();
    }

    int arraySize = elemNumIso.size();

    for(int i=0; i<int(elemNumIso.size()); i++)
    {
        if(elemNumIso[i]==0)
            elemNumIso[i]=1;
    }

    stream.str("");
    stream.clear();

    numConv.str("");
    numConv.clear();

    //double **natIsoAbun = new double *[elemNumIso.size()];
    vector <vector<double>> natIsoAbun;
    vector<double> tempVec;
    natIsoAbun.reserve(119);
    for(int i=0; i<int(elemNumIso.size()); i++)
    {
        /*
        if(elemNumIso[i]!=0)
            natIsoAbun[i] = new double [elemNumIso[i]];
        else
            natIsoAbun[i]=NULL;
        for(int j=0; j<int(elemNumIso[i]); j++)
        {
            natIsoAbun[i][j]=0.;
        }
        */
        tempVec.clear();
        for(int j=0; j<int(elemNumIso[i]); j++)
        {
            tempVec.push_back(0.);
        }
        natIsoAbun.push_back(tempVec);
    }

    GetDataStream(natAbunFile, stream);

    char letter;
    int Z=0, A=0;
    double abun=0.;
    while(stream&&(elemBaseA.size()>Z)&&(natIsoAbun.size()>Z))
    {
        letter = stream.peek();
        if((letter>='0')&&(letter<='9'))
        {
            letter = stream.get();
            numConv.str("");
            numConv.clear();
            while((letter>='0')&&(letter<='9'))
            {
                numConv << letter;
                letter = stream.get();
            }
            numConv >> A;

            while(!((letter>='0')&&(letter<='9')))
            {
                letter = stream.get();
            }
            numConv.str("");
            numConv.clear();
            while(((letter>='0')&&(letter<='9'))||(letter=='.'))
            {
                numConv << letter;
                letter = stream.get();
            }
            numConv >> abun;
            if((elemBaseA[Z]<=A)&&(A-elemBaseA[Z]<natIsoAbun[Z].size()))
                natIsoAbun[Z][A-elemBaseA[Z]] = abun;
            else
            {
                cout << "Warning: isotope Z:" << Z << " A:" << A << " Does not exist in the given MCNP data files, but it does in the given isotope natural abundance file " << endl;
                cout << "Using isotope Z:" << Z << " A:" << elemBaseA[Z] << " instead " << endl;
                if(natIsoAbun[Z].size()!=0)
                    natIsoAbun[Z][0] += abun;
            }
        }
        else if(((letter>='a')&&(letter<='z'))||((letter>='A')&&(letter<='Z')))
        {
            stream >> word;
            if(elemNames.CheckName(word))
            {
                for(Z=1; Z<119; Z++)
                {
                    if(elemNames.CheckName(word, Z))
                        break;
                }
                if(Z>107)
                    break;
            }
        }
        else
        {
            letter = stream.get();
        }
    }

    double sum;
    for(int i=0; i<int(elemNumIso.size()); i++)
    {
        sum=0;
        for(int j=0; j<int(elemNumIso[i]); j++)
        {
            sum+=natIsoAbun[i][j];
        }
        if(sum!=0)
        {
            for(int j=0; j<int(elemNumIso[i]); j++)
            {
                natIsoAbun[i][j]/=sum;
            }
        }
    }

    stream.str("");
    stream.clear();

    stream << "\n" << endl;

    stream << "void IsotopeMass::SetIsotopeMass()\n{\n\telemNumIso = new int [" << arraySize << "];\n\telemBaseA = new int [" <<
                arraySize << "];\n\tisotopeMass = new double *[" << arraySize << "];\n\tisoNatAbun = new double *[" << arraySize << "];\n"  << endl;

    for(int i=0; i<int(elemNumIso.size()); i++)
    {
        stream << "\telemNumIso[" << i << "] = " << elemNumIso[i] << ";" << endl;
    }

    stream << "\n" << endl;

    for(int i=0; i<int(elemBaseA.size()); i++)
    {
        stream << "\telemBaseA[" << i << "] = " << elemBaseA[i] << ";" << endl;
    }

    stream << "\n" << endl;

    stream << "\tfor(int i=0; i<" << arraySize << "; i++)\n\t{\n\t\tisotopeMass[i] = new double [elemNumIso[i]];\n\t}\n" << endl;

    int count=0;
    for(int i=0; i<int(elemNumIso.size()); i++)
    {
        for(int j=0; j<int(elemNumIso[i]); j++, count++)
        {
            stream << "\tisotopeMass[" << i << "]" << "[" << j << "] = " << isotopeMass[count] << ";" << endl;
        }
    }

    stream << "\n" << endl;

    stream << "\tfor(int i=0; i<" << arraySize << "; i++)\n\t{\n\t\tisoNatAbun[i] = new double [elemNumIso[i]];\n\t}" << endl;

    for(int i=0; i<int(elemNumIso.size()); i++)
    {
        for(int j=0; j<int(elemNumIso[i]); j++)
        {
            stream << "\tisoNatAbun[" << i << "]" << "[" << j << "] = " << natIsoAbun[i][j] << ";" << endl;
        }
    }

    stream << "\n}" << endl;

    string outFileName=CreateMacroName(inFileName, outDirName);
    SetDataStream(outFileName, stream);

    /*
    for(int i=0; i<int(elemNumIso.size()); i++)
    {
        if(natIsoAbun[i]!=NULL)
            delete [] natIsoAbun[i];
    }

    delete [] natIsoAbun;
    */
    elemNames.ClearStore();

    if(result>1)
        result=1;
    return result;
}

//extracts the isotope mass data contianed in the given MCNP CS data file
int ProcessFile(stringstream &stream, std::vector<double> &isotopeMass, std::vector<int> &elemNumIso, std::vector<int> &elemBaseA,
                std::vector<int> &elemIsoIndex, char lib)
{
    string word;
    char check1, check2, check3;
    char line[256];
    int result=0, offset, Z, A;
    double mass;

    while(stream)
    {
        stream >> word;
        check1=word[int(word.find_last_of('.')+1)];
        check2=word[int(word.find_last_of('.')+2)];
        check3=word[int(word.find_last_of('.')+3)];
        if((check1==lib)&&(check2=='0'))
        {
            if((check3=='c')||(check3=='d'))
            {
                result += GetIsoMass(stream, Z, A, mass);

                offset=0;
                if(elemNumIso[Z]==0)
                {
                    elemBaseA[Z]=A;
                    isotopeMass[elemIsoIndex[Z]]=mass;
                    elemNumIso[Z]++;
                }
                else if(elemBaseA[Z]>A)
                {
                    int oldA = elemBaseA[Z];
                    elemBaseA[Z]=A;
                    isotopeMass.insert(isotopeMass.begin()+elemIsoIndex[Z], mass);
                    offset++;
                    for(int i=0; i<(oldA-elemBaseA[Z]-1); i++, offset++)
                    {
                        isotopeMass.insert(isotopeMass.begin()+elemIsoIndex[Z]+1, 0.);
                    }
                }
                else
                {
                    for(int i=0; i<(A-elemNumIso[Z]-elemBaseA[Z]+1); i++, offset++)
                    {
                        isotopeMass.insert(isotopeMass.begin()+elemNumIso[Z]+elemIsoIndex[Z]+i, 0.);
                    }
                    isotopeMass[A-elemBaseA[Z]+elemIsoIndex[Z]]=mass;
                }
                if(offset>0)
                {
                    for(int i=Z+1; i<119; i++)
                    {
                        elemIsoIndex[i]+=offset;
                    }
                }
                elemNumIso[Z]+=offset;
            }
            else
            {
                stream.getline(line,256);
            }
        }
        else
        {
            stream.getline(line,256);
        }
    }

    return result;
}

//Extracts the mass of the isotope
int GetIsoMass(stringstream &stream, int &Z, int &A, double &mass)
{
    char line[256];
    stringstream numConv;
    int dummy, isoNum;

    stream >> mass;
    mass *= nMass;

    for(int i=0; i<6; i++)
    {
        stream.getline(line,256);
    }

    stream >> dummy >> isoNum;

    Z=floor(isoNum/1000);
    A=isoNum-Z*1000;
    return 0;
}


void GetDataStream( string filename, std::stringstream& ss)
{
   string* data=NULL;

// Use regular text file
      std::ifstream thefData( filename.c_str() , std::ios::in | std::ios::ate );
      if ( thefData.good() )
      {
         int file_size = thefData.tellg();
         thefData.seekg( 0 , std::ios::beg );
         char* filedata = new char[ file_size ];
         while ( thefData )
         {
            thefData.read( filedata , file_size );
         }
         thefData.close();
         data = new string ( filedata , file_size );
         delete [] filedata;
      }
      else
      {
// found no data file
//                 set error bit to the stream
         ss.setstate( std::ios::badbit );
         cout << endl << "### failed to open ascii file " << filename << " ###" << endl;
      }
   if (data != NULL)
   {
        ss.str(*data);
        if(data->back()!='\n')
            ss << "\n";
        ss.seekg( 0 , std::ios::beg );
    }

    if(data!=NULL)
        delete data;
}

string CreateMacroName(string geoFileName, string outDirName)
{
    size_t pos = geoFileName.find_last_of('/');
    size_t pos2 = std::string::npos;
    if(pos == std::string::npos)
        pos=0;
    else
        pos++;

    return (outDirName+"IsoMassData"+geoFileName.substr(pos, pos2-pos));
}


void SetDataStream( string filename , std::stringstream& ss)
{
    //bool cond=true;
// Use regular text file
    string compfilename(filename);

    if(compfilename.substr((compfilename.length()-2),2)==".z")
    {
        compfilename.pop_back();
        compfilename.pop_back();
    }

      std::ofstream out( compfilename.c_str() , std::ios::out | std::ios::trunc );
      if ( ss.good() )
      {
         ss.seekg( 0 , std::ios::end );
         int file_size = ss.tellg();
         ss.seekg( 0 , std::ios::beg );
         char* filedata = new char[ file_size ];
         while ( ss ) {
            ss.read( filedata , file_size );
            if(!file_size)
            {
                cout << "\n #### Error the size of the stringstream is invalid ###" << endl;
                break;
            }
         }
         out.write(filedata, file_size);
         if (out.fail())
        {
            cout << endl << "writing the ascii data to the output file " << compfilename << " failed" << endl
                 << " may not have permission to delete an older version of the file" << endl;
        }
         out.close();
         delete [] filedata;
      }
      else
      {
// found no data file
//                 set error bit to the stream
         ss.setstate( std::ios::badbit );

         cout << endl << "### failed to write to ascii file " << compfilename << " ###" << endl;
      }
   ss.str("");
}
