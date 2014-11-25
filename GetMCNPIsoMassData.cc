#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <dirent.h>
#include <iomanip>

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
    string inFileName, outDirName, fileName , test;
    int result=0;

    std::vector<double> isotopeMass(119, 0);
    std::vector<int> elemNumIso(119, 0);
    std::vector<int> elemBaseA(119, 0);
    std::vector<int> elemIsoIndex;
    elemIsoIndex.reserve(119);
    stringstream numConv;

    for(int i=0; i<119; i++)
    {
        elemIsoIndex.push_back(i);
    }

    elemNumIso[0]=1;

    stringstream stream;

    if(argc==4)
    {
        stream << argv[1] << ' ' << argv[2] << ' ' << argv[3] ;
        stream >> inFileName >> outDirName >> version;
    }
    else if(argc==3)
    {
        stream << argv[1] << ' ' << argv[2];
        stream >> inFileName >> outDirName;
    }
    else
    {
        cout << "Incorrect number of inputs; give the MCNP files for the data to be extracted from and the output directory name" << endl;
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

    stream << "\n" << endl;

    stream << "void IsotopeMass::SetIsotopeMass()\n{\n\telemNumIso = new int [" << arraySize << "];\n\telemBaseA = new int [" << arraySize << "];\n\tisotopeMass = new double *[" << arraySize << "];\n" << endl;

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

    stream << "\n}" << endl;

    string outFileName=CreateMacroName(inFileName, outDirName);
    SetDataStream(outFileName, stream);

    if(result>1)
        result=1;
    return result;
}

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
