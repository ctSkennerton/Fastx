/*
 *  Created by Connor Skennerton on 10/08/11.
 *  Copyright 2011 Connor Skennerton. All rights reserved. 
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 *
 *                     A B R A K A D A B R A
 *                      A B R A K A D A B R
 *                       A B R A K A D A B
 *                        A B R A K A D A       	
 *                         A B R A K A D
 *                          A B R A K A
 *                           A B R A K
 *                            A B R A
 *                             A B R
 *                              A B
 *                               A
 */

#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include "Fastx.h"

std::istream& Fastx::read(std::istream& stream)
{
    this->clear();
    try 
    {
        char c = stream.get();
        if ('>' == c )
        {
            FX_Fasta = true;
        }
        else if ('\n' == c)
        {
            std::cout<<"new line"<<std::endl;
        }
        else if ('@' != c) 
        {
            throw "badFormat! The first char in a file should be a '>' or '@' for fasta or fastq files";
        }
    } catch (char * c) {
        std::cerr<<c<<std::endl;
    }
    
 
    getline(stream,FX_Header);

    if (FX_Fasta) 
    {        
        char c;
        while (stream >> c) 
        {
            if (c == '>' || stream.eof()) 
            {
                stream.putback(c);
                break;
            }
            else if (c != '\n')
            {
                FX_Sequence += c; 
            }
        }
    }
    // it's Fastq
    // only thinks about sequencing type 4 line files
    else
    {
        getline(stream,FX_Sequence);
        getline(stream,FX_Comment);
        getline(stream,FX_Quality);

    }   
    FX_SeqLength = FX_Sequence.length();
    return stream;
}

std::ostream& Fastx::print (std::ostream& s)
{
   if (FX_Fasta) {
       s<<'>'<<FX_Header<<std::endl<<FX_Sequence;
   } else {
       s<<'@'<<FX_Header<<std::endl<<FX_Sequence<<std::endl<<FX_Comment<<std::endl<<FX_Quality;
   }
    return s;
}
               
float Fastx::GCContent(void)
{
    std::string::iterator seq_iter = FX_Sequence.begin();
    float count = 0.0;
    while (seq_iter != FX_Sequence.end()) 
    {
        switch (*seq_iter) 
        {
            case 'G':
            case 'g':
            case 'C':
            case 'c':
                count++;
                break;
            default:
                break;
        }
        seq_iter++;
    }
    return count / (float)FX_SeqLength;
}


float Fastx::GCSkew(void)
{
    float g = 0.0;
    float c = 0.0;
    std::string::iterator seq_iter = FX_Sequence.begin();
    while (seq_iter != FX_Sequence.end()) 
    {
        switch (*seq_iter) 
        {
            case 'G':
            case 'g':
                g++;
                break;
            case 'C':
            case 'c':
                c++;
                break;
            default:
                break;
        }
        seq_iter++;
    }
    return ((g-c)/(g+c));
}

float Fastx::ATContent(void)
{
    std::string::iterator seq_iter = FX_Sequence.begin();
    float count = 0.0;
    while (seq_iter != FX_Sequence.end()) 
    {
        switch (*seq_iter) 
        {
            case 'A':
            case 'a':
            case 'T':
            case 't':
                count++;
                break;
            default:
                break;
        }
        seq_iter++;
    }
    return count / (float)FX_SeqLength;
}


float Fastx::ATSkew(void)
{
    float a = 0.0;
    float t = 0.0;
    std::string::iterator seq_iter = FX_Sequence.begin();
    while (seq_iter != FX_Sequence.end()) 
    {
        switch (*seq_iter) 
        {
            case 'A':
            case 'a':
                a++;
                break;
            case 'T':
            case 't':
                t++;
                break;
            default:
                break;
        }
        seq_iter++;
    }
    return ((a-t)/(a+t));
}
Fastx Fastx::subseq(int b, int l)
{
    Fastx tmp;
    tmp.FX_Sequence = this->FX_Sequence.substr(b,l);
    tmp.FX_Header = this->FX_Header;
    tmp.FX_SeqLength = tmp.FX_Sequence.length();
    if (!FX_Fasta) {
        tmp.FX_Quality = this->FX_Quality.substr(b,l);
        tmp.FX_Comment = this->FX_Comment;
    }
    
    return tmp;
    
}

Fastx Fastx::subseq(int b)
{
    Fastx tmp;
    tmp.FX_Sequence = this->FX_Sequence.substr(b);
    tmp.FX_Header = this->FX_Header;
    tmp.FX_SeqLength = tmp.FX_Sequence.length();
    if (!FX_Fasta) {
        tmp.FX_Quality = this->FX_Quality.substr(b);
        tmp.FX_Comment = this->FX_Comment;
    }
    return tmp;
}

Fastx Fastx::reverseComplement()
{
    Fastx tmp;
    tmp.FX_Header = this->FX_Header;
    std::stringstream revcomp_string;
    std::string::reverse_iterator rit;
    for ( rit = FX_Sequence.rbegin() ; rit < FX_Sequence.rend(); rit++ )
    {
		switch ((*rit)) 
        {
            case 'A':
            case 'a':
                revcomp_string << 'T';
                break;
            case 'C':
            case 'c':
                revcomp_string << 'G';
                break;
            case 'G':
            case 'g':
                revcomp_string << 'C';
                break;
            case 'T':
            case 't':
            case 'U':
            case 'u':
                revcomp_string << 'A';
                break;
            case 'N':
            case 'n':
                revcomp_string << 'N';
                break;
            case 'M':
            case 'm':
                revcomp_string <<'K'; 
                break;
            case 'R':
            case 'r':
                revcomp_string << 'Y';
                break;
            case 'W':
            case 'w':
                revcomp_string << 'W';
                break;
            case 'S':
            case 's':
                revcomp_string << 'S';
                break;
            case 'Y':
            case 'y':
                revcomp_string << 'R';
                break;
            case 'K':
            case 'k':
                revcomp_string << 'M';
                break;
            case 'V':
            case 'v':
                revcomp_string << 'B';
                break;
            case 'H':
            case 'h':
                revcomp_string << 'D';
                break;
            case 'D':
            case 'd':
                revcomp_string << 'H';
                break;
            case 'B':
            case 'b':
                revcomp_string << 'V';
                break;
            default:
                revcomp_string << 'N';
                break;
		}
	}
    tmp.FX_Sequence =revcomp_string.str();
    if (!FX_Fasta) {
        std::string::reverse_iterator rit = FX_Quality.rbegin();
        while (rit != FX_Quality.rend()) {
            tmp.FX_Quality += *rit;
            rit++;
        }
        tmp.FX_Comment = this->FX_Comment;
    }
    return tmp;
}

Fastx Fastx::laurenize()
{
    Fastx seq2 = this->reverseComplement();
    if ((*this) < seq2)
    {
        return (*this);
    }
    return seq2;
}
// change the single ascii values into their corespnding Phread quality score
std::string Fastx::qualityToPhreadScoreAsString(qualType qt)
{
    std::stringstream phred;
    std::string::iterator iter = FX_Quality.begin();
    int Q;
    while (iter != FX_Quality.end()) 
    {
        switch (qt) 
        {
            case sanger:
                Q = (int)(*iter - 33);
                phred << Q;
                break;
            case illumina:
                Q = (int)(*iter - 64);
                phred << Q;
                break;
            default:
                break;
        }
        iter++;
    }
    return phred.str();
    
}

// change the single ascii values into their corespnding Phread quality score
std::vector<int> Fastx::qualityToPhreadScoreAsVector(qualType qt)
{
    std::vector<int> phred;
    std::string::iterator iter = FX_Quality.begin();
    int Q;
    while (iter != FX_Quality.end()) 
    {
        switch (qt) 
        {
            case sanger:
                Q = (int)(*iter - 33);
                phred.push_back(Q);
                break;
            case illumina:
                Q = (int)(*iter - 64);
                phred.push_back(Q);
                break;
            default:
                break;
        }
        iter++;
    }
    return phred;
    
}


void Fastx::convertPhreadScoreToQuality(qualType qt)
{
    //std::stringstream qual;
    std::string::iterator iter = FX_Quality.begin();
    while (iter != FX_Quality.end()) 
    {
        //$q = chr(($Q<=93? $Q : 93) + 33);
        switch (qt) 
        {
            case sanger:
                *iter = (char)(((*iter) <= 93 ? (*iter) : 93) + 33);
                break;
            case illumina:
                *iter = (char)(((*iter) <= 62 ? (*iter) : 62) + 64);
                break;
            default:
                break;
        }
        iter++;
    }
    
}



// change the current ascii quality character into the corresponding Sanger type quality character
void Fastx::convertQualityToSangerAsciiValues(qualType qt)
{
    std::string::iterator qual_iter;
    switch (qt) 
    {
        case illumina:
            qual_iter = FX_Quality.begin();
            while (qual_iter != FX_Quality.end()) 
            {
                *qual_iter =(char) (*qual_iter - 31);
                qual_iter++;
            }
            break;
        case sanger:
            break;
        default:
            break;
    }
    
}


// change the current quality values into their corresponding illumina values
void Fastx::convertQualityToIlluminaAsciiValues(qualType qt)
{
    std::string::iterator qual_iter;
    switch (qt) 
    {
        case sanger:
            qual_iter = FX_Quality.begin();
            while (qual_iter != FX_Quality.end()) 
            {
                *qual_iter = (char)(*qual_iter + 31);
                qual_iter++;
            }
            break;
        case illumina:
            break;
        default:
            break;
    }
}


std::ostream & operator<< (std::ostream& s, Fastx * c)
{
    return c->print (s);
    
}

std::istream & operator>> (std::istream& s, Fastx * c)
{
    return c->read (s);
}


std::ostream& operator<< (std::ostream& s, Fastx& c)
{
    return c.print(s);
    
}

std::istream& operator>> (std::istream & s, Fastx& c)
{
    return c.read(s);
}

bool operator< (Fastx& f1, Fastx& f2)
{
    if (f1.FX_Sequence < f2.FX_Sequence) 
    {
        return true;
    }
    return false;
}

bool operator> (Fastx& f1, Fastx& f2)
{
    if (f1.FX_Sequence > f2.FX_Sequence) 
    {
        return true;
    }
    return false;
}