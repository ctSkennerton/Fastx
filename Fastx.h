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

#ifndef crass_Fastx_h
#define crass_Fastx_h

#include <string>
#include <vector>

class Fastx {
    
    
public:
    
    enum qualType {sanger, illumina};

    //members
    std::string FX_Header;
    std::string FX_Sequence;
    std::string FX_Quality;
    std::string FX_Comment;
    size_t      FX_SeqLength;
    bool        FX_Fasta;
    
    // constructor
    Fastx(void)
    {}
    
    // destructor
    virtual ~Fastx(void)
    {}
    
    void clear(void)
    {
        FX_Comment.clear();
        FX_Header.clear();
        FX_Quality.clear();
        FX_Fasta = false;
        FX_SeqLength = 0;
        FX_Sequence.clear();
    }
    // Getters & Setters
    inline std::string seq(void)
    {
        return FX_Sequence;
    }
    inline void seq(std::string s)
    {
        FX_Sequence = s;
        FX_SeqLength = s.length();
    }
    
    inline void seq(const char * s)
    {
        FX_Sequence = s;
        FX_SeqLength = strlen(s);
    }    
    inline std::string header(void)
    {
        return FX_Header;
    }
    
    inline void header(std::string h)
    {
        FX_Header = h;
    }
    
    inline void header(const char * h)
    {
        FX_Header = h;
    } 
    inline void comment(std::string c)
    {
        FX_Comment = c;
    }
    inline void comment(const char * c)
    {
        FX_Comment = c;
    }
    inline std::string comment(void)
    {
        return FX_Comment;
    }
    inline void quality(std::string q)
    {
        FX_Quality = q;
    }
    inline void quality(const char * q)
    {
        FX_Quality = q;
    }
    inline std::string quality(void)
    {
        return FX_Quality;
    }
    inline size_t length(void)
    {
        return FX_SeqLength;
    }

    
    
    // member functions
    
    std::istream& read(std::istream& stream);
    std::ostream& print (std::ostream& s);

    
    inline std::string substr(int begin, int length)
    {
        return FX_Sequence.substr(begin, length);
    }
    
    inline std::string substr(int begin)
    {
        return FX_Sequence.substr(begin);
    }
    
    inline char nucleotideAt(int ref)
    {
        return FX_Sequence.at(ref);
    }
    
    inline char qualityAt(int ref)
    {
        return FX_Quality.at(ref);
    }
    
    Fastx subseq(int begin, int length);
    
    Fastx subseq(int begin);
    
    inline Fastx truncate(int begin)
    {
        return subseq(begin);
    }
    
    
    float GCContent(void);
    
    inline int GCPercent(void)
    {
        return (int)(GCContent() * 100);
    }
    
    float GCSkew(void);
    
    float ATContent(void);
    
    inline int ATPercent(void)
    {
        return (int)(ATContent() * 100);
    }
    
    float ATSkew(void);
    
    
    Fastx reverseComplement(void);
    
    Fastx laurenize(void);
    
    inline Fastx lowestLexicographicalForm(void)
    {
        return laurenize();
    }
    
    void convertQualityToSangerAsciiValues(qualType qt);
    
    inline void convertQualityToSangerAsciiValues(void)
    {
        convertQualityToSangerAsciiValues(illumina);
    }     
    
    void convertQualityToIlluminaAsciiValues(qualType qt);
    
    inline void convertQualityToIlluminaAsciiValues(void)
    {
        convertQualityToIlluminaAsciiValues(sanger);
    } 
    
    std::vector<int> qualityToPhreadScoreAsVector(qualType qt);
    
    inline std::vector<int> qualityToPhreadScoreAsVector(void)
    {
        qualityToPhreadScoreAsVector(sanger);
    }    
    std::string qualityToPhreadScoreAsString(qualType qt);
    
    inline std::string qualityToPhreadScoreAsString(void)
    {
        qualityToPhreadScoreAsString(sanger);
    }
    void convertPhreadScoreToQuality(qualType qt);
    
    inline void convertPhreadScoreToQuality(void)
    {
        return convertPhreadScoreToQuality(sanger);
    }

};



// overloaded operators to deal with normal variables
std::ostream& operator<< (std::ostream& s,  Fastx& c);

std::istream& operator>> (std::istream& s,  Fastx& c);

bool operator<(Fastx& f1, Fastx& f2);

bool operator>(Fastx& f1, Fastx& f2);

#endif
