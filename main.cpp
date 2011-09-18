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
#include "Fastx.h"
#include "gzstream.h"
#include <iostream>
#include <fstream>



int main(int argc, char * argv[])
{
  
    igzstream f(argv[1]);
    
    if(f.good())
    {
        do {
            Fastx ff;
            f >> ff;
            std::cout<<ff<<std::endl;
        } while (f.good());

        f.close();
    }
    else
    {
        std::cerr<<"error in opening file "<<argv[1]<<std::endl;
    }
  return 0;         
}