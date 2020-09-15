//
//  Exception.cpp
//  DigitalFlex
//
//  Created by 杨丰 on 2019/6/20.
//  Copyright © 2019 杨丰. All rights reserved.
//

#include "vox_exception.hpp"
#include "pch.h"

namespace vox {

Exception::Exception(int num, const std::string& desc, const std::string& src) :
line( 0 ),
number( num ),
description( desc ),
source( src )
{
    // Log this error - not any more, allow catchers to do it
    //LogManager::getSingleton().logMessage(this->getFullDescription());
}

Exception::Exception(int num, const std::string& desc, const std::string& src,
                     const char* typ, const char* fil, long lin) :
line( lin ),
number( num ),
typeName(typ),
description( desc ),
source( src ),
file( fil )
{
    // Log this error, mask it from debug though since it may be caught and ignored
    LOG(ERROR) << this->getFullDescription();
}

Exception::Exception(const Exception& rhs)
: line( rhs.line ),
number( rhs.number ),
typeName( rhs.typeName ),
description( rhs.description ),
source( rhs.source ),
file( rhs.file )
{
}

Exception& Exception::operator = ( const Exception& rhs )
{
    description = rhs.description;
    number = rhs.number;
    source = rhs.source;
    file = rhs.file;
    line = rhs.line;
    typeName = rhs.typeName;
    
    return *this;
}

const std::string& Exception::getFullDescription(void) const
{
    if (fullDesc.empty())
    {
        
        std::stringstream desc;
        
        desc <<  "OGRE EXCEPTION(" << number << ":" << typeName << "): "
        << description
        << " in " << source;
        
        if( line > 0 )
        {
            desc << " at " << file << " (line " << line << ")";
        }
        
        fullDesc = desc.str();
    }
    
    return fullDesc;
}

int Exception::getNumber(void) const throw()
{
    return number;
}

}
