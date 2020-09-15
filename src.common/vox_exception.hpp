//
//  Exception.hpp
//  DigitalFlex
//
//  Created by 杨丰 on 2019/6/20.
//  Copyright © 2019 杨丰. All rights reserved.
//

#ifndef Exception_hpp
#define Exception_hpp

#include <exception>
#include <string>

// Check for OGRE assert mode
#ifndef OGRE_ASSERT_MODE
#define OGRE_ASSERT_MODE 1
#endif


// RELEASE_EXCEPTIONS mode
#if OGRE_ASSERT_MODE == 1
#   if OGRE_DEBUG_MODE
#       define OgreAssert( a, b ) assert( (a) && (b) )
#   else
#       if OGRE_COMP != OGRE_COMPILER_BORL
#           define OgreAssert( a, b ) if( !(a) ) OGRE_EXCEPT( Exception::ERR_RT_ASSERTION_FAILED, (b), "no function info")
#       else
#           define OgreAssert( a, b ) if( !(a) ) OGRE_EXCEPT( Exception::ERR_RT_ASSERTION_FAILED, (b), __func__ )
#       endif
#   endif

// EXCEPTIONS mode
#elif OGRE_ASSERT_MODE == 2
#   if OGRE_COMP != OGRE_COMPILER_BORL
#       define OgreAssert( a, b ) if( !(a) ) OGRE_EXCEPT( Exception::ERR_RT_ASSERTION_FAILED, (b), "no function info")
#   else
#       define OgreAssert( a, b ) if( !(a) ) OGRE_EXCEPT( Exception::ERR_RT_ASSERTION_FAILED, (b), __func__ )
#   endif

// STANDARD mode
#else
#   define OgreAssert( a, b ) assert( (a) && (b) )
#endif

#if OGRE_DEBUG_MODE
#   define OgreAssertDbg( a, b ) OgreAssert( a, b )
#else
#   define OgreAssertDbg( a, b )
#endif

namespace vox {
/** When thrown, provides information about an error that has occurred inside the engine.
 @remarks
 OGRE never uses return values to indicate errors. Instead, if an
 error occurs, an exception is thrown, and this is the object that
 encapsulates the detail of the problem. The application using
 OGRE should always ensure that the exceptions are caught, so all
 OGRE engine functions should occur within a
 try{} catch(Exception& e) {} block.
 @par
 The user application should never create any instances of this
 object unless it wishes to unify its error handling using the
 same object.
 */
class Exception : public std::exception
{
protected:
    long line;
    int number;
    std::string typeName;
    std::string description;
    std::string source;
    std::string file;
    mutable std::string fullDesc;
public:
    /** Static definitions of error codes.
     @todo
     Add many more exception codes, since we want the user to be able
     to catch most of them.
     */
    enum ExceptionCodes {
        ERR_CANNOT_WRITE_TO_FILE,
        ERR_INVALID_STATE,
        ERR_INVALIDPARAMS,
        ERR_RENDERINGAPI_ERROR,
        ERR_DUPLICATE_ITEM,
        ERR_ITEM_NOT_FOUND,
        ERR_FILE_NOT_FOUND,
        ERR_INTERNAL_ERROR,
        ERR_RT_ASSERTION_FAILED,
        ERR_NOT_IMPLEMENTED,
        ERR_INVALID_CALL
    };
    
    /** Default constructor.
     */
    Exception( int number, const std::string& description, const std::string& source );
    
    /** Advanced constructor.
     */
    Exception( int number, const std::string& description, const std::string& source, const char* type, const char* file, long line );
    
    /** Copy constructor.
     */
    Exception(const Exception& rhs);
    
    /// Needed for compatibility with std::exception
    ~Exception() throw() {}
    
    /** Assignment operator.
     */
    Exception & operator = (const Exception& rhs);
    
    /** Returns a string with the full description of this error.
     @remarks
     The description contains the error number, the description
     supplied by the thrower, what routine threw the exception,
     and will also supply extra platform-specific information
     where applicable. For example - in the case of a rendering
     library error, the description of the error will include both
     the place in which OGRE found the problem, and a text
     description from the 3D rendering library, if available.
     */
    virtual const std::string& getFullDescription(void) const;
    
    /** Gets the error code.
     */
    virtual int getNumber(void) const throw();
    
    /** Gets the source function.
     */
    virtual const std::string &getSource() const { return source; }
    
    /** Gets source file name.
     */
    virtual const std::string &getFile() const { return file; }
    
    /** Gets line number.
     */
    virtual long getLine() const { return line; }
    
    /** Returns a string with only the 'description' field of this exception. Use
     getFullDescriptionto get a full description of the error including line number,
     error number and what function threw the exception.
     */
    virtual const std::string &getDescription(void) const { return description; }
    
    /// Override std::exception::what
    const char* what() const throw() { return getFullDescription().c_str(); }
    
};


/** Template struct which creates a distinct type for each exception code.
 @note
 This is useful because it allows us to create an overloaded method
 for returning different exception types by value without ambiguity.
 From 'Modern C++ Design' (Alexandrescu 2001).
 */

// Specialised exceptions allowing each to be caught specifically
// backwards-compatible since exception codes still used

class UnimplementedException : public Exception
{
public:
    UnimplementedException(int inNumber, const std::string& inDescription, const std::string& inSource, const char* inFile, long inLine)
    : Exception(inNumber, inDescription, inSource, "UnimplementedException", inFile, inLine) {}
};
class FileNotFoundException : public Exception
{
public:
    FileNotFoundException(int inNumber, const std::string& inDescription, const std::string& inSource, const char* inFile, long inLine)
    : Exception(inNumber, inDescription, inSource, "FileNotFoundException", inFile, inLine) {}
};
class IOException : public Exception
{
public:
    IOException(int inNumber, const std::string& inDescription, const std::string& inSource, const char* inFile, long inLine)
    : Exception(inNumber, inDescription, inSource, "IOException", inFile, inLine) {}
};
class InvalidStateException : public Exception
{
public:
    InvalidStateException(int inNumber, const std::string& inDescription, const std::string& inSource, const char* inFile, long inLine)
    : Exception(inNumber, inDescription, inSource, "InvalidStateException", inFile, inLine) {}
};
class InvalidParametersException : public Exception
{
public:
    InvalidParametersException(int inNumber, const std::string& inDescription, const std::string& inSource, const char* inFile, long inLine)
    : Exception(inNumber, inDescription, inSource, "InvalidParametersException", inFile, inLine) {}
};
class ItemIdentityException : public Exception
{
public:
    ItemIdentityException(int inNumber, const std::string& inDescription, const std::string& inSource, const char* inFile, long inLine)
    : Exception(inNumber, inDescription, inSource, "ItemIdentityException", inFile, inLine) {}
};
class InternalErrorException : public Exception
{
public:
    InternalErrorException(int inNumber, const std::string& inDescription, const std::string& inSource, const char* inFile, long inLine)
    : Exception(inNumber, inDescription, inSource, "InternalErrorException", inFile, inLine) {}
};
class RenderingAPIException : public Exception
{
public:
    RenderingAPIException(int inNumber, const std::string& inDescription, const std::string& inSource, const char* inFile, long inLine)
    : Exception(inNumber, inDescription, inSource, "RenderingAPIException", inFile, inLine) {}
};
class RuntimeAssertionException : public Exception
{
public:
    RuntimeAssertionException(int inNumber, const std::string& inDescription, const std::string& inSource, const char* inFile, long inLine)
    : Exception(inNumber, inDescription, inSource, "RuntimeAssertionException", inFile, inLine) {}
};
class InvalidCallException : public Exception
{
public:
    InvalidCallException(int inNumber, const std::string& inDescription, const std::string& inSource, const char* inFile, long inLine)
    : Exception(inNumber, inDescription, inSource, "InvalidCallException", inFile, inLine) {}
};

/** Class implementing dispatch methods in order to construct by-value
 exceptions of a derived type based just on an exception code.
 @remarks
 This nicely handles construction of derived Exceptions by value (needed
 for throwing) without suffering from ambiguity - each code is turned into
 a distinct type so that methods can be overloaded. This allows OGRE_EXCEPT
 to stay small in implementation (desirable since it is embedded) whilst
 still performing rich code-to-type mapping.
 */
class ExceptionFactory
{
private:
    /// Private constructor, no construction
    ExceptionFactory() {}
public:
    static __attribute__((noreturn)) void throwException(
                                                         Exception::ExceptionCodes code, int number,
                                                         const std::string& desc,
                                                         const std::string& src, const char* file, long line)
    {
        switch (code)
        {
            case Exception::ERR_CANNOT_WRITE_TO_FILE:   throw IOException(number, desc, src, file, line);
            case Exception::ERR_INVALID_STATE:          throw InvalidStateException(number, desc, src, file, line);
            case Exception::ERR_INVALIDPARAMS:          throw InvalidParametersException(number, desc, src, file, line);
            case Exception::ERR_RENDERINGAPI_ERROR:     throw RenderingAPIException(number, desc, src, file, line);
            case Exception::ERR_DUPLICATE_ITEM:         throw ItemIdentityException(number, desc, src, file, line);
            case Exception::ERR_ITEM_NOT_FOUND:         throw ItemIdentityException(number, desc, src, file, line);
            case Exception::ERR_FILE_NOT_FOUND:         throw FileNotFoundException(number, desc, src, file, line);
            case Exception::ERR_INTERNAL_ERROR:         throw InternalErrorException(number, desc, src, file, line);
            case Exception::ERR_RT_ASSERTION_FAILED:    throw RuntimeAssertionException(number, desc, src, file, line);
            case Exception::ERR_NOT_IMPLEMENTED:        throw UnimplementedException(number, desc, src, file, line);
            case Exception::ERR_INVALID_CALL:           throw InvalidCallException(number, desc, src, file, line);
            default:                                    throw Exception(number, desc, src, "Exception", file, line);
        }
    }
    
};

}

#ifndef OGRE_EXCEPT
#define OGRE_EXCEPT(code, desc, src)         ExceptionFactory::throwException(code, code, desc, src, __FILE__, __LINE__)
#define OGRE_EXCEPT_EX(code, num, desc, src) ExceptionFactory::throwException(code, num, desc, src, __FILE__, __LINE__)
#else
#define OGRE_EXCEPT_EX(code, num, desc, src) OGRE_EXCEPT(code, desc, src)
#endif


#endif /* Exception_hpp */
