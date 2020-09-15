//
//  singleton.h
//  DigitalFlex
//
//  Created by 杨丰 on 2019/6/6.
//  Copyright © 2019 杨丰. All rights reserved.
//

#ifndef singleton_h
#define singleton_h

#include <cassert>

namespace vox {

/** Template class for creating single-instance global classes.
 */
template <typename T> class Singleton
{
private:
    /** @brief Explicit private copy constructor. This is a forbidden operation.*/
    Singleton(const Singleton<T> &);
    
    /** @brief Private operator= . This is a forbidden operation. */
    Singleton& operator=(const Singleton<T> &);
    
protected:
    
    static T* msSingleton;
    
public:
    Singleton( void )
    {
        assert( !msSingleton );
        msSingleton = static_cast< T* >( this );
    }
    ~Singleton( void )
    {  assert( msSingleton );  msSingleton = 0;  }
    static T& getSingleton( void )
    {   assert( msSingleton );  return ( *msSingleton ); }
    static T* getSingletonPtr( void )
    { return msSingleton; }
};

}
#endif /* singleton_h */
