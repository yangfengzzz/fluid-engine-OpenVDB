//
//  globalClock.cpp
//  DigitalFlex
//
//  Created by 杨丰 on 2019/7/10.
//  Copyright © 2019 杨丰. All rights reserved.
//

#include "global_clock.hpp"

namespace vox {

template<> GlobalClock* Singleton<GlobalClock>::msSingleton = 0;
//-----------------------------------------------------------------------
GlobalClock* GlobalClock::getSingletonPtr(void)
{
    return msSingleton;
}
GlobalClock& GlobalClock::getSingleton(void)
{
    assert( msSingleton );  return ( *msSingleton );
}

GlobalClock::GlobalClock(){}
GlobalClock::~GlobalClock(){}


void GlobalClock::start(double m){
    t = m;
}

double GlobalClock::accumulate(double dt){
    t += dt;
    return t;
}

double GlobalClock::getTime(){
    return t;
}

}
