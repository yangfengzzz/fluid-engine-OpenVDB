//
//  globalClock.h
//  DigitalFlex
//
//  Created by 杨丰 on 2019/7/10.
//  Copyright © 2019 杨丰. All rights reserved.
//

#ifndef globalClock_h
#define globalClock_h

#include "singleton.h"

namespace vox {

class GlobalClock : public Singleton<GlobalClock> {
public:
    GlobalClock();
    ~GlobalClock();
    
    static GlobalClock& getSingleton(void);
    static GlobalClock* getSingletonPtr(void);
public:
    void start(double m = 0);
    
    double accumulate(double dt);
    
    double getTime();
    
private:
    double t;
};

}

#endif /* globalClock_h */
