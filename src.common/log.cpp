//
//  Log.cpp
//  src.common
//
//  Created by 杨丰 on 2020/6/5.
//  Copyright © 2020 杨丰. All rights reserved.
//

#include "log.hpp"
#include <iomanip>
#include <iostream>
#include <algorithm>

namespace vox {
//-----------------------------------------------------------------------
Log::Log( const std::string& name, bool debuggerOuput, bool suppressFile ) :
mLogLevel(LL_NORMAL), mDebugOut(debuggerOuput),
mSuppressFile(suppressFile), mTimeStamp(true), mLogName(name)
{
    if (!mSuppressFile)
    {
        mLog.open(name.c_str());
    }
}
//-----------------------------------------------------------------------
Log::~Log()
{
    LOCK_AUTO_MUTEX;
    if (!mSuppressFile)
    {
        mLog.close();
    }
}
//-----------------------------------------------------------------------
void Log::logMessage( const std::string& message, LogMessageLevel lml, bool maskDebug )
{
    LOCK_AUTO_MUTEX;
    if ((mLogLevel + lml) >= OGRE_LOG_THRESHOLD)
    {
        bool skipThisMessage = false;
        for( mtLogListener::iterator i = mListeners.begin(); i != mListeners.end(); ++i )
            (*i)->messageLogged( message, lml, maskDebug, mLogName, skipThisMessage);
        
        if (!skipThisMessage)
        {

            if (mDebugOut && !maskDebug)
            {
                if (lml == LML_CRITICAL)
                    std::cerr << message << std::endl;
                else
                    std::cout << message << std::endl;
            }
            
            // Write time into log
            if (!mSuppressFile)
            {
                if (mTimeStamp)
                {
                    struct tm *pTime;
                    time_t ctTime; time(&ctTime);
                    pTime = localtime( &ctTime );
                    mLog << std::setw(2) << std::setfill('0') << pTime->tm_hour
                    << ":" << std::setw(2) << std::setfill('0') << pTime->tm_min
                    << ":" << std::setw(2) << std::setfill('0') << pTime->tm_sec
                    << ": ";
                }
                mLog << message << std::endl;
                
                // Flush stcmdream to ensure it is written (incase of a crash, we need log to be up to date)
                mLog.flush();
            }
        }
    }
}

//-----------------------------------------------------------------------
void Log::setTimeStampEnabled(bool timeStamp)
{
    LOCK_AUTO_MUTEX;
    mTimeStamp = timeStamp;
}

//-----------------------------------------------------------------------
void Log::setDebugOutputEnabled(bool debugOutput)
{
    LOCK_AUTO_MUTEX;
    mDebugOut = debugOutput;
}

//-----------------------------------------------------------------------
void Log::setLogDetail(LoggingLevel ll)
{
    LOCK_AUTO_MUTEX;
    mLogLevel = ll;
}

//-----------------------------------------------------------------------
void Log::addListener(LogListener* listener)
{
    LOCK_AUTO_MUTEX;
    mListeners.push_back(listener);
}

//-----------------------------------------------------------------------
void Log::removeListener(LogListener* listener)
{
    LOCK_AUTO_MUTEX;
    mListeners.erase(std::find(mListeners.begin(), mListeners.end(), listener));
}
//---------------------------------------------------------------------
Log::Stream Log::stream(LogMessageLevel lml, bool maskDebug)
{
    return Stream(this, lml, maskDebug);
}

}
