//
//  LogManager.cpp
//  src.common
//
//  Created by 杨丰 on 2020/6/5.
//  Copyright © 2020 杨丰. All rights reserved.
//

#include "log_manager.hpp"

#include "vox_exception.hpp"
#include <algorithm>
#include <iostream>

namespace vox {
//-----------------------------------------------------------------------
template<> LogManager* Singleton<LogManager>::msSingleton = 0;
LogManager* LogManager::getSingletonPtr(void)
{
    return msSingleton;
}
LogManager& LogManager::getSingleton(void)
{
    assert( msSingleton );  return ( *msSingleton );
}
//-----------------------------------------------------------------------
LogManager::LogManager()
{
    mDefaultLog = NULL;
}
//-----------------------------------------------------------------------
LogManager::~LogManager()
{
    LOCK_AUTO_MUTEX;
    // Destroy all logs
    LogList::iterator i;
    for (i = mLogs.begin(); i != mLogs.end(); ++i)
    {
        delete i->second;
    }
}
//-----------------------------------------------------------------------
Log* LogManager::createLog( const std::string& name, bool defaultLog, bool debuggerOutput,
                           bool suppressFileOutput)
{
    LOCK_AUTO_MUTEX;
    
    Log* newLog = new Log(name, debuggerOutput, suppressFileOutput);
    
    if( !mDefaultLog || defaultLog )
    {
        mDefaultLog = newLog;
    }
    
    mLogs.insert( LogList::value_type( name, newLog ) );
    
    return newLog;
}
//-----------------------------------------------------------------------
Log* LogManager::getDefaultLog()
{
    LOCK_AUTO_MUTEX;
    return mDefaultLog;
}
//-----------------------------------------------------------------------
Log* LogManager::setDefaultLog(Log* newLog)
{
    LOCK_AUTO_MUTEX;
    Log* oldLog = mDefaultLog;
    mDefaultLog = newLog;
    return oldLog;
}
//-----------------------------------------------------------------------
Log* LogManager::getLog( const std::string& name)
{
    LOCK_AUTO_MUTEX;
    LogList::iterator i = mLogs.find(name);
    if (i != mLogs.end())
        return i->second;
    else
        OGRE_EXCEPT(Exception::ERR_INVALIDPARAMS, "Log not found. ", "LogManager::getLog");
}

Log* LogManager::createOrRetrieve(const std::string& name, bool defaultLog, bool debuggerOutput,
                                  bool suppressFileOutput){
    LOCK_AUTO_MUTEX;
    LogList::iterator i = mLogs.find(name);
    if (i != mLogs.end())
        return i->second;
    else
        return createLog(name, defaultLog, debuggerOutput, suppressFileOutput);
}

//-----------------------------------------------------------------------
void LogManager::destroyLog(const std::string& name)
{
    LogList::iterator i = mLogs.find(name);
    if (i != mLogs.end())
    {
        if (mDefaultLog == i->second)
        {
            mDefaultLog = 0;
        }
        delete i->second;
        mLogs.erase(i);
    }
    
    // Set another default log if this one removed
    if (!mDefaultLog && !mLogs.empty())
    {
        mDefaultLog = mLogs.begin()->second;
    }
}
//-----------------------------------------------------------------------
void LogManager::destroyLog(Log* log)
{
    if(!log)
        OGRE_EXCEPT(Exception::ERR_INVALIDPARAMS, "Cannot destroy a null log.", "LogManager::destroyLog");
    
    destroyLog(log->getName());
}
//-----------------------------------------------------------------------
void LogManager::logMessage( const std::string& message, LogMessageLevel lml, bool maskDebug)
{
    LOCK_AUTO_MUTEX;
    if (mDefaultLog)
    {
        mDefaultLog->logMessage(message, lml, maskDebug);
    }
}
//-----------------------------------------------------------------------
void LogManager::setLogDetail(LoggingLevel ll)
{
    LOCK_AUTO_MUTEX;
    if (mDefaultLog)
    {
        mDefaultLog->setLogDetail(ll);
    }
}
//---------------------------------------------------------------------
Log::Stream LogManager::stream(LogMessageLevel lml, bool maskDebug)
{
    LOCK_AUTO_MUTEX;
    if (mDefaultLog)
        return mDefaultLog->stream(lml, maskDebug);
    else{
        mDefaultLog = new Log("default_log", true, false);
        return mDefaultLog->stream(lml, maskDebug);
    }
//        OGRE_EXCEPT(Exception::ERR_INVALIDPARAMS, "Default log not found. ", "LogManager::stream");
    
}

}
