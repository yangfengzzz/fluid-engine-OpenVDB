//
//  LogManager.hpp
//  src.common
//
//  Created by 杨丰 on 2020/6/5.
//  Copyright © 2020 杨丰. All rights reserved.
//

#ifndef LogManager_hpp
#define LogManager_hpp

#include <map>
#include "log.hpp"
#include "singleton.h"

namespace vox {
/** The log manager handles the creation and retrieval of logs for the
 application.
 @remarks
 This class will create new log files and will retrieve instances
 of existing ones. Other classes wishing to log output can either
 create a fresh log or retrieve an existing one to output to.
 One log is the default log, and is the one written to when the
 logging methods of this class are called.
 @par
 By default, Root will instantiate a LogManager (which becomes the
 Singleton instance) on construction, and will create a default log
 based on the Root construction parameters. If you want more control,
 for example redirecting log output right from the start or suppressing
 debug output, you need to create a LogManager yourself before creating
 a Root instance, then create a default log. Root will detect that
 you've created one yourself and won't create one of its own, thus
 using all your logging preferences from the first instance.
 */
class LogManager : public Singleton<LogManager>
{
protected:
    typedef std::map<std::string, Log*> LogList;
    
    /// A list of all the logs the manager can access
    LogList mLogs;
    
    /// The default log to which output is done
    Log* mDefaultLog;
    
public:
    AUTO_MUTEX; // public to allow external locking
    
    LogManager();
    ~LogManager();
    
    /** Creates a new log with the given name.
     @param
     name The name to give the log e.g. 'Ogre.log'
     @param
     defaultLog If true, this is the default log output will be
     sent to if the generic logging methods on this class are
     used. The first log created is always the default log unless
     this parameter is set.
     @param
     debuggerOutput If true, output to this log will also be
     routed to the debugger's output window.
     @param
     suppressFileOutput If true, this is a logical rather than a physical
     log and no file output will be written. If you do this you should
     register a LogListener so log output is not lost.
     */
    Log* createLog( const std::string& name, bool defaultLog = false, bool debuggerOutput = true,
                   bool suppressFileOutput = false);
    
    /** Retrieves a log managed by this class.
     */
    Log* getLog( const std::string& name);
    
    /** Retrieves or Create a log managed by this class.
     */
    Log* createOrRetrieve(const std::string& name, bool defaultLog = false, bool debuggerOutput = true,
                          bool suppressFileOutput = false);
    
    /** Returns a pointer to the default log.
     */
    Log* getDefaultLog();
    
    /** Closes and removes a named log. */
    void destroyLog(const std::string& name);
    /** Closes and removes a log. */
    void destroyLog(Log* log);
    
    /** Sets the passed in log as the default log.
     @return The previous default log.
     */
    Log* setDefaultLog(Log* newLog);
    
    /** Log a message to the default log.
     */
    void logMessage( const std::string& message, LogMessageLevel lml = LML_NORMAL,
                    bool maskDebug = false);
    
    /** Log a message to the default log (signature for backward compatibility).
     */
    void logMessage( LogMessageLevel lml, const std::string& message,
                    bool maskDebug = false) { logMessage(message, lml, maskDebug); }
    
    /** Get a stream on the default log. */
    Log::Stream stream(LogMessageLevel lml = LML_NORMAL,
                       bool maskDebug = false);
    
    /** Sets the level of detail of the default log.
     */
    void setLogDetail(LoggingLevel ll);
    /** Override standard Singleton retrieval.
     @remarks
     Why do we do this? Well, it's because the Singleton
     implementation is in a .h file, which means it gets compiled
     into anybody who includes it. This is needed for the
     Singleton template to work, but we actually only want it
     compiled into the implementation of the class based on the
     Singleton, not all of them. If we don't change this, we get
     link errors when trying to use the Singleton-based class from
     an outside dll.
     @par
     This method just delegates to the template version anyway,
     but the implementation stays in this single compilation unit,
     preventing link errors.
     */
    static LogManager& getSingleton(void);
    /** Override standard Singleton retrieval.
     @remarks
     Why do we do this? Well, it's because the Singleton
     implementation is in a .h file, which means it gets compiled
     into anybody who includes it. This is needed for the
     Singleton template to work, but we actually only want it
     compiled into the implementation of the class based on the
     Singleton, not all of them. If we don't change this, we get
     link errors when trying to use the Singleton-based class from
     an outside dll.
     @par
     This method just delegates to the template version anyway,
     but the implementation stays in this single compilation unit,
     preventing link errors.
     */
    static LogManager* getSingletonPtr(void);
};

}

#endif /* LogManager_hpp */
