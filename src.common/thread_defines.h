//
//  ThreadDefines.h
//  DigitalFlex
//
//  Created by 杨丰 on 2019/6/20.
//  Copyright © 2019 杨丰. All rights reserved.
//

#ifndef ThreadDefines_h
#define ThreadDefines_h

#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>

#define AUTO_MUTEX_NAME mutex

#define TOKEN_PASTE(x, y) x ## y
#define TOKEN_PASTE_EXTRA(x, y) TOKEN_PASTE(x, y)

#define LOCK_AUTO_MUTEX std::unique_lock<std::recursive_mutex> AutoMutexLock(AUTO_MUTEX_NAME)
#define LOCK_MUTEX(name) std::unique_lock<std::recursive_mutex> TOKEN_PASTE_EXTRA(nameLock, __LINE__) (name)
#define LOCK_MUTEX_NAMED(mutexName, lockName) std::unique_lock<std::recursive_mutex> lockName(mutexName)
#define THREAD_SYNCHRONISER(sync) std::condition_variable_any sync
#define THREAD_SLEEP(ms) std::this_thread::sleep_for(std::chrono::milliseconds(ms))

#define AUTO_MUTEX mutable std::recursive_mutex AUTO_MUTEX_NAME
#define MUTEX(name) mutable std::recursive_mutex name
#define STATIC_MUTEX(name) static std::recursive_mutex name
#define STATIC_MUTEX_INSTANCE(name) std::recursive_mutex name
// like AUTO_MUTEX but mutex held by pointer
#define AUTO_SHARED_MUTEX mutable std::recursive_mutex *AUTO_MUTEX_NAME
#define LOCK_AUTO_SHARED_MUTEX assert(AUTO_MUTEX_NAME); std::recursive_mutex::scoped_lock AutoMutexLock(*AUTO_MUTEX_NAME)
#define NEW_AUTO_SHARED_MUTEX assert(!AUTO_MUTEX_NAME); AUTO_MUTEX_NAME = new std::recursive_mutex()
#define DELETE_AUTO_SHARED_MUTEX do { assert(AUTO_MUTEX_NAME); delete AUTO_MUTEX_NAME; } while (0)
#define COPY_AUTO_SHARED_MUTEX(from) assert(!AUTO_MUTEX_NAME); AUTO_MUTEX_NAME = from
#define SET_AUTO_SHARED_MUTEX_NULL AUTO_MUTEX_NAME = 0
#define MUTEX_CONDITIONAL(mutex) if (mutex)
#define THREAD_WAIT(sync, mutex, lock) sync.wait(lock)
#define THREAD_NOTIFY_ONE(sync) sync.notify_one()
#define THREAD_NOTIFY_ALL(sync) sync.notify_all()

// Read-write mutex
#define RW_MUTEX(name) mutable std::recursive_mutex name
#define LOCK_RW_MUTEX_READ(name) std::unique_lock<std::recursive_mutex> TOKEN_PASTE_EXTRA(nameLock, __LINE__) (name)
#define LOCK_RW_MUTEX_WRITE(name) std::unique_lock<std::recursive_mutex> TOKEN_PASTE_EXTRA(nameLock, __LINE__) (name)

// Thread objects and related functions
#define THREAD_TYPE std::thread
#define THREAD_CREATE(name, worker) std::thread* name = new std::thread(worker)
#define THREAD_DESTROY(name) delete name
#define THREAD_HARDWARE_CONCURRENCY std::thread::hardware_concurrency()
#define THREAD_CURRENT_ID std::this_thread::get_id()
#define THREAD_WORKER_INHERIT

// Utility
#define THREAD_ID_TYPE std::thread::id
#define THREAD_YIELD std::this_thread::yield()

#endif /* ThreadDefines_h */

