/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2021 Chase Xu

 This file is NOT an official part of QuantLib. It inherits from LazyObject as
 to be a wrapper to support MPI and multiple threads in QuantLib usages and it
 is for the purpose of test instead of official version for QuantLib users.
 The details can be referred in README.md.
*/

/*! \file threaded_lazyobject.hpp
    \brief wrapper to support LazyObject in a parallelized style
*/

#ifndef quantlib_threaded_lazy_object_h
#define quantlib_threaded_lazy_object_h

#include <ql/patterns/lazyobject.hpp>
#include <boost/unordered_set.hpp>
#include <boost/bind/bind.hpp>
#include <boost/thread.hpp>
#include <boost/any.hpp>
#include <boost/atomic.hpp>
#include <boost/signals2.hpp>
#include <ql/parallel/Strategy.hpp>
#include <iostream>
#include <ql/pricingengine.hpp>


namespace QuantLib {

    using namespace boost::placeholders;
    class ObjectWrapper;
    //! Framework for calculation on demand and result caching.
    /*! \ingroup patterns */
    class ThreadedLazyObject : public virtual LazyObject
    {
        friend class ObjectWrapper;
      public:
          typedef boost::function<void(const COMM_FLOAT&)> MPIHandle;
          typedef boost::any any_type;
          using Runnable = boost::function<void()>;

          template <typename SignalType>
          class dispatcher : boost::noncopyable {
          public:
              explicit dispatcher(ThreadedLazyObject& sig) : sig_(sig), abortRequested_(false) {}

              ~dispatcher() {
                  mpi_.disconnect_all_slots();
                  if (thread_.joinable()) {
                      thread_.join();
                  }
              }

              template <typename F>
              void registerWith(const unsigned int& pid, const unsigned int& sid, const std::string name, F&& callback) {
                  pid_ = pid;
                  sid_ = sid;
                  name_ = name;
                  mpi_.connect(std::forward<F>(callback));
              }

              void post(ThreadedLazyObject::Runnable&& aRunnable) {
                  thread_ = boost::thread(std::move(aRunnable));
              }

              void stop() { 
                  abortAndJoin(); 
                  mpi_.disconnect_all_slots();
              }

              template <typename SignalType>
              void reset() {
                  additionalResults_.clear();
                  pid_ = 0;
                  sid_ = 0;
                  memset(&name_, 0, sizeof(std::string));
                  memset(&value_, 0, sizeof(SignalType));
              }

              template <typename SignalType>
              void notify(SignalType&& value) {
                  value.sid = sid_;
                  value.pid = pid_;
                  strncpy(value.name, name_.c_str(), 9);
                  value.name[9] = '\0';
                  mpi_(std::forward<SignalType>(value));
              }

          private:
              void abortAndJoin() {
                  abortRequested_.store(true);
                  if (thread_.joinable()) {
                      thread_.join();
                  }
              }

              ThreadedLazyObject& sig_;
              boost::thread thread_;
              boost::atomic_bool abortRequested_;
              boost::signals2::signal<void(SignalType)> mpi_;
              unsigned int pid_;
              unsigned int sid_;
              std::string name_;
              std::map<std::string, any_type> additionalResults_;
              SignalType value_;
          };

        
        ThreadedLazyObject();

        virtual ~ThreadedLazyObject() {}       
        ext::shared_ptr<dispatcher<COMM_FLOAT> > getSlave() const { return slave_; }
        ext::shared_ptr<PricingEngine> getEngine() const { return engine_; }        
        //@}
        //! \name Modifiers
        //@{
        //! set the pricing engine to be used.
        /*! \warning calling this method will have no effects in
                     case the <b>performCalculation</b> method
                     was overridden in a derived class.
        */
        void setPricingEngine(const ext::shared_ptr<PricingEngine>&);

      protected:
        virtual void calculate() const;
        void stop() const;
        virtual void performCalculations() const = 0;
        ext::shared_ptr<dispatcher<COMM_FLOAT> > slave_;
        ext::shared_ptr<PricingEngine> engine_;
    };

    class ObjectWrapper : public Strategy<COMM_FLOAT> {
    public:
        typedef boost::function<void(const COMM_FLOAT&)> MPIHandle;
        typedef boost::function<void()> ERRHandle;

        ObjectWrapper() {}

        ObjectWrapper(void (Strategy<COMM_FLOAT>::* func)(const COMM_FLOAT&) const) {
            hmpi_ = boost::bind(func, boost::ref(*this), _1);
            herr_ = boost::bind(&ObjectWrapper::errorHandle, boost::ref(*this));
        }

        ~ObjectWrapper() {
            stopSignals();
        }

        boost::unordered_set<boost::shared_ptr<ThreadedLazyObject> > signals() const {
            return Signals_;
        };

        size_t numofsignals() const { return Signals_.size(); };

        void start() {
            //setupSignals(e);
            startSignals();
            try {
                takeall();
            }
            catch (...) {
                stopSignals();
            }
            stopSignals();
        }

        auto subscribeSignal(const boost::shared_ptr<ThreadedLazyObject>& s, const std::string name) {
            if (s) {
                unsigned int processID = rank();
                unsigned int signalID = Signals_.size();
                s->getSlave()->registerWith(processID, signalID, name, hmpi_);
                return Signals_.insert(s);
            }
            return std::make_pair(Signals_.end(), false);
        }

        std::size_t unregisterSignal(const boost::shared_ptr<ThreadedLazyObject>& s) {
            return Signals_.erase(s);
        }

        void reset() { 
            valueAttribute_ = "value";
            Signals_.clear(); }

        void reset(const std::string& valueAttribute) { 
            valueAttribute_ = valueAttribute;
            Signals_.clear();
        }

        std::string getValueAttribute() const { return valueAttribute_; }

    private:
        boost::unordered_set<boost::shared_ptr<ThreadedLazyObject> > Signals_;
        MPIHandle hmpi_;
        ERRHandle herr_;
        std::string valueAttribute_;

    protected:
        virtual void takeall() = 0;

        void setupSignals(const ext::shared_ptr<PricingEngine>& e) {
            for (auto s : Signals_) 
            {
                s->setPricingEngine(e);
            }
        }
        void startSignals() {
            for (auto s : Signals_)
            {
                s->calculate();
            }
        }

        void stopSignals() {
            for (auto s : Signals_)
            {
                s->stop();
            }
        }

        void errorHandle() {
            throw "NOT implemented";
        }

    };

    // inline definitions

    inline ThreadedLazyObject::ThreadedLazyObject() : LazyObject(),
        slave_(ext::shared_ptr<dispatcher<COMM_FLOAT> >(new dispatcher<COMM_FLOAT>(*this))) {
    }    

    inline void ThreadedLazyObject::calculate() const {
        if (!calculated_ && !frozen_) {
            calculated_ = true;   // prevent infinite recursion in
                                  // case of bootstrapping
            try {
                slave_->post(boost::bind(&ThreadedLazyObject::performCalculations, boost::ref(*this)));
            } catch (...) {
                calculated_ = false;
                throw;
            }
        }
    }

    inline void ThreadedLazyObject::stop() const {
        try {
            slave_->stop();
        } catch (...) {
            calculated_ = false;
            throw;
        }
    }

    inline void ThreadedLazyObject::setPricingEngine(
        const ext::shared_ptr<PricingEngine>& e) {
        if (engine_ != 0)
            unregisterWith(engine_);
        engine_ = e;
        if (engine_ != 0)
            registerWith(engine_);
        // trigger (lazy) recalculation and notify observers
        update();
    }

}

#endif
