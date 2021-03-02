/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2021 Chase Xu

 This file is NOT an official part of QuantLib. It is for the purpose of test
 instead of official version for QuantLib users.
 The details can be referred in README.md.
*/

/*! \file strategy.hpp
    \Functions to support diffrent broadcasting mode via MPI.
*/

#pragma once
#include "Communicator.hpp"

template<typename ValueType>
class Strategy : public Communicator
{
public:
    Strategy(MPI_Comm mpi = MPI_COMM_WORLD, int argc = 0, char** argv = NULL) :
        Communicator(mpi) {
        int mpi_initialized;
        MPI_Initialized(&mpi_initialized);
        
        if (mpi_initialized) {
            connect();
            return;
        }
        int provided = -1;
        MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
        if (provided < MPI_THREAD_FUNNELED)
        {
            printf("The threading support level is lesser than that demanded.\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        connect();

    }

    ~Strategy() {
        MPI_Finalize();
    }

    void broadcastFromMaster(const ValueType&) const;

    void gatherFromSlaves(const ValueType&) const;

    void allToAll(const ValueType&) const;
    
    void CustomerStrategy(const ValueType&) const;

  protected:

    virtual void notify(const std::vector<ValueType>&) const = 0;

};

template<typename ValueType>
void Strategy<ValueType>::broadcastFromMaster(const ValueType& value) const {
    std::vector<ValueType> values;

    if (Communicator::rank() == 0)
    {
        values.push_back(value);
        Communicator::broadcast(values, 0);
    }
    if (Communicator::rank() > 0) {
        Communicator::broadcast(values, 0);
        notify(values);
    }
}

template<typename ValueType>
void Strategy<ValueType>::allToAll(const ValueType& value) const {
    const std::size_t comm_size = size();
    std::vector<std::vector<ValueType>> invalues(comm_size);
    for (unsigned int i = 0; i < comm_size; ++i) {
        invalues[i].push_back(value);
    }
    std::vector<ValueType> unpackedvalues;
    Communicator::all_to_all(invalues, unpackedvalues);

    notify(unpackedvalues);

}


template<typename ValueType>
void Strategy<ValueType>::gatherFromSlaves(const ValueType& value) const {
    std::vector<ValueType> invalues;
    if (Communicator::size() == 1) {
        std::cout << "Number of Processes for gather mode SHOULD be equal with or more than two." << std::endl;
        return;
    }
    if (Communicator::rank() == 0 && Communicator::size() > 1)
    {
        std::vector<ValueType> unpackedvalues;
        Communicator::gather(invalues, unpackedvalues, 0);
        notify(unpackedvalues);
    }
    if (Communicator::rank() > 0) {
        invalues.push_back(value);
        std::vector<ValueType> unpackedvalues;
        Communicator::gather(invalues, unpackedvalues, 0);
    }
}

