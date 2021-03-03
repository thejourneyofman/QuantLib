/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2021 Chase Xu

 This file is NOT an official part of QuantLib. It is for the purpose of test 
 instead of official version for QuantLib users.
 The details can be referred in README.md.
*/

/*! \file communicator.hpp
    \brief wrapper of Microsoft MPI interface
*/

#pragma once

#include <vector>
#include <numeric>
#include <assert.h>
#include <MSMPI/include/mpi.h>


struct COMM_FLOAT {
    unsigned int pid;
    unsigned int sid;
    double value;
    char name[10];
};

//A Wrapper for MPI to offer different communications methodology between processes and threads.

class Communicator {
public:
    Communicator(MPI_Comm mpi) : mpi_(mpi) {}

    ~Communicator()
    {
        MPI_Comm_free(&mpi_);
    }   

    MPI_Comm mpi() const { return mpi_; }

    void connect() {
        MPI_Comm_rank(mpi_, &rank_);
        MPI_Comm_size(mpi_, &size_);
    }

    /// Assignment operator
    const Communicator& operator() (const MPI_Comm& mpi)
    {
        this->mpi_ = mpi;
        return *this;
    }

    unsigned int rank() const { return rank_; }

    unsigned int size() const { return size_; }

protected:

    template<typename ValueType>
    void broadcast(std::vector<ValueType>& values, unsigned int broadcaster) const
    {
        values.resize(1);
        MPI_Bcast(const_cast<ValueType*>(values.data()), 1, mpi_type<ValueType>(), broadcaster,
                Communicator::mpi());
    }

    template<typename ValueType>
    void gather(const std::vector<ValueType>& in_values, std::vector<ValueType>& out_values,
        unsigned int receiving_process) const
    {
        const std::size_t comm_size = size();

        // Get data size on each process
        std::vector<int> pcounts(comm_size);
        const int local_size = in_values.size();
        MPI_Gather(const_cast<int*>(&local_size), 1, mpi_type<int>(),
            pcounts.data(), 1, mpi_type<int>(),
            receiving_process, Communicator::mpi());

        // Build offsets
        std::vector<int> offsets(comm_size + 1, 0);
        for (std::size_t i = 1; i <= comm_size; ++i)
            offsets[i] = offsets[i - 1] + pcounts[i - 1];

        const std::size_t n = std::accumulate(pcounts.begin(), pcounts.end(), 0);
        out_values.resize(n);
        MPI_Gatherv(const_cast<ValueType*>(in_values.data()), in_values.size(),
            mpi_type<ValueType>(),
            out_values.data(), pcounts.data(), offsets.data(),
            mpi_type<ValueType>(), receiving_process, Communicator::mpi());
    }

    //fixed the bug of original codes on unsupporting multiple threads for all_to_all mode
    template<typename ValueType>
    void all_to_all(std::vector<std::vector<ValueType>>& in_values,
        std::vector<ValueType>& unpacked_values) const
    {
        const std::size_t comm_size = Communicator::size();

        // Data size per destination
        assert(in_values.size() == comm_size);

        std::vector<int> recv_count(comm_size);
        std::vector<int> send_count(comm_size);
        std::vector<int> send_disp(comm_size + 1, 0);
        for (std::size_t p = 0; p < comm_size; ++p)
        {
            recv_count[p] = in_values[p].size();
            send_count[p] = in_values[p].size();
            send_disp[p + 1] = send_disp[p] + send_count[p];
        }

        std::vector<int> recv_disp(comm_size + 1, 0);
        std::vector<ValueType> send_buff(send_disp[comm_size]);
        for (std::size_t p = 0; p < comm_size; ++p)
        {
            recv_disp[p + 1] = recv_disp[p] + send_count[p];
            std::copy(in_values[p].begin(), in_values[p].end(),
                send_buff.begin() + send_disp[p]);
        }

        // Send/receive data
        unpacked_values.resize(recv_disp[comm_size]);

        MPI_Alltoallv(send_buff.data(), send_count.data(),
            send_disp.data(), mpi_type<ValueType>(),
            unpacked_values.data(), recv_count.data(), recv_disp.data(),
                      mpi_type<ValueType>(), Communicator::mpi());

        MPI_Barrier(Communicator::mpi());
    }

private:
    MPI_Comm mpi_;
    int rank_, size_;

    static inline MPI_Datatype comm_float = -1;

    template<typename ValueType>
    struct dependent_false : std::false_type {};
    template<typename ValueType> static int mpi_type()
    {        
        static_assert(dependent_false<ValueType>::value, "Unknown MPI Data Type");
        return MPI_CHAR;
    }

};

//template<> inline int Communicator::mpi_type<float>() { return MPI_FLOAT; }
//Rewrite MPI Datatypes to manage a Process ID and a thread ID
template<> inline int Communicator::mpi_type<COMM_FLOAT>() {
    if (comm_float == -1) {
        const int lens[4] = { 1, 1, 1, 10 };
        const MPI_Aint displace[4] = { 0, sizeof(unsigned int),
                                      sizeof(unsigned int) + sizeof(unsigned int),
                                      sizeof(unsigned int) + sizeof(unsigned int) + sizeof(double)};
        const MPI_Datatype types[4] = {MPI_UNSIGNED, MPI_UNSIGNED, MPI_DOUBLE, MPI_CHAR};
        MPI_Type_create_struct(4, lens, displace, types, &comm_float);
        MPI_Type_commit(&comm_float);
    }    
    return comm_float;
}



template<> inline int Communicator::mpi_type<double>() { return MPI_DOUBLE; }
template<> inline int Communicator::mpi_type<short int>() { return MPI_SHORT; }
template<> inline int Communicator::mpi_type<int>() { return MPI_INT; }
template<> inline int Communicator::mpi_type<long int>() { return MPI_LONG; }
template<> inline int Communicator::mpi_type<unsigned int>() { return MPI_UNSIGNED; }
template<> inline int Communicator::mpi_type<unsigned long int>() { return MPI_UNSIGNED_LONG; }
template<> inline int Communicator::mpi_type<long long>() { return MPI_LONG_LONG; }


