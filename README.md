# Parallelized Quantlib (For test use)

## Motivation
#### Challenged some ad-hoc way to make giant projects (in case of Quantlib, about half million LOC) work in a parallel computing style that supports both multiple processes and multiple threads by adding as little codes as possible.

## Restrictions
#### It's NOT an integral HPC version-up for Quantlib users but a test bed of the integrity and extensibility how a fulcrum can pry to a whole project. In this case, I added about 100-200 lines of code and started working from the base classes which are inherited by many top-level objects and functions.

## Features
#### For now it supports multiple threads in a particle size of Instrument, that is, one pricing engine (in some case it's bounded with a certain time-consuming stochastic process) works in an independent single thread for a single instrument.
#### Meanwhile, it provides a protocol to use MPI that can support multiple processes to manage the main thread and the instruments can be grouped to different processes in a in a particle size of "Portfolio" as to be a future work. 

## A use case (EquityOption)
#### You have a basket of options to be monitored and performed by calculations using different pricing engines for different option types. And each option may use different engines (BSM, Binomial, MC for an example) to compare the NPV or other indicators for the purpose of your strategy.

#### In its original synchronous mode, it works like

```europeanOption->setPricingEngine(bsEngine);
std::cout << europeanOption->NPV() << std::endl;
bermudanOption->setPricingEngine(binomialEngine);
std::cout << bermudanOption->NPV() << std::endl;
americanOption->setPricingEngine(mcEngine);
std::cout << americanOption->NPV() << std::endl;
```

In a new asynchronous mode, it works like

```europeanOption->setPricingEngine(bsEngine);
bermudanOption->setPricingEngine(binomialEngine);
americanOption->setPricingEngine(mcEngine);
portfolio->reset();
portfolio->subscribeSignal(europeanOption);
portfolio->subscribeSignal(bermudanOption);
portfolio->subscribeSignal(americanOption);
portfolio->start();
```

## Broadcast Strategies
#### It supports three kinds of broadcasting strategies to communicate all the signals from different processes and threads.
#### allToAll
All signals from each process and its threads will be broadcasted to all other processes (full connection).

#### gatherFromSlaves
Main processes will gather all signals from other processes and their threads.
(inward start connection)

#### broadcastFromMaster
Main processes will broadcast its signals to other processes and their threads.
(outward start connection)

## Objects/Classes Added
#### ThreadedLazyObject (based on LazyObject), which manages the thread, signal and slot.
#### ObjectWrapper, which handle the threading processes and errors.
#### Strategy, which provide different modes of broadcasting between processes.
#### Communicator, which is a simple wrapper of MPI.

## How to use it

#### go to the Solution Directory(where sln file is) \MSMPI\Bin
#### Run " .\mpiexec.exe -n [Number of Processes] Solution Directory\Examples\EquityOption\bin\EquityOption-x64-mt-gd.exe"

## How to build it
#### Following the guidance of official site of Quanlib with some additional configuration,
#### add the path of boost thread and MSMPI libraries to the project
#### add "USE_MPI" to the preprocessors
#### go to the file of userconfig.hpp, uncomment the following line
```ifndef QL_ENABLE_THREAD_SAFE_OBSERVER_PATTERN
#    define QL_ENABLE_THREAD_SAFE_OBSERVER_PATTERN
#endif
```

## Tools
#### -- Boost (http://www.boost.org)
#### -- MSMPI (https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi)

## Environment
#### -- Windows 10
#### -- Visual Studio 2019 (v142)
#### -- Microsoft .NET Framework (4.8.03752)
#### -- ISO C++17 Standard (/std:c++17)
#### -- boost (1.75.0)
#### -- MSMPI (10.1.2) 64 bit

## The following parts are cloned from original repository,
https://github.com/lballabio/QuantLib

# QuantLib: the free/open-source library for quantitative finance

[![Download](https://img.shields.io/github/v/release/lballabio/QuantLib?label=Download&sort=semver)](https://github.com/lballabio/QuantLib/releases/latest)
[![Licensed under the BSD 3-Clause License](https://img.shields.io/badge/License-BSD--3--Clause-blue.svg)](https://github.com/lballabio/QuantLib/blob/master/LICENSE.TXT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1440997.svg)](https://doi.org/10.5281/zenodo.1440997)
[![PRs Welcome](https://img.shields.io/badge/PRs%20-welcome-brightgreen.svg)](#contributing)

[![Linux build status](https://github.com/lballabio/QuantLib/workflows/Linux%20build/badge.svg?branch=master)](https://github.com/lballabio/QuantLib/actions?query=workflow%3A%22Linux+build%22)
[![Windows build status](https://ci.appveyor.com/api/projects/status/bmpiucu74eldfkm0/branch/master?svg=true)](https://ci.appveyor.com/project/lballabio/quantlib/branch/master)
[![Mac OS build status](https://github.com/lballabio/QuantLib/workflows/Mac%20OS%20build/badge.svg?branch=master)](https://github.com/lballabio/QuantLib/actions?query=workflow%3A%22Mac+OS+build%22)
[![CMake build status](https://github.com/lballabio/QuantLib/workflows/CMake%20build/badge.svg?branch=master)](https://github.com/lballabio/QuantLib/actions?query=workflow%3A%22CMake+build%22)

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/b4bc1058db994f24aa931b119a885eea)](https://www.codacy.com/app/lballabio/QuantLib)
[![Code Quality: Cpp](https://img.shields.io/lgtm/grade/cpp/g/lballabio/QuantLib.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/lballabio/QuantLib/context:cpp)
[![Coverage Status](https://coveralls.io/repos/github/lballabio/QuantLib/badge.svg?branch=master)](https://coveralls.io/github/lballabio/QuantLib?branch=master)

---

The QuantLib project (<http://quantlib.org>) is aimed at providing a
comprehensive software framework for quantitative finance. QuantLib is
a free/open-source library for modeling, trading, and risk management
in real-life.

QuantLib is Non-Copylefted Free Software and OSI Certified Open Source
Software.


## Download and usage

QuantLib can be downloaded from <http://quantlib.org/download.shtml>;
installation instructions are available at
<http://quantlib.org/install.shtml> for most platforms.

Documentation for the usage and the design of the QuantLib library is
available from <http://quantlib.org/docs.shtml>.

A list of changes for each past versions of the library can be
browsed at <http://quantlib.org/reference/history.html>.


## Questions and feedback

The preferred channel for questions (and the one with the largest
audience) is the quantlib-users mailing list.  Instructions for
subscribing are at <http://quantlib.org/mailinglists.shtml>.

Bugs can be reported as a GitHub issue at
<https://github.com/lballabio/QuantLib/issues>; if you have a patch
available, you can open a pull request instead (see "Contributing"
below).


## Contributing

The preferred way to contribute is through pull requests on GitHub.
Get a GitHub account if you don't have it already and clone the
repository at <https://github.com/lballabio/QuantLib> with the "Fork"
button in the top right corner of the page. Check out your clone to
your machine, code away, push your changes to your clone and submit a
pull request; instructions are available at
<https://help.github.com/articles/fork-a-repo>.

In case you need them, more detailed instructions for creating pull
requests are at
<https://help.github.com/articles/using-pull-requests>, and a basic
guide to GitHub is at
<https://guides.github.com/activities/hello-world/>.  GitHub also
provides interactive learning at <https://lab.github.com/>.

It's likely that we won't merge your code right away, and we'll ask
for some changes instead. Don't be discouraged! That's normal; the
library is complex, and thus it might take some time to become
familiar with it and to use it in an idiomatic way.

We're looking forward to your contributions.

