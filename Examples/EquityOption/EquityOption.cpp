/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*!
 Copyright (C) 2005, 2006, 2007, 2009 StatPro Italia srl

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include <ql/qldefines.hpp>
#ifdef BOOST_MSVC
#  include <ql/auto_link.hpp>
#endif
#include <ql/instruments/vanillaoption.hpp>
#include <ql/pricingengines/vanilla/binomialengine.hpp>
#include <ql/pricingengines/vanilla/analyticeuropeanengine.hpp>
#include <ql/pricingengines/vanilla/analytichestonengine.hpp>
#include <ql/pricingengines/vanilla/baroneadesiwhaleyengine.hpp>
#include <ql/pricingengines/vanilla/bjerksundstenslandengine.hpp>
#include <ql/pricingengines/vanilla/batesengine.hpp>
#include <ql/pricingengines/vanilla/integralengine.hpp>
#include <ql/pricingengines/vanilla/fdblackscholesvanillaengine.hpp>
#include <ql/pricingengines/vanilla/mceuropeanengine.hpp>
#include <ql/pricingengines/vanilla/mcamericanengine.hpp>
#include <ql/pricingengines/vanilla/analyticeuropeanvasicekengine.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/utilities/dataformatters.hpp>
#include <ql/models/shortrate/onefactormodels/vasicek.hpp>

#include <iostream>
#include <iomanip>

using namespace QuantLib;

#if defined(QL_ENABLE_SESSIONS)
namespace QuantLib {

    ThreadKey sessionId() { return 0; }

}
#endif

#if defined(USE_MPI)
class OptionWrapper : public ObjectWrapper {
  public:
    OptionWrapper(void (Strategy<COMM_FLOAT>::*func)(const COMM_FLOAT&) const)
    : ObjectWrapper(func) {}
    OptionWrapper() : ObjectWrapper(&Strategy::gatherFromSlaves) {}

    void takeall() { 
        //Not implemented
    }

  protected:
    void notify(const std::vector<COMM_FLOAT>& unpackedvalues) const {
        for (std::size_t i = 0; i < unpackedvalues.size(); ++i) {
            std::cout << "[PID " << unpackedvalues[i].pid << ", SID " << unpackedvalues[i].sid
                      << ", " << unpackedvalues[i].name << "][" << getValueAttribute() << ":"
                      << unpackedvalues[i].value << "]" << std::endl;
        }
    }
};
#endif

int main(int, char* []) {

    try {        

        std::cout << std::endl;

        // set up dates
        Calendar calendar = TARGET();
        Date todaysDate(15, May, 1998);
        Date settlementDate(17, May, 1998);
        Settings::instance().evaluationDate() = todaysDate;

        // our options
        Option::Type type(Option::Put);
        Real underlying = 36;
        Real strike = 40;
        Spread dividendYield = 0.00;
        Rate riskFreeRate = 0.06;
        Volatility volatility = 0.20;
        Date maturity(17, May, 1999);
        DayCounter dayCounter = Actual365Fixed();

        std::cout << "Option type = "  << type << std::endl;
        std::cout << "Maturity = "        << maturity << std::endl;
        std::cout << "Underlying price = "        << underlying << std::endl;
        std::cout << "Strike = "                  << strike << std::endl;
        std::cout << "Risk-free interest rate = " << io::rate(riskFreeRate)
                  << std::endl;
        std::cout << "Dividend yield = " << io::rate(dividendYield)
                  << std::endl;
        std::cout << "Volatility = " << io::volatility(volatility)
                  << std::endl;
        std::cout << std::endl;
        std::string method;
        std::cout << std::endl ;

        // write column headings
        Size widths[] = { 35, 14, 14, 14 };
        std::cout << std::setw(widths[0]) << std::left << "Method"
                  << std::setw(widths[1]) << std::left << "European"
                  << std::setw(widths[2]) << std::left << "Bermudan"
                  << std::setw(widths[3]) << std::left << "American"
                  << std::endl;

        std::vector<Date> exerciseDates;
        for (Integer i=1; i<=4; i++)
            exerciseDates.push_back(settlementDate + 3*i*Months);

        ext::shared_ptr<Exercise> europeanExercise(
                                         new EuropeanExercise(maturity));

        ext::shared_ptr<Exercise> bermudanExercise(
                                         new BermudanExercise(exerciseDates));

        ext::shared_ptr<Exercise> americanExercise(
                                         new AmericanExercise(settlementDate,
                                                              maturity));

        Handle<Quote> underlyingH(
            ext::shared_ptr<Quote>(new SimpleQuote(underlying)));

        // bootstrap the yield/dividend/vol curves
        Handle<YieldTermStructure> flatTermStructure(
            ext::shared_ptr<YieldTermStructure>(
                new FlatForward(settlementDate, riskFreeRate, dayCounter)));
        Handle<YieldTermStructure> flatDividendTS(
            ext::shared_ptr<YieldTermStructure>(
                new FlatForward(settlementDate, dividendYield, dayCounter)));
        Handle<BlackVolTermStructure> flatVolTS(
            ext::shared_ptr<BlackVolTermStructure>(
                new BlackConstantVol(settlementDate, calendar, volatility,
                                     dayCounter)));
        ext::shared_ptr<StrikedTypePayoff> payoff(
                                        new PlainVanillaPayoff(type, strike));
        ext::shared_ptr<BlackScholesMertonProcess> bsmProcess(
                 new BlackScholesMertonProcess(underlyingH, flatDividendTS,
                                               flatTermStructure, flatVolTS));

        #if defined(USE_MPI)
        // options
        boost::shared_ptr<VanillaOption> europeanOption(
            new VanillaOption(payoff, europeanExercise));
        boost::shared_ptr<VanillaOption> bermudanOption(
            new VanillaOption(payoff, bermudanExercise));
        boost::shared_ptr<VanillaOption> americanOption(
            new VanillaOption(payoff, americanExercise));

        #else

        // options
        VanillaOption europeanOption(payoff, europeanExercise);
        VanillaOption bermudanOption(payoff, bermudanExercise);
        VanillaOption americanOption(payoff, americanExercise);

        #endif

        // Analytic formulas:

        // Black-Scholes for European
        method = "Black-Scholes";
        #if defined(USE_MPI)
        europeanOption->setPricingEngine(
            ext::shared_ptr<PricingEngine>(new AnalyticEuropeanEngine(bsmProcess)));
        auto portfolio = std::make_unique<OptionWrapper>();
        portfolio->reset(method);
        portfolio->subscribeSignal(europeanOption, "European");
        portfolio->start(); // Start pricing

        #else
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                                     new AnalyticEuropeanEngine(bsmProcess)));
        
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << "N/A"
                  << std::setw(widths[3]) << std::left << "N/A"
                  << std::endl;
        #endif

        //Vasicek rates model for European
        method = "Black Vasicek Model";
        Real r0 = riskFreeRate;
        Real a = 0.3;
        Real b = 0.3;
        Real sigma_r = 0.15;
        Real riskPremium = 0.0;
        Real correlation = 0.5;
        ext::shared_ptr<Vasicek> vasicekProcess(new Vasicek(r0, a, b, sigma_r, riskPremium));

        #if defined(USE_MPI)
        europeanOption->setPricingEngine(ext::shared_ptr<PricingEngine>(
            new AnalyticBlackVasicekEngine(bsmProcess, vasicekProcess, correlation)));
        portfolio->reset(method);
        portfolio->subscribeSignal(europeanOption, "European");
        portfolio->start(); // Start pricing

        #else 
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                new AnalyticBlackVasicekEngine(bsmProcess, vasicekProcess, correlation)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << "N/A"
                  << std::setw(widths[3]) << std::left << "N/A"
                  << std::endl;
        #endif

        // semi-analytic Heston for European
        method = "Heston semi-analytic";
        ext::shared_ptr<HestonProcess> hestonProcess(
            new HestonProcess(flatTermStructure, flatDividendTS,
                              underlyingH, volatility*volatility,
                              1.0, volatility*volatility, 0.001, 0.0));
        ext::shared_ptr<HestonModel> hestonModel(
                                              new HestonModel(hestonProcess));

        #if defined(USE_MPI)
        europeanOption->setPricingEngine(
            ext::shared_ptr<PricingEngine>(new AnalyticHestonEngine(hestonModel)));
        portfolio->reset(method);
        portfolio->subscribeSignal(europeanOption, "European");
        portfolio->start(); // Start pricing

        #else
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                                     new AnalyticHestonEngine(hestonModel)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << "N/A"
                  << std::setw(widths[3]) << std::left << "N/A"
                  << std::endl;
        #endif

        // semi-analytic Bates for European        
        method = "Bates semi-analytic";
        ext::shared_ptr<BatesProcess> batesProcess(
            new BatesProcess(flatTermStructure, flatDividendTS,
                             underlyingH, volatility*volatility,
                             1.0, volatility*volatility, 0.001, 0.0,
                             1e-14, 1e-14, 1e-14));
        ext::shared_ptr<BatesModel> batesModel(new BatesModel(batesProcess));

        #if defined(USE_MPI)
        europeanOption->setPricingEngine(
            ext::shared_ptr<PricingEngine>(new BatesEngine(batesModel)));
        portfolio->reset(method);
        portfolio->subscribeSignal(europeanOption, "European");
        portfolio->start(); // Start pricing

        #else
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                                                new BatesEngine(batesModel)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << "N/A"
                  << std::setw(widths[3]) << std::left << "N/A"
                  << std::endl;
        #endif

        // Barone-Adesi and Whaley approximation for American
        method = "Barone-Adesi/Whaley";
        #if defined(USE_MPI)
        americanOption->setPricingEngine(
            ext::shared_ptr<PricingEngine>(new BaroneAdesiWhaleyApproximationEngine(bsmProcess)));
        portfolio->reset(method);
        portfolio->subscribeSignal(americanOption, "American");
        portfolio->start(); // Start pricing
        #else
        americanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                       new BaroneAdesiWhaleyApproximationEngine(bsmProcess)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << "N/A"
                  << std::setw(widths[2]) << std::left << "N/A"
                  << std::setw(widths[3]) << std::left << americanOption.NPV()
                  << std::endl;
        #endif

        // Bjerksund and Stensland approximation for American
        method = "Bjerksund/Stensland";
        #if defined(USE_MPI)
        americanOption->setPricingEngine(
            ext::shared_ptr<PricingEngine>(new BjerksundStenslandApproximationEngine(bsmProcess)));
        portfolio->reset(method);
        portfolio->subscribeSignal(americanOption, "American");
        portfolio->start(); // Start pricing
        #else
        americanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                      new BjerksundStenslandApproximationEngine(bsmProcess)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << "N/A"
                  << std::setw(widths[2]) << std::left << "N/A"
                  << std::setw(widths[3]) << std::left << americanOption.NPV()
                  << std::endl;
        #endif

        // Integral
        method = "Integral";
        #if defined(USE_MPI)
        europeanOption->setPricingEngine(
            ext::shared_ptr<PricingEngine>(new IntegralEngine(bsmProcess)));
        portfolio->reset(method);
        portfolio->subscribeSignal(americanOption, "American");
        portfolio->start(); // Start pricing
        #else
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                                             new IntegralEngine(bsmProcess)));
        
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << "N/A"
                  << std::setw(widths[3]) << std::left << "N/A"
                  << std::endl;
        #endif

        // Finite differences
        Size timeSteps = 801;
        method = "Finite differences";

        #if defined(USE_MPI)
        ext::shared_ptr<BlackScholesMertonProcess> bsmProcess1(new BlackScholesMertonProcess(
            underlyingH,
            Handle<YieldTermStructure>(ext::shared_ptr<YieldTermStructure>(
                new FlatForward(settlementDate, dividendYield, dayCounter))),
            Handle<YieldTermStructure>(ext::shared_ptr<YieldTermStructure>(
                new FlatForward(settlementDate, riskFreeRate, dayCounter))),           
            flatVolTS));
        ext::shared_ptr<BlackScholesMertonProcess> bsmProcess2(new BlackScholesMertonProcess(
            underlyingH,
            Handle<YieldTermStructure>(ext::shared_ptr<YieldTermStructure>(
                new FlatForward(settlementDate, dividendYield, dayCounter))),
            Handle<YieldTermStructure>(ext::shared_ptr<YieldTermStructure>(
                new FlatForward(settlementDate, riskFreeRate, dayCounter))),
            flatVolTS));
        ext::shared_ptr<BlackScholesMertonProcess> bsmProcess3(new BlackScholesMertonProcess(
            underlyingH,
            Handle<YieldTermStructure>(ext::shared_ptr<YieldTermStructure>(
                new FlatForward(settlementDate, dividendYield, dayCounter))),
            Handle<YieldTermStructure>(ext::shared_ptr<YieldTermStructure>(
                new FlatForward(settlementDate, riskFreeRate, dayCounter))),
            flatVolTS));

        ext::shared_ptr<PricingEngine> fdengine1 =
            ext::make_shared<FdBlackScholesVanillaEngine>(bsmProcess1, timeSteps, timeSteps - 1);
        ext::shared_ptr<PricingEngine> fdengine2 =
            ext::make_shared<FdBlackScholesVanillaEngine>(bsmProcess2, timeSteps, timeSteps - 1);
        ext::shared_ptr<PricingEngine> fdengine3 =
            ext::make_shared<FdBlackScholesVanillaEngine>(bsmProcess3, timeSteps, timeSteps - 1);

        europeanOption->setPricingEngine(fdengine1);
        bermudanOption->setPricingEngine(fdengine2);
        americanOption->setPricingEngine(fdengine3);

        portfolio->reset(method);
        portfolio->subscribeSignal(europeanOption, "European");
        portfolio->subscribeSignal(bermudanOption, "Bermudan");
        portfolio->subscribeSignal(americanOption, "American");
        portfolio->start(); // Start pricing

        #else
        ext::shared_ptr<PricingEngine> fdengine =
            ext::make_shared<FdBlackScholesVanillaEngine>(bsmProcess,
                                                          timeSteps,
                                                          timeSteps-1);
        europeanOption.setPricingEngine(fdengine);
        bermudanOption.setPricingEngine(fdengine);
        americanOption.setPricingEngine(fdengine);
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << bermudanOption.NPV()
                  << std::setw(widths[3]) << std::left << americanOption.NPV()
                  << std::endl;
        #endif

        // Binomial method: Jarrow-Rudd
        method = "Binomial Jarrow-Rudd";
        #if defined(USE_MPI)
        europeanOption->setPricingEngine(ext::shared_ptr<PricingEngine>(
            new BinomialVanillaEngine<JarrowRudd>(bsmProcess1, timeSteps)));
        bermudanOption->setPricingEngine(ext::shared_ptr<PricingEngine>(
            new BinomialVanillaEngine<JarrowRudd>(bsmProcess2, timeSteps)));
        americanOption->setPricingEngine(ext::shared_ptr<PricingEngine>(
            new BinomialVanillaEngine<JarrowRudd>(bsmProcess3, timeSteps)));
        portfolio->reset(method);
        portfolio->subscribeSignal(europeanOption, "European");
        portfolio->subscribeSignal(bermudanOption, "Bermudan");
        portfolio->subscribeSignal(americanOption, "American");
        portfolio->start(); // Start pricing

        #else
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                new BinomialVanillaEngine<JarrowRudd>(bsmProcess,timeSteps)));
        bermudanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                new BinomialVanillaEngine<JarrowRudd>(bsmProcess,timeSteps)));
        americanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                new BinomialVanillaEngine<JarrowRudd>(bsmProcess,timeSteps)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << bermudanOption.NPV()
                  << std::setw(widths[3]) << std::left << americanOption.NPV()
                  << std::endl;

        #endif


        method = "Binomial Cox-Ross-Rubinstein";
        #if defined(USE_MPI)
        europeanOption->setPricingEngine(ext::shared_ptr<PricingEngine>(
            new BinomialVanillaEngine<CoxRossRubinstein>(bsmProcess1, timeSteps)));
        bermudanOption->setPricingEngine(ext::shared_ptr<PricingEngine>(
            new BinomialVanillaEngine<CoxRossRubinstein>(bsmProcess2, timeSteps)));
        americanOption->setPricingEngine(ext::shared_ptr<PricingEngine>(
            new BinomialVanillaEngine<CoxRossRubinstein>(bsmProcess3, timeSteps)));
        portfolio->reset(method);
        portfolio->subscribeSignal(europeanOption, "European");
        portfolio->subscribeSignal(bermudanOption, "Bermudan");
        portfolio->subscribeSignal(americanOption, "American");
        portfolio->start(); // Start pricing

        #else
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                      new BinomialVanillaEngine<CoxRossRubinstein>(bsmProcess,
                                                                   timeSteps)));
        bermudanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                      new BinomialVanillaEngine<CoxRossRubinstein>(bsmProcess,
                                                                   timeSteps)));
        americanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                      new BinomialVanillaEngine<CoxRossRubinstein>(bsmProcess,
                                                                   timeSteps)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << bermudanOption.NPV()
                  << std::setw(widths[3]) << std::left << americanOption.NPV()
                  << std::endl;
        #endif
        

        // Binomial method: Additive equiprobabilities
        method = "Additive equiprobabilities";
        #if defined(USE_MPI)
        europeanOption->setPricingEngine(ext::shared_ptr<PricingEngine>(
            new BinomialVanillaEngine<AdditiveEQPBinomialTree>(bsmProcess1, timeSteps)));
        bermudanOption->setPricingEngine(ext::shared_ptr<PricingEngine>(
            new BinomialVanillaEngine<AdditiveEQPBinomialTree>(bsmProcess2, timeSteps)));
        americanOption->setPricingEngine(ext::shared_ptr<PricingEngine>(
            new BinomialVanillaEngine<AdditiveEQPBinomialTree>(bsmProcess3, timeSteps)));
        portfolio->reset(method);
        portfolio->subscribeSignal(europeanOption, "European");
        portfolio->subscribeSignal(bermudanOption, "Bermudan");
        portfolio->subscribeSignal(americanOption, "American");
        portfolio->start(); // Start pricing

        #else
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                new BinomialVanillaEngine<AdditiveEQPBinomialTree>(bsmProcess,
                                                                   timeSteps)));
        bermudanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                new BinomialVanillaEngine<AdditiveEQPBinomialTree>(bsmProcess,
                                                                   timeSteps)));
        americanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                new BinomialVanillaEngine<AdditiveEQPBinomialTree>(bsmProcess,
                                                                   timeSteps)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << bermudanOption.NPV()
                  << std::setw(widths[3]) << std::left << americanOption.NPV()
                  << std::endl;
        #endif

        // Binomial method: Binomial Trigeorgis
        method = "Binomial Trigeorgis";
        #if defined(USE_MPI)
        europeanOption->setPricingEngine(ext::shared_ptr<PricingEngine>(
            new BinomialVanillaEngine<Trigeorgis>(bsmProcess, timeSteps)));
        bermudanOption->setPricingEngine(ext::shared_ptr<PricingEngine>(
            new BinomialVanillaEngine<Trigeorgis>(bsmProcess, timeSteps)));
        americanOption->setPricingEngine(ext::shared_ptr<PricingEngine>(
            new BinomialVanillaEngine<Trigeorgis>(bsmProcess, timeSteps)));
        portfolio->reset(method);
        portfolio->subscribeSignal(europeanOption, "European");
        portfolio->subscribeSignal(bermudanOption, "Bermudan");
        portfolio->subscribeSignal(americanOption, "American");
        portfolio->start(); // Start pricing
        #else
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                new BinomialVanillaEngine<Trigeorgis>(bsmProcess,timeSteps)));
        bermudanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                new BinomialVanillaEngine<Trigeorgis>(bsmProcess,timeSteps)));
        americanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                new BinomialVanillaEngine<Trigeorgis>(bsmProcess,timeSteps)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << bermudanOption.NPV()
                  << std::setw(widths[3]) << std::left << americanOption.NPV()
                  << std::endl;
        #endif

        // Binomial method: Binomial Tian
        method = "Binomial Tian";
        #if defined(USE_MPI)
        europeanOption->setPricingEngine(
            ext::shared_ptr<PricingEngine>(new BinomialVanillaEngine<Tian>(bsmProcess1, timeSteps)));
        bermudanOption->setPricingEngine(
            ext::shared_ptr<PricingEngine>(new BinomialVanillaEngine<Tian>(bsmProcess2, timeSteps)));
        americanOption->setPricingEngine(
            ext::shared_ptr<PricingEngine>(new BinomialVanillaEngine<Tian>(bsmProcess3, timeSteps)));
        portfolio->reset(method);
        portfolio->subscribeSignal(europeanOption, "European");
        portfolio->subscribeSignal(bermudanOption, "Bermudan");
        portfolio->subscribeSignal(americanOption, "American");
        portfolio->start(); // Start pricing
        #else
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                      new BinomialVanillaEngine<Tian>(bsmProcess,timeSteps)));
        bermudanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                      new BinomialVanillaEngine<Tian>(bsmProcess,timeSteps)));
        americanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                      new BinomialVanillaEngine<Tian>(bsmProcess,timeSteps)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << bermudanOption.NPV()
                  << std::setw(widths[3]) << std::left << americanOption.NPV()
                  << std::endl;
        #endif

        // Binomial method: Binomial Leisen-Reimer
        method = "Binomial Leisen-Reimer";
        #if defined(USE_MPI)
        europeanOption->setPricingEngine(ext::shared_ptr<PricingEngine>(
            new BinomialVanillaEngine<LeisenReimer>(bsmProcess1, timeSteps)));
        bermudanOption->setPricingEngine(ext::shared_ptr<PricingEngine>(
            new BinomialVanillaEngine<LeisenReimer>(bsmProcess2, timeSteps)));
        americanOption->setPricingEngine(ext::shared_ptr<PricingEngine>(
            new BinomialVanillaEngine<LeisenReimer>(bsmProcess3, timeSteps)));
        portfolio->reset(method);
        portfolio->subscribeSignal(europeanOption, "European");
        portfolio->subscribeSignal(bermudanOption, "Bermudan");
        portfolio->subscribeSignal(americanOption, "American");
        portfolio->start(); // Start pricing
        #else
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
              new BinomialVanillaEngine<LeisenReimer>(bsmProcess,timeSteps)));
        bermudanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
              new BinomialVanillaEngine<LeisenReimer>(bsmProcess,timeSteps)));
        americanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
              new BinomialVanillaEngine<LeisenReimer>(bsmProcess,timeSteps)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << bermudanOption.NPV()
                  << std::setw(widths[3]) << std::left << americanOption.NPV()
                  << std::endl;
        #endif
        // Binomial method: Binomial Joshi
        method = "Binomial Joshi";
        #if defined(USE_MPI)
        europeanOption->setPricingEngine(ext::shared_ptr<PricingEngine>(
            new BinomialVanillaEngine<Joshi4>(bsmProcess1, timeSteps)));
        bermudanOption->setPricingEngine(ext::shared_ptr<PricingEngine>(
            new BinomialVanillaEngine<Joshi4>(bsmProcess2, timeSteps)));
        americanOption->setPricingEngine(ext::shared_ptr<PricingEngine>(
            new BinomialVanillaEngine<Joshi4>(bsmProcess3, timeSteps)));
        portfolio->reset(method);
        portfolio->subscribeSignal(europeanOption, "European");
        portfolio->subscribeSignal(bermudanOption, "Bermudan");
        portfolio->subscribeSignal(americanOption, "American");
        portfolio->start(); // Start pricing
        #else
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                    new BinomialVanillaEngine<Joshi4>(bsmProcess,timeSteps)));
        bermudanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                    new BinomialVanillaEngine<Joshi4>(bsmProcess,timeSteps)));
        americanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                    new BinomialVanillaEngine<Joshi4>(bsmProcess,timeSteps)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << bermudanOption.NPV()
                  << std::setw(widths[3]) << std::left << americanOption.NPV()
                  << std::endl;
        #endif

        // Monte Carlo Method: MC (crude)
        timeSteps = 1;
        method = "MC (crude)";
        Size mcSeed = 42;
        ext::shared_ptr<PricingEngine> mcengine1;
        mcengine1 = MakeMCEuropeanEngine<PseudoRandom>(bsmProcess)
            .withSteps(timeSteps)
            .withAbsoluteTolerance(0.02)
            .withSeed(mcSeed);

        #if defined(USE_MPI)
        europeanOption->setPricingEngine(mcengine1);
        portfolio->reset(method);
        portfolio->subscribeSignal(europeanOption, "European");
        portfolio->start(); // Start pricing
        #else
        europeanOption.setPricingEngine(mcengine1);
        // Real errorEstimate = europeanOption.errorEstimate();
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << "N/A"
                  << std::setw(widths[3]) << std::left << "N/A"
                  << std::endl;
        #endif

        // Monte Carlo Method: QMC (Sobol)
        method = "QMC (Sobol)";
        Size nSamples = 32768;  // 2^15

        ext::shared_ptr<PricingEngine> mcengine2;
        mcengine2 = MakeMCEuropeanEngine<LowDiscrepancy>(bsmProcess)
            .withSteps(timeSteps)
            .withSamples(nSamples);
        #if defined(USE_MPI)
        europeanOption->setPricingEngine(mcengine2);
        portfolio->reset(method);
        portfolio->subscribeSignal(europeanOption, "European");
        portfolio->start(); // Start pricing
        #else
        europeanOption.setPricingEngine(mcengine2);
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << "N/A"
                  << std::setw(widths[3]) << std::left << "N/A"
                  << std::endl;
        #endif

        // Monte Carlo Method: MC (Longstaff Schwartz)
        method = "MC (Longstaff Schwartz)";
        ext::shared_ptr<PricingEngine> mcengine3;
        mcengine3 = MakeMCAmericanEngine<PseudoRandom>(bsmProcess)
            .withSteps(100)
            .withAntitheticVariate()
            .withCalibrationSamples(4096)
            .withAbsoluteTolerance(0.02)
            .withSeed(mcSeed);
        #if defined(USE_MPI)
        americanOption->setPricingEngine(mcengine3);
        portfolio->reset(method);
        portfolio->subscribeSignal(americanOption, "American");
        portfolio->start(); // Start pricing
        #else
        americanOption.setPricingEngine(mcengine3);
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << "N/A"
                  << std::setw(widths[2]) << std::left << "N/A"
                  << std::setw(widths[3]) << std::left << americanOption.NPV()
                  << std::endl;
        #endif
        // End test
        //MPI_Finalize();
        return 0;

    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "unknown error" << std::endl;
        return 1;
    }
}
