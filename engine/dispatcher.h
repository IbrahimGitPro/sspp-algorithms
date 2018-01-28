/*
 *  Copyright (C) 2011 Universidad Simon Bolivar
 *  Modified by Ibrahim Abdoulahi Lastname.Firstname@gmail.com
 *
 * 
 *  Permission is hereby granted to distribute this software for
 *  non-commercial research purposes, provided that this copyright
 *  notice is included with any such distribution.
 *  
 *  THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 *  EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE
 *  SOFTWARE IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU
 *  ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.
 *  
 *  Blai Bonet, bonet@ldc.usb.ve
 *
 */

#ifndef DISPATCHER_H
#define DISPATCHER_H

#include "problem.h"
#include "hash.h"
#include "algorithm.h"
#include "parameters.h"
#include "rollout.h"
#include "parameters.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cassert>
#include <strings.h>

//#define DEBUG

namespace Dispatcher {
	
	inline const char *algorithm_name(int index) {
		switch( index ) {
        case  0: return "vi";
        case  1: return "slrtdp";
        case  2: return "ilao";
        case  3: return "hdp";
        case  4: return "lfds+";
        case  5: return "fvi";
		}
		return 0;
	}

template<typename T> struct algorithm_table_t {
    typedef size_t (*type)(const Problem::problem_t<T>&, const T&, Problem::hash_t<T>&, const Algorithm::parameters_t&);
    type operator[](int i) const {
        switch( i ) {
            case  0: return Algorithm::value_iteration<T>;
            case  1: return Algorithm::standard_lrtdp<T>;
            case  2: return Algorithm::improved_lao<T>;
            case  3: return Algorithm::hdp_driver<T>;
            case  4: return Algorithm::ldfs_plus_driver<T>;
            case  5: return Algorithm::fvi_driver<T>;
            default: return 0;
        }
        return 0;
    }
};

template<typename T> struct result_t {
    int algorithm_;
    const char *algorithm_name_;
    unsigned seed_;
    float value_;
    unsigned trials_;
    unsigned updates_;
    unsigned expansions_;
    unsigned psize_;
    float atime_;
    float htime_;
    unsigned numSccs_;
    Problem::hash_t<T> *hash_;
};

template<typename T>
inline void solve(const Problem::problem_t<T> &problem, const Heuristic::heuristic_t<T> *heuristic, const T &s, unsigned bitmap, const Algorithm::parameters_t &parameters, std::vector<result_t<T> > &results) {

    // solve problem with algorithms
    unsigned index = 0;
    while( bitmap != 0 ) {
        for(; bitmap % 2 == 0; ++index, bitmap = bitmap >> 1 );

        typename algorithm_table_t<T>::type algorithm = algorithm_table_t<T>()[index];
        if( algorithm != 0 ) {
            result_t<T> result;
            result.algorithm_ = index;
            result.algorithm_name_ = algorithm_name(index);
            result.seed_ = parameters.seed_;
            Random::set_seed(parameters.seed_);

            float start_time = Utils::read_time_in_seconds();
            result.hash_ = new Problem::hash_t<T>(problem, new Heuristic::wrapper_t<T>(heuristic));
            problem.clear_expansions();
            if( heuristic != 0 ) heuristic->reset_stats();
            result.trials_ = (*algorithm)(problem, s, *result.hash_, parameters);
            result.value_ = result.hash_->value(s);
            result.updates_ = result.hash_->updates();
            result.expansions_ = problem.expansions();
            result.psize_ = std::numeric_limits<unsigned>::max();
            if( algorithm != static_cast<typename algorithm_table_t<T>::type>(Algorithm::simple_astar<T>) )
                result.psize_ = problem.policy_size(*result.hash_, s);

            float end_time = Utils::read_time_in_seconds();
            result.htime_ = !heuristic ? 0 : heuristic->total_time();
            float dtime = !heuristic ? 0 : heuristic->eval_time();
            result.atime_ = end_time - start_time - dtime;
            result.numSccs_ = problem.tarjan(s, *result.hash_);
            results.push_back(result);
        }
        bitmap = bitmap >> 1;
        ++index;
    }
}

template<typename T>
inline void print_result(std::ostream &os, const result_t<T> *result) {
    os << std::fixed;
    if( result == 0 ) {
	   
		os << std::setw(4) << "#" << " "
		   << std::setw(7) << "alg" << " "
		   << std::setw(12) << "V*(s0)" << " "
		   << std::setw(12) << "trials" << " "
		   << std::setw(12) << "updates" << " "
		   << std::setw(12) << "expansions" << " "
		   << std::setw(12) << "hashsz" << " "
		   << std::setw(12) << "psize" << " "
		   << std::setw(12) << "atime" << " "
		   << std::setw(12) << "htime" << " "
           << std::setw(15) << "#policy-Sccs"
		   << std::endl;
    
    } else {
	    
		os << std::setw(4) << result->algorithm_ << " "
		   << std::setw(7) << result->algorithm_name_ << " "
		   << std::setw(12) << std::setprecision(10) << result->value_ << std::setprecision(2) << " "
		   << std::setw(12) << result->trials_ << " "
		   << std::setw(12) << result->updates_ << " "
		   << std::setw(12) << result->expansions_ << " "
		   << std::setw(12) << result->hash_->size() << " "
		   << std::setw(12) << result->psize_ << " "
		   << std::setw(12) << result->atime_ << " "
		   << std::setw(12) << result->htime_ << " "
            << std::setw(12) << result->numSccs_
		   << std::endl;
	    
    }
}

}; // namespace Dispatcher

#undef DEBUG

#endif

