/*
 *  Copyright (C) 2011 Universidad Simon Bolivar
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

#ifndef PROBLEM_H
#define PROBLEM_H

#include "hash.h"
#include "random.h"
#include "utils.h"

#include <iostream>
#include <cassert>
#include <list>
#include <limits>
#include <vector>
#include <float.h>

//#define DEBUG

//static bool FLAG = false;

namespace Problem {

#ifndef __ACTION_TYPE
#define __ACTION_TYPE
    typedef int action_t;
    const action_t noop = -1;
#endif


struct labeling_flags {
	labeling_flags():newpolicy(false), solvable(true),unsolvable(true),deadend(true), residual(0), solvablescc(false){}
	bool newpolicy;
	bool solvable;
	bool unsolvable;
	bool deadend;
	float residual;
	bool solvablescc;
	//bool deadendscc;
};
	
template<typename T> class problem_t;

// The hash class implements a hash table that stores information related
// to the states of the problem which is used by different algorithms.

template<typename T> class hash_t : public Hash::hash_map_t<T> {

  public:
    typedef Hash::hash_map_t<T> base_type;

  protected:
    const problem_t<T> &problem_;
    unsigned updates_;

  public:
    hash_t(const problem_t<T> &problem, typename base_type::eval_function_t *heuristic = 0)
      : Hash::hash_map_t<T>(heuristic),
        problem_(problem), updates_(0) {
    }
    virtual ~hash_t() { }

    const problem_t<T>& problem() const { return problem_; }
    unsigned updates() const { return updates_; }
    void inc_updates() { ++updates_; }
    void update(const T &s, float value) {
        Hash::hash_map_t<T>::update(s, value);
        inc_updates();
    }

    virtual float QValue(const T &s, action_t a) const;
    std::pair<action_t, float> bestQValue(const T &s) const {
        action_t best_action = noop;
        float best_value = std::numeric_limits<float>::max();
        for( action_t a = 0; a < problem_.number_actions(s); ++a ) {
            if( problem_.applicable(s, a) ) {
                float value = QValue(s, a);
                if( value < best_value ) {
                    best_value = value;
                    best_action = a;
                }
            }
        }
        return std::make_pair(best_action, best_value);
   }

    //do backup and save best action
	float Backup(const T &s){
		Problem::action_t best_action = -1;
		std::vector<std::pair<T, float> > outcomes;
        float value, hv, qv, best_value = std::numeric_limits<float>::max();
		int a;
		unsigned osize;
		hv = this->value(s);
        for(a = 0; a < problem_.number_actions(s); ++a ) {
            if( problem_.applicable(s, a) ) {
				
				problem_.next(s, a, outcomes);
				osize = outcomes.size();
				qv = 0.0;
				for( unsigned i = 0; i < osize; ++i ) {
					qv += outcomes[i].second * this->value(outcomes[i].first);
				}
                value = problem_.cost(s, a) + problem_.discount() * qv;
				
                if( value < best_value ) {
                    best_value = value;
                    best_action = a;
                }
            }
        }
		update(s, best_value);
		this->set_action(s, best_action);

		
		float residual = best_value - hv;
		if(residual < 0)
			residual = -residual;
		return residual;
	}


    
    //do backup and save best action
    float SBackup(const T &s, bool &pichange){
        Problem::action_t best_action = -1;
        Problem::action_t old_action = this->action(s);
        std::vector<std::pair<T, float> > outcomes;
        float value, hv, qv, best_value = std::numeric_limits<float>::max();
        int a;
        unsigned osize;
        hv = this->value(s);
        for(a = 0; a < problem_.number_actions(s); ++a ) {
            if( problem_.applicable(s, a) ) {
                
                problem_.next(s, a, outcomes);
                osize = outcomes.size();
                qv = 0.0;
                for( unsigned i = 0; i < osize; ++i ) {
                    qv += outcomes[i].second * this->value(outcomes[i].first);
                }
                value = problem_.cost(s, a) + problem_.discount() * qv;
                
                if( value < best_value ) {
                    best_value = value;
                    best_action = a;
                }
            }
        }
        update(s, best_value);
        this->set_action(s, best_action);
        //check policy change
        if(old_action != best_action)
        pichange = true;
        
        float residual = best_value - hv;
        if(residual < 0)
        residual = -residual;
        return residual;
    }
    


 //do backup and save best action and check for unsolvability
	labeling_flags Backupprime(const T &s, Problem::hash_t<T> &hash, Hash::data_t* dptr){
		Problem::action_t old_best_action = this->action(s); //old best action
		Problem::action_t best_action = -1;
		std::vector<std::pair<T, float> > outcomes;
        float value, hv, qv, best_value = std::numeric_limits<float>::max();
		int a;
		unsigned osize;
		labeling_flags flags;

		//flags for checking unsolvability
		bool solv = true; //to check for solvability
		bool unsolv = true; //to check for unsolvability
		bool deadend = true ; //to check for deadend
		bool pchange = false; //to track policy change
		
		hv = this->value(s);
		bool unsolaction;
		bool solvablescc = false;//for solved, we need one action to decide
		//bool deadendscc = false;
		bool solveaction;
		bool deadendaction;
		bool oneissolved;
		bool oneisdeadend;
	  
        for(a = 0; a < problem_.number_actions(s); ++a ) { //for every action
            if( problem_.applicable(s, a) ) {
				solveaction = true;
				unsolaction = false;
				deadendaction = true;
				oneissolved = false;
				oneisdeadend = false;
				problem_.next(s, a, outcomes);
				osize = outcomes.size();
				qv = 0.0;
				for( unsigned i = 0; i < osize; ++i ) { //for every successor state of an action
					Hash::data_t *ptr = hash.data_ptr(outcomes[i].first);
					qv += outcomes[i].second * this->value(outcomes[i].first);

					if(ptr->solved())
						oneissolved = true;
					if(ptr->deadended())
						oneisdeadend = true;

					//check for solvable
					solveaction = solveaction && (ptr->solved() || (dptr->scc_low() == ptr->scc_low()));

					//check for unsol
					unsolaction = unsolaction || (ptr->unsolvabled() || (oneissolved && oneisdeadend));
					
					//check for deadend
					deadendaction = deadendaction && ((dptr->scc_low() == ptr->scc_low()) || ptr->deadended());
					
				}
				//one successor state is solved: See if this SCC contains a solved state
				if(oneissolved)
					solvablescc = true;
				//this action is a solving action, need only one for a state to be solvable
				if(solveaction)
					solv = true;
		    
				
				unsolv = unsolv && unsolaction;
				deadend = deadend && deadendaction;
                value = problem_.cost(s, a) + problem_.discount() * qv;
				
                if( value < best_value ) {
                    best_value = value;
                    best_action = a;
                }
            }
        }
	    
		
		update(s, best_value);
		this->set_action(s, best_action);
		if(old_best_action != best_action) // if flag then this is a postorder backup
		{
			pchange = true;
			
		}
		float residual = best_value - hv;
		if(residual < 0)
			residual = -residual;


	    
		flags.newpolicy = pchange;
		flags.solvable = solv;
		flags.unsolvable = unsolv;
		flags.deadend = deadend;
		flags.residual = residual;
		flags.solvablescc = solvablescc;
		//flags.deadendscc = deadendscc;
    
		return flags;
	}

	

};

template<typename T> class min_hash_t : public hash_t<T> {
    typedef Hash::hash_map_t<T> hash_base_type;

  public:
    min_hash_t(const problem_t<T> &problem,
               typename hash_base_type::eval_function_t *heuristic = 0)
      : hash_t<T>(problem, heuristic) {
    }
    virtual ~min_hash_t() { }
    virtual float QValue(const T &s, action_t a) const;
};


// A instance of problem_t represents an MDP problem. It contains all the 
// necessary information to run the different algorithms.

template<typename T> class problem_t {

  protected:
    float discount_;
    float dead_end_value_;
    mutable size_t expansions_;

  public:
    problem_t(float discount = 1.0, float dead_end_value = 1e6)
      : discount_(discount), dead_end_value_(dead_end_value), expansions_(0) { }
    virtual ~problem_t() { }

    float discount() const { return discount_; }
    float dead_end_value() const { return dead_end_value_; }

    size_t expansions() const {
        return expansions_;
    }
    void clear_expansions() const {
        expansions_ = 0;
    }

    virtual action_t number_actions(const T &s) const = 0;
    virtual const T& init() const = 0;
    virtual bool terminal(const T &s) const = 0;
    virtual bool dead_end(const T &s) const = 0;
    virtual bool applicable(const T &s, action_t a) const = 0;
    virtual float cost(const T &s, action_t a) const = 0;
    virtual void next(const T &s, action_t a, std::vector<std::pair<T, float> > &outcomes) const = 0;

    // sample next state given action using problem's dynamics
    std::pair<T, bool> sample(const T &s, action_t a) const {
        std::vector<std::pair<T, float> > outcomes;
        next(s, a, outcomes);
        unsigned osize = outcomes.size();
        assert(osize > 0);

        float r = Random::real();
        for( unsigned i = 0; i < osize; ++i ) {
            if( r < outcomes[i].second ) {
                return std::make_pair(outcomes[i].first, true);
            }
            r -= outcomes[i].second;
        }
        return std::make_pair(outcomes[0].first, true);
    }

    // sample next state given action uniformly among all possible next states
    std::pair<T, bool> usample(const T &s, action_t a) const {
        std::vector<std::pair<T, float> > outcomes;
        next(s, a, outcomes);
        unsigned osize = outcomes.size();
        return std::make_pair(outcomes[Random::uniform(osize)].first, true);
    }

    // sample next (unlabeled) state given action; probabilities are re-weighted
    std::pair<T, bool> nsample(const T &s, action_t a, const hash_t<T> &hash) const {
        std::vector<std::pair<T, float> > outcomes;
        next(s, a, outcomes);
        unsigned osize = outcomes.size();
        std::vector<bool> label(osize, false);

        size_t n = 0;
        float mass = 0;
        for( unsigned i = 0; i < osize; ++i ) {
            if( (label[i] = hash.solved(outcomes[i].first)) ) {
                mass += outcomes[i].second;
                ++n;
            }
        }

        mass = 1.0 - mass;
        n = osize - n;
        if( n == 0 ) return std::make_pair(s, false);

        float d = Random::real();
        for( unsigned i = 0; i < osize; ++i ) {
            if( !label[i] && ((n == 1) || (d <= outcomes[i].second / mass)) ) {
                return std::make_pair(outcomes[i].first, true);
            } else if( !label[i] ) {
                --n;
                d -= outcomes[i].second / mass;
            }
        }
        assert(0);
        return std::make_pair(s, false);
    }

    // type of sampling function
    typedef std::pair<T, bool> (problem_t<T>::*sample_function)(int) const;

    // compute the size of policy stored at hash table for given state
    size_t policy_size(hash_t<T> &hash, const T &s) const {
        hash.unmark_all();
        size_t size = policy_size_aux(hash, s);
        return size;
    }

    size_t policy_size_aux(hash_t<T> &hash, const T &s) const {
        std::vector<std::pair<T, float> > outcomes;
        size_t size = 0;
        if( !terminal(s) && !hash.marked(s) ) {
            hash.mark(s);
            std::pair<action_t, float> p = hash.bestQValue(s);
            next(s, p.first, outcomes);
            unsigned osize = outcomes.size();
            for( unsigned i = 0; i < osize; ++i ) {
                size += policy_size_aux(hash, outcomes[i].first);
            }
            ++size;
        }
        return size;
    }



   size_t tarjan( const T &s, hash_t<T> &hash) const {
        hash.unmark_all();
        size_t index = 0;
       size_t numSccs = 0;
	std::list<Hash::data_t*> stack;
        Hash::data_t *dptr = hash.data_ptr(s);
        tarjan_dfs(s, hash, dptr, index, stack, numSccs);
        return numSccs;
    }

   void tarjan_dfs(const T &s, hash_t<T> &hash, Hash::data_t* dptr, size_t &index, std::list<Hash::data_t*> &stack, size_t &numSccs) const {
	   size_t idx = index++;
	   dptr->set_scc_low(idx);
	   dptr->set_scc_idx(idx);
	   stack.push_front(dptr);
	   dptr->mark();
	   //BestAction successors:
	   if(dptr->action() >= 0){
	   std::vector<std::pair<T, float> > outcomes;
	   next(s, dptr->action(), outcomes);
	   unsigned osize = outcomes.size();
	   //Foreach sucessor state:
	   for( unsigned i = 0; i < osize; ++i ) {
           Hash::data_t *ptr = hash.data_ptr(outcomes[i].first);
		   if(ptr->scc_idx() == std::numeric_limits<unsigned>::max()) {
			   tarjan_dfs(outcomes[i].first, hash, ptr, index, stack, numSccs);
			   dptr->set_scc_low(Utils::min(dptr->scc_low(), ptr->scc_low()));
		   } else if( ptr->marked() ) {
			   dptr->set_scc_low(Utils::min(dptr->scc_low(), ptr->scc_idx()));
		   }
	   }
	   }
	   if( dptr->scc_low() == dptr->scc_idx() )
	   {
		   numSccs++; //found a new SCC
		   while( !stack.empty() && (stack.front() != dptr)/*->scc_idx() >= idx)*/ ) {
			   stack.pop_front();
		   } 
	   }
   }



    // print problem description
    virtual void print(std::ostream &os) const = 0;
};

template<typename T>
inline float hash_t<T>::QValue(const T &s, action_t a) const {
    if( problem_.terminal(s) ) return 0;

    std::vector<std::pair<T, float> > outcomes;
    problem_.next(s, a, outcomes);
    unsigned osize = outcomes.size();

    float qv = 0.0;
    for( unsigned i = 0; i < osize; ++i ) {
        qv += outcomes[i].second * this->value(outcomes[i].first);
    }
    return problem_.cost(s, a) + problem_.discount() * qv;
}

template<typename T>
inline float min_hash_t<T>::QValue(const T &s, action_t a) const {
    if( hash_t<T>::problem_.terminal(s) ) return 0;
    std::vector<std::pair<T, float> > outcomes;
    hash_t<T>::problem_.next(s, a, outcomes);
    unsigned osize = outcomes.size();

    float qv = std::numeric_limits<float>::max();
    for( unsigned i = 0; i < osize; ++i ) {
        qv = Utils::min(qv, this->value(outcomes[i].first));
    }
    return qv == std::numeric_limits<float>::max() ? std::numeric_limits<float>::max() : hash_t<T>::problem_.cost(s, a) + qv;
}

}; // namespace Problem

template<typename T>
inline std::ostream& operator<<(std::ostream &os, const Problem::problem_t<T> &problem) {
    problem.print(os);
    return os;
}

#undef DEBUG

#endif

