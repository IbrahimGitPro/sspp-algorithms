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

#ifndef ALGORITHM_H
#define ALGORITHM_H

#include "problem.h"
#include "parameters.h"
#include<iomanip>
#include <cassert>
#include <list>
#include <queue>
#include <vector>
#include <math.h>





namespace Algorithm {

template<typename T>
struct min_priority_t {
    bool operator()(const std::pair<T, Hash::data_t*> &p1, const std::pair<T, Hash::data_t*> &p2) const {
        return p1.second->f() > p2.second->f();
    }
};



	

template<typename T>
size_t simple_astar(const Problem::problem_t<T> &problem,
                    const T &s,
                    Problem::hash_t<T> &hash,
                    const parameters_t &parameters) {

    typedef typename std::priority_queue<std::pair<T, Hash::data_t*>,
                                         std::vector<std::pair<T, Hash::data_t*> >,
                                         min_priority_t<T> >
            priority_queue;

    priority_queue open;
    std::vector<std::pair<T, float> > outcomes;

    Hash::data_t *dptr = hash.data_ptr(s);
    dptr->set_g(0);
    dptr->set_parent(0);
    dptr->set_action(Problem::noop);
    dptr->mark();
    open.push(std::make_pair(s, dptr));

#ifdef DEBUG
    std::cout << "PUSH " << "state" ///s
              << " w/ g=" << dptr->g()
              << " and h=" << dptr->h()
              << " => f=" << dptr->f()
              << std::endl;
    std::cout << "queue.sz=" << open.size() << std::endl;
#endif

    while( !open.empty() ) {
        std::pair<T, Hash::data_t*> n = open.top();
        open.pop();

#ifdef DEBUG
        std::cout << "POP  " << "state: " << n.first
                  << " w/ g=" << n.second->g()
                  << " and h=" << n.second->h()
                  << " => f=" << n.second->f()
                  << std::endl;
        std::cout << "queue.sz=" << open.size() << std::endl;

        std::vector<const Hash::data_t*> path;
        path.reserve(n.second->g());
        for( const Hash::data_t *dptr = n.second; dptr != 0; dptr = dptr->parent() )
            path.push_back(dptr);

        std::cout << "path=<";
        for( size_t i = path.size() - 1; i > 0; --i )
            std::cout << path[i-1]->action() << ",";
        std::cout << ">" << std::endl;
#endif

        // check for termination
        if( problem.terminal(n.first) ) {
#ifdef DEBUG
            std::cout << "GOAL FOUND!" << std::endl;
#endif

            std::vector<const Hash::data_t*> plan;
            plan.reserve(n.second->g());
            for( const Hash::data_t *dptr = n.second; dptr != 0; dptr = dptr->parent() )
                plan.push_back(dptr);

            std::cout << "plan=<";
            for( size_t i = plan.size() - 1; i > 0; --i )
                std::cout << plan[i-1]->action() << ",";
            std::cout << ">:" << plan.size() - 1 << std::endl;

            return (int)n.second->g();
        }

        // expand state
        for( Problem::action_t a = 0; a < problem.number_actions(n.first); ++a ) {
            if( problem.applicable(n.first, a) ) {
                problem.next(n.first, a, outcomes);
                assert(outcomes.size() == 1);
                Hash::data_t *ptr = hash.data_ptr(outcomes[0].first);
                float g = n.second->g() + problem.cost(n.first, a);
                if( !ptr->marked() || (g + ptr->h() < ptr->f()) ) {
                    ptr->set_g(g);
                    ptr->set_parent(n.second);
                    ptr->set_action(a);
                    ptr->mark();
                    open.push(std::make_pair(outcomes[0].first, ptr));
#ifdef DEBUG
                    std::cout << "PUSH " << "state: " <<outcomes[0].first
                              << " w/ g=" << ptr->g()
                              << " and h=" << ptr->h()
                              << " => f=" << ptr->f()
                              << std::endl;
                    std::cout << "queue.sz=" << open.size() << std::endl;
#endif
                }
            }
        }
    }
    return std::numeric_limits<size_t>::max();
}


////////////////////////////Code by Ibrahim ///////////////////////////////////////////////////

template<typename T>
class Hansen_LAO {
protected:
    const Problem::problem_t<T> &problem_;
    Problem::hash_t<T> &hash_;
	const parameters_t &parameters;
	
public:
Hansen_LAO(const Problem::problem_t<T> &problem, Problem::hash_t<T> &hash, const parameters_t &para)
	: problem_(problem), hash_(hash),parameters(para), NumExpandedStatesIter(0), NumExpandedStates(0),
		NumSolutionStates(0),ExitFlag(false), Residual(0.0) {}
    ~Hansen_LAO() { }
	size_t Iteration;
	size_t NumExpandedStatesIter; 
	size_t NumExpandedStates; 
	size_t NumSolutionStates;
	bool ExitFlag;
	float Residual;

	size_t iters(void){
		return Iteration;
	}

	//called on non terminal states
	float Backup(const T &s){
		Hash::data_t *dptr0, *dptr = hash_.data_ptr(s);
		Problem::action_t best_action = -1;
		std::vector<std::pair<T, float> > outcomes;
        float value, hv, qv, best_value = std::numeric_limits<float>::max();
		int a;
		unsigned osize;
		hv = dptr->value();
        for(a = 0; a < problem_.number_actions(s); ++a ) {
            if( problem_.applicable(s, a) ) {

				problem_.next(s, a, outcomes);
				osize = outcomes.size();
				qv = 0.0;
				for( unsigned i = 0; i < osize; ++i ) {
					dptr0 = hash_.data_ptr(outcomes[i].first);
					qv += outcomes[i].second * dptr0->value();
				}
                value = problem_.cost(s, a) + problem_.discount() * qv;
				
                if( value < best_value ) {
                    best_value = value;
                    best_action = a;
                }
            }
        }
		dptr->update(best_value);
		dptr->set_action(best_action);
		return (float)fabs(best_value - hv);
	}
	
	void ExpandSolution(const T &s) {
		
		std::vector<std::pair<T, float> > outcomes;
		std::pair<Problem::action_t, float> p;
		Hash::data_t *dptr = hash_.data_ptr(s);
		
		if(problem_.terminal(s)){
			dptr->update(0.0);
			return;
		}
		
		if(dptr->ExpandValue() == 0){
			dptr->SetExpand(Iteration);
			NumExpandedStatesIter++;
			NumExpandedStates++;
		}

		float Diff = Backup(s);
		if(Diff > Residual)
			Residual = Diff;
    		
		problem_.next(s, dptr->action(), outcomes);
		unsigned osize = outcomes.size();
		for( unsigned i = 0; i < osize; ++i ) {
			dptr = hash_.data_ptr(outcomes[i].first);
			// Stop if already visited this iteration or terminal state
			if((dptr->ExpandValue() < Iteration) && (!problem_.terminal(outcomes[i].first))) {
				// If already expanded, just mark as visited this iteration 
				if (dptr->ExpandValue() > 0)
					dptr->SetExpand(Iteration);
				ExpandSolution(outcomes[i].first);
			}
		}
		// Perform backup when backtrack to this node 
	    Backup(s);
	}

	size_t ConvergenceTest(const T &Start)
	{
		ExitFlag = 0;
		do {
			Iteration++;
			NumSolutionStates = 0;
			Residual = 0.0;
		   	ConvergenceTestRecursive(Start);  
	        
			
		} while (!ExitFlag && (Residual > parameters.epsilon_));
		if (ExitFlag){
			return( 0 );
		}
		else{
			return (1);
		}
	}

	void ConvergenceTestRecursive(const T &s)
	{
		std::vector<std::pair<T, float> > outcomes;
		std::pair<Problem::action_t, float> p;
		Hash::data_t *dptr = hash_.data_ptr(s);
		
		// Count this node as part of best solution 
		NumSolutionStates++;
		
		// Terminate recursion at goal node
		if(problem_.terminal(s)){
			dptr->update(0.0);
			return;
		}
		// Exit convergence test if best solution has unexpanded node
		if(dptr->ExpandValue() == 0)
			ExitFlag = true;

		
		//Backup
		float Diff = Backup(s);
		if(Diff > Residual)
			Residual = Diff;

		problem_.next(s, dptr->action(), outcomes);
		// Recursively consider unvisited successor nodes
		unsigned osize = outcomes.size();
		for( unsigned i = 0; i < osize; ++i ) {
			dptr = hash_.data_ptr(outcomes[i].first);
			
			if (dptr->Update() < Iteration){
				dptr->SetUpdate(Iteration);
				ConvergenceTestRecursive(outcomes[i].first);
			}
		}
	}


	void lLAOstar_driver(const T &s)
	{

		NumExpandedStates = 0;
		NumExpandedStatesIter = 0;
	    
		for (Iteration = 1; ;Iteration++) {
    
  
			// Expand best partial solution
			NumExpandedStatesIter = 0;
			ExpandSolution(s);
    
			// Call convergence test if solution graph has no unexpanded tip states
			if ((NumExpandedStatesIter == 0) && ConvergenceTest(s))
				break;
			
		}
	}

};
/////////////////////////////////End code LAO by Ibrahim////////////////////////




template<typename T>
class policy_graph_t {
  protected:
    typedef typename std::list<std::pair<T, Hash::data_t*> > pair_list;
    typedef typename std::vector<T>::iterator vector_iterator;
    typedef typename pair_list::iterator list_iterator;
    const Problem::problem_t<T> &problem_;
    Problem::hash_t<T> &hash_;
    size_t size_;
    std::vector<T> roots_;
    pair_list tips_;
    pair_list nodes_;

  public:
    policy_graph_t(const Problem::problem_t<T> &problem, Problem::hash_t<T> &hash)
      : problem_(problem), hash_(hash), size_(0) {
    }
    ~policy_graph_t() { }

    size_t size() const { return size_; }
    void add_root(const T &s) { roots_.push_back(s); }
    const pair_list& tips() const { return tips_; }
    const pair_list& nodes() const { return nodes_; }

    void recompute() {
        std::vector<std::pair<T, float> > outcomes;
        size_ = 0;
        tips_.clear();
        nodes_.clear();

        pair_list open;
        for( vector_iterator si = roots_.begin(); si != roots_.end(); ++si ) {
            Hash::data_t *dptr = hash_.data_ptr(*si);
            open.push_back(std::make_pair(*si, dptr));
            dptr->solve();
            dptr->mark();
        }

        while( !open.empty() ) {
            std::pair<T, Hash::data_t*> n = open.front();
            open.pop_front();

            if( problem_.terminal( n.first ) ) {
                n.second->set_action(Problem::noop);
                nodes_.push_front(n);
                ++size_;
            } else {
                bool unsolved = false;
                std::pair<Problem::action_t, float> p = hash_.bestQValue(n.first);
                assert(p.first != Problem::noop);
                n.second->set_action(p.first);

                problem_.next(n.first, p.first, outcomes);
                unsigned osize = outcomes.size();

                for( unsigned i = 0; i < osize; ++i ) {
                    if( !hash_.solved(outcomes[i].first) ) {
                        unsolved = true;
                        break;
                    }
                }

                ++size_;
                if( !unsolved ) {
                    for( unsigned i = 0; i < osize; ++i ) {
                        Hash::data_t *dptr = hash_.data_ptr(outcomes[i].first);
                        if( !dptr->marked() ) {
                            open.push_back(std::make_pair(outcomes[i].first, dptr));
                            dptr->solve();
                            dptr->mark();
                        }
                    }
                    nodes_.push_front(n);
                } else {
                    tips_.push_back(n);
                    n.second->set_action(Problem::noop);
                }
            }
        }

        // unmark nodes
        for( list_iterator si = nodes_.begin(); si != nodes_.end(); ++si )
            si->second->unmark();
        for( list_iterator si = tips_.begin(); si != tips_.end(); ++si )
            si->second->unmark();
    }

    void postorder_dfs(const T &s, pair_list &visited) {
        pair_list open;

        std::vector<std::pair<T, float> > outcomes;
        Hash::data_t *dptr = hash_.data_ptr(s);
        open.push_back(std::make_pair(s, dptr ));
        dptr->mark();

        while( !open.empty() ) {
            std::pair<T, Hash::data_t*> n = open.back();

	    // do a preorder backup here:
	    //std::pair<Problem::action_t, float> p = hash_.bestQValue(n.first);
            //n.second->update(p.second);
            visited.push_front(n);
            open.pop_back();
            assert(n.second->solved());

            if( n.second->action() != Problem::noop ) {
                problem_.next(n.first, n.second->action(), outcomes);
                unsigned osize = outcomes.size();
                for( unsigned i = 0; i < osize; ++i ) {
                    Hash::data_t *dptr = hash_.data_ptr(outcomes[i].first);
                    if( !dptr->marked() ) {
                        open.push_back(std::make_pair(outcomes[i].first, dptr));
                        dptr->mark();
                    }
                }
            } else {
		
                std::pair<Problem::action_t, float> p = hash_.bestQValue(n.first);
                problem_.next(n.first, p.first, outcomes);
                unsigned osize = outcomes.size();
                for( unsigned i = 0; i < osize; ++i )
                    hash_.solve(outcomes[i].first);
            }
        }

        // unmark nodes
        for( list_iterator si = visited.begin(); si != visited.end(); ++si )
            si->second->unmark();
    }

    void update(pair_list &visited) {
        for( list_iterator si = visited.begin(); si != visited.end(); ++si ) {
            std::pair<Problem::action_t, float> p = hash_.bestQValue(si->first);
            si->second->update(p.second);
            hash_.inc_updates();
        }
    }
};

template<typename T>
void generate_space(const Problem::problem_t<T> &problem,
                    const T &s,
                    Problem::hash_t<T> &hash) {
    std::list<std::pair<T, Hash::data_t*> > open;

    std::vector<std::pair<T, float> > outcomes;
    Hash::data_t *dptr = hash.data_ptr(s);
    open.push_back(std::make_pair(s, dptr));
    dptr->mark();

#ifdef DEBUG
    std::cout << "marking " << s << std::endl;
#endif

    while( !open.empty() ) {
        std::pair<T, Hash::data_t*> n = open.front();
        open.pop_front();
        if( problem.terminal(n.first) ) continue;

        for( Problem::action_t a = 0; a < problem.number_actions(n.first); ++a ) {
            if( problem.applicable(n.first, a) ) {
                problem.next(n.first, a, outcomes);
                unsigned osize = outcomes.size();
                for( unsigned i = 0; i < osize; ++i ) {
                    Hash::data_t *ptr = hash.data_ptr(outcomes[i].first);
                    if( !ptr->marked() ) {
                        open.push_back(std::make_pair(outcomes[i].first, ptr));
                        ptr->mark();
#ifdef DEBUG
                        std::cout << "marking " << outcomes[i].first << std::endl;
#endif
                    }
                }
            }
        }
    }

    hash.unmark_all();
}

template<typename T>
size_t value_iteration(const Problem::problem_t<T> &problem,
                       const T &s,
                       Problem::hash_t<T> &hash,
                       const parameters_t &parameters) {
    typedef typename Problem::hash_t<T>::iterator hash_iterator;
    generate_space(problem, s, hash);

#ifdef DEBUG
    std::cout << "state space = " << hash.size() << std::endl;
#endif

    size_t iters = 0;
    float residual = 1 + parameters.epsilon_;
    while( residual > parameters.epsilon_ ) {
        if( iters > parameters.vi.max_number_iterations_ ) break;
        residual = 0;
        for( hash_iterator hi = hash.begin(); hi != hash.end(); ++hi ) {
            float hv = hi->second->value();
            std::pair<Problem::action_t, float> p = hash.bestQValue(hi->first);
            float res = (float)fabs(p.second - hv);
            residual = Utils::max(residual, res);
            hi->second->update(p.second);
            hash.inc_updates();

#ifdef DEBUG
            if( res > parameters.epsilon_ ) {
                std::cout << "value for " << hi->first
                          << " changed from " << hv << " to "
                          << p.second << std::endl;
            }
#endif
        }
        ++iters;

#ifdef DEBUG
        std::cout << "residual=" << residual << std::endl;
#endif
    }
    return iters;
}

template<typename T>
bool check_solved(const Problem::problem_t<T> &problem,
                  const T &s,
                  Problem::hash_t<T> &hash,
                  float epsilon) {
    std::list<std::pair<T, Hash::data_t*> > open, closed;

    std::vector<std::pair<T, float> > outcomes;
    Hash::data_t *dptr = hash.data_ptr(s);
    if( !dptr->solved() ) {
        open.push_back(std::make_pair(s, dptr));
        dptr->mark();
    }

    bool rv = true;
    while( !open.empty() ) {
        std::pair<T, Hash::data_t*> n = open.back();
        closed.push_back(n);
        open.pop_back();
        if( problem.terminal(n.first) ) continue;

        std::pair<Problem::action_t, float> p = hash.bestQValue(n.first);
        if( fabs(p.second - n.second->value()) > epsilon ) {
            rv = false;
            continue;
        }

        problem.next(n.first, p.first, outcomes);
        unsigned osize = outcomes.size();

        for( unsigned i = 0; i < osize; ++ i ) {
            Hash::data_t *dptr = hash.data_ptr(outcomes[i].first);
            if( !dptr->solved() && !dptr->marked() ) {
                open.push_back(std::make_pair(outcomes[i].first, dptr));
                dptr->mark();
            }
        }
    }

    if( rv ) {
        while( !closed.empty() ) {
            closed.back().second->solve();
            closed.pop_back();
        }
    } else {
        while( !closed.empty() ) {
            std::pair<Problem::action_t, float> p = hash.bestQValue(closed.back().first);
            closed.back().second->update(p.second);
            closed.back().second->unmark();
            hash.inc_updates();
            closed.pop_back();
        }
    }
    return rv;
}

template<typename T, int V>
size_t lrtdp_trial(const Problem::problem_t<T> &problem,
                   const T &s,
                   Problem::hash_t<T> &hash,
                   const parameters_t &parameters) {
    std::list<T> states;
    std::pair<T, bool> n;

    Hash::data_t *dptr = hash.data_ptr(s);
    states.push_back(s);
    dptr->inc_count();
    T t = s;

#ifdef DEBUG
    std::cout << "trial: begin" << std::endl;
#endif

    size_t steps = 1;
    while( !problem.terminal(t) && !dptr->solved() && (dptr->count() <= parameters.rtdp.bound_) ) {

#ifdef DEBUG
        std::cout << "  " << t << " = " << dptr->value() << std::endl;
#endif

        std::pair<Problem::action_t, float> p = hash.bestQValue(t);
        dptr->update(p.second);
        hash.inc_updates();

        if( Random::real() < parameters.rtdp.epsilon_greedy_ ) {
            n = problem.usample(t, p.first);
        } else {
            if( V == 0 ) {
                n = problem.sample(t, p.first);
            } else if( V == 1 ) {
                n = problem.usample(t, p.first);
            } else if( V == 2 ) {
                n = problem.nsample(t, p.first, hash);
            }
        }
        if( !n.second ) break;

        t = n.first;
        dptr = hash.data_ptr(t);
        states.push_back(t);
        dptr->inc_count();
        ++steps;
    }

#ifdef DEBUG
    std::cout << "trial: end" << std::endl;
#endif

    while( !states.empty() ) {
        hash.clear_count(states.back());
        bool solved = check_solved(problem, states.back(), hash, parameters.epsilon_);
        states.pop_back();
        if( !solved ) break;
    }

    while( !states.empty() ) {
        hash.clear_count(states.back());
        states.pop_back();
    }
    return steps;
}

template<typename T, int V>
size_t lrtdp(const Problem::problem_t<T> &problem,
             const T &s,
             Problem::hash_t<T> &hash,
             const parameters_t &parameters) {
    size_t trials = 0, max_steps = 0;
    while( !hash.solved(s) ) {
        size_t steps = lrtdp_trial<T, V>(problem, s, hash, parameters);
        max_steps = Utils::max(max_steps, steps);
        ++trials;
    }
    return trials;
}

template<typename T>
size_t standard_lrtdp(const Problem::problem_t<T> &problem,
                      const T &s,
                      Problem::hash_t<T> &hash,
                      const parameters_t &parameters) {
    return lrtdp<T, 0>(problem, s, hash, parameters);
}

template<typename T>
size_t uniform_lrtdp(const Problem::problem_t<T> &problem,
                     const T &s,
                     Problem::hash_t<T> &hash,
                     const parameters_t &parameters) {
    return lrtdp<T, 1>(problem, s, hash, parameters);
}

template<typename T>
size_t bounded_lrtdp(const Problem::problem_t<T> &problem,
                     const T &s,
                     Problem::hash_t<T> &hash,
                     const parameters_t &parameters) {
    return lrtdp<T, 2>(problem, s, hash, parameters);
}

template<typename T>
size_t improved_lao(const Problem::problem_t<T> &problem,
                    const T &s,
                    Problem::hash_t<T> &hash,
                    const parameters_t &parameters) {

        
	Hansen_LAO<T> LLAO(problem, hash, parameters);
	LLAO.lLAOstar_driver(s);
	return LLAO.Iteration;
	
	/*
    typedef typename std::list<std::pair<T, Hash::data_t*> > pair_list;
    typedef typename pair_list::const_iterator const_list_iterator;
    pair_list visited;

    policy_graph_t<T> graph(problem, hash);
    graph.add_root(s);
    graph.recompute();

    size_t iterations = 0;

  loop:
    ++iterations;
    while( !graph.tips().empty() ) {
        visited.clear();
        graph.postorder_dfs(s, visited);
        //std::cout<<"Iteration: "<< iterations<<" Visited set size: "<<visited.size()<<std::endl;
        graph.update(visited);
        graph.recompute();
    }

    // convergence test
    float residual = 1.0 + parameters.epsilon_;
    while( residual > parameters.epsilon_ ) {
        residual = 0.0;
        for( const_list_iterator si = graph.nodes().begin(); si != graph.nodes().end(); ++si ) {
            float hv = si->second->value();
            std::pair<Problem::action_t, float> p = hash.bestQValue(si->first);
            residual = Utils::max(residual, (float)fabs(p.second - hv));
            si->second->update(p.second);
            hash.inc_updates();
        }
        graph.recompute();
        if( !graph.tips().empty() ) goto loop;
    }
    return iterations;
*/
}

template<typename T>
size_t plain_check(const Problem::problem_t<T> &problem,
                   const T &s,
                   Problem::hash_t<T> &hash,
                   const parameters_t &parameters) {
    size_t trials = 0;
    while( !hash.solved(s) ) {
        check_solved(problem, s, hash, parameters.epsilon_);
        ++trials;
    }
    return trials;
}

#if 0
template<typename T>
bool lenforce(const Problem::problem_t<T> &problem,
              Problem::hash_t<T> &hash,
              const T &s,
              Hash::data_t *dptr,
              float epsilon,
              size_t &index,
              size_t kappa,
              std::list<Hash::data_t*> &stack,
              std::list<Hash::data_t*> &visited) {
    unsigned osize = 0;
    std::pair<T, float> outcomes[MAXOUTCOMES];

    // base cases
    if( dptr->solved() || problem.terminal(s) ) {
        dptr->solve();
        return true;
    }
    assert(!dptr->marked());

    // if residual > epsilon, update and return
    std::pair<Problem::action_t, float> p = hash.bestQValue(problem, s);
    if( fabs(p.second - dptr->value()) > epsilon ) {
        dptr->update(p.second);
        hash.inc_updates();
        return false;
    }

    // Tarjan's
    visited.push_front(dptr);
    stack.push_front(dptr);
    size_t idx = index++;
    dptr->set_scc_low(idx);
    dptr->set_scc_idx(idx);
    dptr->mark();

    // compute normalization constant
    size_t normalization = std::numeric_limits<unsigned>::max();
    problem.next(s, p.first, outcomes, osize);
    for( unsigned i = 0; i < osize; ++i ) {
        size_t kv = Utils::kappa_value(outcomes[i].second);
        normalization = Utils::min(normalization, kv);
    }

    // recursive call
    bool flag = true;
    for( unsigned i = 0; i < osize; ++i ) {
        size_t normalized_kappa = Utils::kappa_value(outcomes[i].second) - normalization;
        if( normalized_kappa <= kappa ) {
            size_t nkappa = kappa - normalized_kappa;
            Hash::data_t *ptr = hash.data_ptr(outcomes[i].first);
            if( ptr->scc_idx() == std::numeric_limits<unsigned>::max() ) {
                bool rv = lenforce(problem, hash, outcomes[i].first, ptr, epsilon, index, nkappa, stack, visited);
                flag = flag && rv;
                dptr->set_scc_low(Utils::min(dptr->scc_low(), ptr->scc_low()));
            } else if( ptr->marked() ) {
                dptr->set_scc_low(Utils::min(dptr->scc_low(), ptr->scc_idx()));
            }
        }
    }

    // update
    if( !flag ) {
        std::pair<Problem::action_t, float> p = hash.bestQValue(s);
        dptr->update( p.second );
        hash.inc_updates();
        while( !stack.empty() && (stack.front()->scc_idx() >= idx) ) {
            stack.front()->set_scc_low(std::numeric_limits<unsigned>::max());
            stack.front()->set_scc_idx(std::numeric_limits<unsigned>::max());
            stack.pop_front();
        }
    } else if( dptr->scc_low() == dptr->scc_idx() ) {
        while( !stack.empty() && (stack.front()->scc_idx() >= idx) ) {
            stack.front()->solve();
            stack.front()->unmark();
            stack.front()->set_scc_low(std::numeric_limits<unsigned>::max());
            stack.front()->set_scc_idx(std::numeric_limits<unsigned>::max());
            stack.pop_front();
        }
    }
    return flag;
}
#endif

#if 0
template<typename T>
bool enforce(const Problem::problem_t<T> &problem,
             Problem::hash_t<T> &hash,
             const T &s,
             Hash::data_t *dptr,
             float epsilon,
             size_t kappa,
             std::list<Hash::data_t*> &visited) {
    unsigned osize = 0;
    std::pair<T, float> outcomes[MAXOUTCOMES];

    // base cases
    if( dptr->marked() || dptr->solved() || problem.terminal(s) ) {
        return true;
    }

    // mark node
    visited.push_front(dptr);
    dptr->mark();

    // if residual > epsilon, update and return
    std::pair<Problem::action_t, float> p = hash.bestQValue(s);
    if( fabs(p.second - dptr->value()) > epsilon ) {
        dptr->update(p.second);
        hash.inc_updates();
        return false;
    }

    // compute normalization constant
    size_t normalization = std::numeric_limits<unsigned>::max();
    problem.next(s, p.first, outcomes, osize);
    for( unsigned i = 0; i < osize; ++i ) {
        size_t kv = Utils::kappa_value(outcomes[i].second);
        normalization = Utils::min(normalization, kv);
    }

    // recursive call
    bool flag = true;
    for( unsigned i = 0; i < osize; ++i ) {
        size_t normalized_kappa = Utils::kappa_value(outcomes[i].second) - normalization;
        if( normalized_kappa <= kappa ) {
            size_t nkappa = kappa - normalized_kappa;
            bool rv = enforce(problem, hash, outcomes[i].first, hash.data_ptr( outcomes[i].first ), epsilon, nkappa, visited);
            flag = flag && rv;
        }
    }

    if( !flag ) {
        std::pair<Problem::action_t, float> p = hash.bestQValue(s);
        dptr->update(p.second);
        hash.inc_updates();
    }
    return flag;
}
#endif

#if 0
template<typename T>
size_t hdp_i(const Problem::problem_t<T> &problem,
             const T &s,
             Problem::hash_t<T> &hash,
             float epsilon,
             size_t kappa,
             float) {
    std::list<Hash::data_t*> visited;
    bool status = false;
    size_t trials = 0;
    while( !status ) {
        status = enforce(problem, hash, s, hash.data_ptr(s), epsilon, kappa, visited);
        while( !visited.empty() ) {
            visited.front()->unmark();
            visited.pop_front();
        }
        ++trials;
    }
    return trials;
}
#endif




/**
 *Labeling DFS from the root of an SCC
 *
 */

template<typename T>
void dfs_label(const Problem::problem_t<T> &problem,
               const T &s,
               Problem::hash_t<T> &hash,
               Hash::data_t* dptr
            ) {
    std::vector<std::pair<T, float> > outcomes;
    dptr->mark();
    //base case
    if( dptr->solved() || problem.terminal(s) ) {
        return;
    } 
    
    // expansion
    problem.next(s, dptr->action(), outcomes);
    unsigned osize = outcomes.size();

    for( unsigned i = 0; i < osize; ++i ) {
        Hash::data_t *ptr = hash.data_ptr(outcomes[i].first);
        if( !ptr->marked() ) {//state not seen this trial
            dfs_label(problem, outcomes[i].first, hash, ptr);
        }
    }
    dptr->solve();

}



/**
 *Labeling procedure unsed in the new LFVI
 *
 */

template<typename T>
void LabelSolve(const Problem::problem_t<T> &problem,
               const T &s,
               Problem::hash_t<T> &hash,
               Hash::data_t* dptr
            ) {
    std::vector<std::pair<T, float> > outcomes;
    dptr->mark();
    //base case
    if( dptr->solved() || problem.terminal(s)){
        return;
    } 
    
    // expansion
    problem.next(s, dptr->action(), outcomes);
    unsigned osize = outcomes.size();

    for( unsigned i = 0; i < osize; ++i ) {
        Hash::data_t *ptr = hash.data_ptr(outcomes[i].first);
        if( !ptr->marked() && (ptr->scc_low() == dptr->scc_low())) {//state not seen this trial
            LabelSolve(problem, outcomes[i].first, hash, ptr);
        }
    }
	std::cout<<"state marked as solved..."<<std::endl;
    dptr->solve();

}




/**
 *Labeling procedure to label states unsolvable in the LFVI
 *
 */

template<typename T>
void LabelUnsolvable(const Problem::problem_t<T> &problem,
               const T &s,
               Problem::hash_t<T> &hash,
               Hash::data_t* dptr
            ) {
    std::vector<std::pair<T, float> > outcomes;
    dptr->mark();
    //base case
    if(dptr->unsolvabled()) {
        return;
    } 
    
    // expansion
    problem.next(s, dptr->action(), outcomes);
    unsigned osize = outcomes.size();

    for( unsigned i = 0; i < osize; ++i ) {
        Hash::data_t *ptr = hash.data_ptr(outcomes[i].first);
        if( !ptr->marked() && (ptr->scc_low() == dptr->scc_low())) {//state not seen this trial
            LabelUnsolvable(problem, outcomes[i].first, hash, ptr);
        }
    }
	std::cout<<"state marked as unsolvable..."<<std::endl;
    dptr->unsolvable();
	dptr->update(999999.99);

}	




/**
 *Labeling procedure for labeling deadend in the LFVI
 *
 */

template<typename T>
void LabelDeadend(const Problem::problem_t<T> &problem,
               const T &s,
               Problem::hash_t<T> &hash,
               Hash::data_t* dptr
            ) {
    std::vector<std::pair<T, float> > outcomes;
    dptr->mark();
    //base case
    if(dptr->deadended()) {
        return;
    } 
    
    // expansion
    problem.next(s, dptr->action(), outcomes);
    unsigned osize = outcomes.size();

    for( unsigned i = 0; i < osize; ++i ) {
        Hash::data_t *ptr = hash.data_ptr(outcomes[i].first);
        if( !ptr->marked() && (ptr->scc_low() == dptr->scc_low())) {//state not seen this trial
            LabelDeadend(problem, outcomes[i].first, hash, ptr);
        }
    }
	std::cout<<"state marked as deadend..."<<std::endl;
    dptr->deadend();
	dptr->update(999999.99);
	

}	

	

/*
 * Orignial HDP
 */

template<typename T>
bool hdp(const Problem::problem_t<T> &problem,
         const T &s,
         Problem::hash_t<T> &hash,
         Hash::data_t* dptr,
         size_t &index,
         std::list<Hash::data_t*> &stack,
         std::list<Hash::data_t*> &visited,
         const parameters_t &parameters) {

    std::vector<std::pair<T, float> > outcomes;

    // base cases
    if( dptr->solved() || problem.terminal(s) ) {
        dptr->solve();
        return true;
    } else if( dptr->marked() ) {
        return false;
    }

    // if residual > epsilon, update and return
    std::pair<Problem::action_t, float> p = hash.bestQValue(s);
    if( fabs(p.second - dptr->value()) > parameters.epsilon_ ) {
        dptr->update(p.second);
        hash.inc_updates();
        return false;
    }

    // Tarjan's
    visited.push_front(dptr);
    stack.push_front(dptr);
    size_t idx = index++;
    dptr->set_scc_low(idx);
    dptr->set_scc_idx(idx);
    dptr->mark();

    // expansion
    problem.next(s, p.first, outcomes);
    unsigned osize = outcomes.size();

    bool flag = true;
    for( unsigned i = 0; i < osize; ++i ) {
        Hash::data_t *ptr = hash.data_ptr(outcomes[i].first);
        if( ptr->scc_idx() == std::numeric_limits<unsigned>::max() ) {
            bool rv = hdp(problem, outcomes[i].first, hash, ptr, index, stack, visited, parameters);
            flag = flag && rv;
            dptr->set_scc_low(Utils::min(dptr->scc_low(), hash.scc_low(outcomes[i].first)));
        } else if( ptr->marked() ) {
            dptr->set_scc_low(Utils::min(dptr->scc_low(), ptr->scc_idx()));
        }
    }

    // update
    if( !flag ) { //commenting these 3 lines below will give HDP2 = HDP -post order backup if not labeling
        std::pair<Problem::action_t, float> p = hash.bestQValue(s);
        dptr->update(p.second);
        hash.inc_updates();
        while( !stack.empty() && (stack.front()->scc_idx() >= idx) ) {
            stack.front()->set_scc_low(std::numeric_limits<unsigned>::max());
            stack.front()->set_scc_idx(std::numeric_limits<unsigned>::max());
            stack.pop_front();
        }
    } else if( dptr->scc_low() == dptr->scc_idx() ) {
        while( !stack.empty() && (stack.front()->scc_idx() >= idx) ) {
            stack.front()->solve();
            stack.front()->set_scc_low(std::numeric_limits<unsigned>::max());
            stack.front()->set_scc_idx(std::numeric_limits<unsigned>::max());
            stack.pop_front();
        }
    }
    return flag;
}

template<typename T>
size_t hdp_driver(const Problem::problem_t<T> &problem,
                  const T &s,
                  Problem::hash_t<T> &hash,
                  const parameters_t &parameters) {
    std::list<Hash::data_t*> stack, visited;
    Hash::data_t *dptr = hash.data_ptr(s);
    size_t trials = 0;
    std::cout<<"Starting here: "<<std::endl;
    while( !dptr->solved() ) {
        size_t index = 0;
        hdp(problem, s, hash, dptr, index, stack, visited, parameters);
        assert(stack.empty());
        while( !visited.empty() ) {
            visited.front()->unmark();
            assert(visited.front()->scc_low() == std::numeric_limits<unsigned>::max());
            assert(visited.front()->scc_idx() == std::numeric_limits<unsigned>::max());
            visited.pop_front();
        }
        ++trials;
    }
    return trials;
}


    //************************************************LFVI using stack***********************************************
    
    /**
     * Labeled Focused Value Iteration with stack
     *
     **/
    //Driver
    template<typename T>
    size_t lfvi_driver(const Problem::problem_t<T> &problem,
                       const T &s,
                       Problem::hash_t<T> &hash,
                       const parameters_t &parameters) {
        std::deque<Hash::data_t*> stack, visited;
        Hash::data_t *dptr = hash.data_ptr(s);
        size_t trials = 1;
        double threshold = parameters.epsilon_;
        double residual = 0.0;
        size_t numSccs = 0;
        bool pichange;
        while( !dptr->solved() ) {
            size_t index = 0;
            pichange = false;
            lfvi_dfs(problem, s, hash, dptr, stack, visited, index, trials, threshold, residual, numSccs, pichange, parameters);
            assert(stack.empty());
            
            visited.erase(visited.begin(), visited.end());
            ++trials;
        }
        return trials;
    }
    
    
    //DFS
    template<typename T>
    double lfvi_dfs(const Problem::problem_t<T> &problem,
                    const T &s,
                    Problem::hash_t<T> &hash,
                    Hash::data_t* dptr, std::deque<Hash::data_t*> &stack,
                    std::deque<Hash::data_t*> &visited,
                    size_t &index,
                    size_t trials,
                    double threshold,
                    double &residual,
                    size_t &numSccs, bool &pichange,
                    const parameters_t &parameters) {
        
        std::vector<std::pair<T, float> > outcomes;
        
        
        // base cases
        if( dptr->solved() || problem.terminal(s) ) {
            dptr->solve();
            return 0;
        }
        
        //Backup to compute greedy action
        float diff = hash.Backup(s);
        
        if (diff > residual)
        residual = diff;
        
        
        // Tarjan's
        visited.push_front(dptr);
        stack.push_front(dptr);
        size_t idx = index++;
        dptr->set_scc_low(idx);
        dptr->set_scc_idx(idx);
        dptr->SetUpdate(trials);
        
        
        
        
        
        
        // expansion
        problem.next(s, dptr->action(), outcomes);
        unsigned osize = outcomes.size();
        for( unsigned i = 0; i < osize; ++i )
        {
            Hash::data_t *ptr = hash.data_ptr(outcomes[i].first);
            if( ptr->Update() < trials ) {//state not seen this trial
                lfvi_dfs(problem, outcomes[i].first, hash, ptr, stack, visited, index, trials, threshold, residual, numSccs, pichange, parameters);
                dptr->set_scc_low(Utils::min(dptr->scc_low(), ptr->scc_low()));
            } else if( ptr->Update() == trials ) {
                dptr->set_scc_low(Utils::min(dptr->scc_low(), ptr->scc_idx()));
            }
        }
        
        float postres = 0.0;
        //post order backup
        diff = hash.SBackup(s, pichange);
        if(diff > postres)
        postres = diff;
        
        if( pichange ) { //policy changed, no labeling
            //hash.inc_updates();
            while( !stack.empty() && (stack.front()->scc_idx() >= idx) ) {
                stack.pop_front();
            }
        } else if((postres < threshold) && (dptr->scc_low() == dptr->scc_idx()) ) { //root of an scc.
            while( !stack.empty() && (stack.front()->scc_idx() >= idx) ) {
                stack.front()->solve();
                stack.pop_front();
            }
        }else{
            stack.erase(stack.begin(), stack.end());
        }
        
        
        return residual;
    }
    
    
    
    




/**
 * *****************************Labeled FVI without a stack***********************************************
 *
 */	
template<typename T>
bool nslfvi(const Problem::problem_t<T> &problem,
				const T &s,
				Problem::hash_t<T> &hash,
				Hash::data_t* dptr,
				size_t &index,
            size_t &trials, double &threshold,
            double &residual, unsigned &numSccs, bool &pichange,
				const parameters_t &parameters) {
	
    std::vector<std::pair<T, float> > outcomes;
	
    // base cases
    if( dptr->solved() || problem.terminal(s) ) {
        dptr->solve();
        return true;
    }
    
    //Backup to compute greedy action
    float diff = hash.Backup(s);
    
    if (diff > residual)
        residual = diff;
    // Tarjan's
    
    
    dptr->set_scc_low(index);
    dptr->set_scc_idx(index);
    dptr->SetUpdate(trials);
	++index;
    
    // expansion
	problem.next(s, dptr->action(), outcomes);
	unsigned osize = outcomes.size();
	
	bool flag = true;
	for( unsigned i = 0; i < osize; ++i )
	{
        Hash::data_t *ptr = hash.data_ptr(outcomes[i].first);
        if( ptr->Update() < trials ) {//state not seen this trial
			bool rv = nslfvi(problem, outcomes[i].first, hash, ptr, index, trials, threshold, residual, numSccs, pichange, parameters);
			flag = flag && rv;
			dptr->set_scc_low(Utils::min(dptr->scc_low(), ptr->scc_low()));
        } else if( ptr->Update() == trials ) {
            dptr->set_scc_low(Utils::min(dptr->scc_low(), ptr->scc_idx()));
        }
    }
    
    float postres = 0.0;
    //post order backup
    diff = hash.SBackup(s, pichange);
    if(diff > postres)
    postres = diff;
    
    //label if no policy change, residual small enough and found a root of an scc.
    
    if(!pichange && (postres < threshold) && (dptr->scc_low() == dptr->scc_idx()))
    {
       dfs_label(problem, s, hash, dptr);
    }

    return flag;
}



template<typename T>
size_t nslfvi_driver(const Problem::problem_t<T> &problem,
                  const T &s,
                  Problem::hash_t<T> &hash,
                  const parameters_t &parameters) {
    double threshold = parameters.epsilon_;
    double residual = 0.0;
    Hash::data_t *dptr = hash.data_ptr(s);
    unsigned numSccs = 0;
    bool pichange;
    size_t trials = 1;//important to start with 1  for marking purpose
    
    while( !dptr->solved())
	{
        residual = 0.0;
        numSccs = 0;
        pichange = false;
        size_t index = 0;
        nslfvi(problem, s, hash, dptr, index, trials, threshold, residual, numSccs, pichange, parameters);
        ++trials;
    }
	return trials;
}



/**
 * Focused Value Iteration DFS
 *
 **/

template<typename T>
bool fvi_dfs(const Problem::problem_t<T> &problem,
			 const T &s,
			 Problem::hash_t<T> &hash,
			 Hash::data_t* dptr,
			 size_t &trials,
			 double &residual,
			 const parameters_t &parameters) {
	
    std::vector<std::pair<T, float> > outcomes;

    // base cases
    if(problem.terminal(s) ) {
		return true;
    } 
	
    //Backup
	float diff = hash.Backup(s);
    if(diff > residual)
		residual = diff;
	
    // expansion
	problem.next(s, dptr->action(), outcomes);
	unsigned osize = outcomes.size();
	
	for( unsigned i = 0; i < osize; ++i )
	{
        Hash::data_t *ptr = hash.data_ptr(outcomes[i].first);
        if( ptr->Update() < trials ) {//state not seen this trial
			ptr->SetUpdate(trials);
			fvi_dfs(problem, outcomes[i].first, hash, ptr, trials, residual, parameters);
		} 
    }
	hash.Backup(s);
	return true;
}


/**
 *Focused VI algorithm's driver
 */

template<typename T>
size_t fvi_driver(const Problem::problem_t<T> &problem,
                  const T &s,
                  Problem::hash_t<T> &hash,
                  const parameters_t &parameters)
{
	Hash::data_t *dptr = hash.data_ptr(s);
    size_t trials;
	double  error;
	double  upperbound;
	double residual;
	std::cout << std::setprecision(8)<<std::endl;
	
    for( trials = 1; ; trials++)
	{
		residual = 0.0;
		dptr->SetUpdate(trials);
		fvi_dfs(problem, s, hash, dptr, trials, residual, parameters);
	    
		if(residual < 1.0)
		{
			upperbound = (dptr->value() - residual)/(1.0 - residual);
			error = (upperbound - dptr->value());
			
		}else
		{
			error = std::numeric_limits<unsigned>::max();	
		}
		if(error <= parameters.epsilon_)
			break;
		
		
    }
    return trials;
}

//******************************************END FVI****************************************************************

/*
 *Ldfs algorithm
 */

template<typename T, int V>
bool ldfs(const Problem::problem_t<T> &problem,
          const T &s,
          Problem::hash_t<T> &hash,
          Hash::data_t *dptr,
          size_t &index,
          std::list<Hash::data_t*> &stack,
          std::list<Hash::data_t*> &visited,
          const parameters_t &parameters) {

    std::vector<std::pair<T, float> > outcomes;

    // base cases
    if( dptr->solved() || problem.terminal(s) ) {
        dptr->solve();
        return true;
    } else if( dptr->marked() ) {
        return false;
    }

    // Tarjan's
    visited.push_front(dptr);
    stack.push_front(dptr);
    size_t idx = index++;
    dptr->set_scc_low(idx);
    dptr->set_scc_idx(idx);
    //dptr->mark();

    if( V == 1 ) {
        std::pair<Problem::action_t, float> p = hash.bestQValue(s);
        dptr->update(p.second);
        hash.inc_updates();
    }

    // expansion
    bool flag = false;
    float bqv = std::numeric_limits<float>::max();
    for( Problem::action_t a = 0; a < problem.number_actions(s); ++a ) {
        if( problem.applicable(s, a) ) {
            problem.next(s, a, outcomes);
            unsigned osize = outcomes.size();

            float qv = 0.0;
            for( unsigned i = 0; i < osize; ++i )
                qv += outcomes[i].second * hash.value(outcomes[i].first);
            qv = problem.cost(s, a) + problem.discount() * qv;
            bqv = Utils::min(bqv, qv);
            if( fabs(qv - dptr->value()) > parameters.epsilon_ ) continue;
            dptr->mark();
            flag = true;
            for( unsigned i = 0; i < osize; ++i ) {
                Hash::data_t *ptr = hash.data_ptr(outcomes[i].first);
                if( ptr->scc_idx() == std::numeric_limits<unsigned>::max() ) {
                    bool rv = ldfs<T, V>(problem, outcomes[i].first, hash, ptr, index, stack, visited, parameters);
                    flag = flag && rv;
                    dptr->set_scc_low(Utils::min(dptr->scc_low(), ptr->scc_low()));
                } else if( ptr->marked() ) {
                    dptr->set_scc_low(Utils::min(dptr->scc_low(), ptr->scc_idx()));
                }
            }
            if( (V == 1) && flag && (hash.QValue(s, a) - dptr->value() > parameters.epsilon_) ) flag = false;
            if( flag ) break;
            while( stack.front()->scc_idx() > idx ) {
                stack.front()->set_scc_low(std::numeric_limits<unsigned>::max());
                stack.front()->set_scc_idx(std::numeric_limits<unsigned>::max());
                stack.pop_front();
            }
        }
    }

    // update
    if( !flag ) {
        if( !dptr->marked() ) {
            dptr->update(bqv);
        } else {
            std::pair<Problem::action_t, float> p = hash.bestQValue(s);
            dptr->update(p.second);
        }
        hash.inc_updates();
        dptr->set_scc_low(std::numeric_limits<unsigned>::max());
        dptr->set_scc_idx(std::numeric_limits<unsigned>::max());
        stack.pop_front();
    } else if( dptr->scc_low() == dptr->scc_idx() ) {
        while( !stack.empty() && (stack.front()->scc_idx() >= idx) ) {
            stack.front()->solve();
            stack.front()->set_scc_low(std::numeric_limits<unsigned>::max());
            stack.front()->set_scc_idx(std::numeric_limits<unsigned>::max());
            stack.pop_front();
        }
    }
    return flag;
}

template<typename T>
size_t ldfs_plus_driver(const Problem::problem_t<T> &problem,
                        const T &s,
                        Problem::hash_t<T> &hash,
                        const parameters_t &parameters) {
    std::list<Hash::data_t*> stack, visited;
    size_t trials = 0;
    Hash::data_t *dptr = hash.data_ptr(s);
    while( !dptr->solved() ) {
        size_t index = 0;
        ldfs<T, 1>(problem, s, hash, dptr, index, stack, visited, parameters);
        assert(stack.empty());
        while( !visited.empty() ) {
            visited.front()->unmark();
            visited.pop_front();
        }
        ++trials;
    }
    return trials;
}

template<typename T>
size_t ldfs_driver(const Problem::problem_t<T> &problem,
                   const T &s,
                   Problem::hash_t<T> &hash,
                   const parameters_t &parameters) {
    std::list<Hash::data_t*> stack, visited;
    size_t trials = 0;
    Hash::data_t *dptr = hash.data_ptr(s);
    while( !dptr->solved() ) {
        size_t index = 0;
        ldfs<T, 0>(problem, s, hash, dptr, index, stack, visited, parameters);
        assert(stack.empty());
        while( !visited.empty() ) {
            visited.front()->unmark();
            visited.pop_front();
        }
        ++trials;
    }
    return trials;
}










}; // namespace Algorithm

#undef DEBUG

#endif

