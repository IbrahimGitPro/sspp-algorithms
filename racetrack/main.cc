#include <iostream>
#include <fstream>
#include <vector>

#include "race.h"
#include <dispatcher.h>

using namespace std;

void usage(ostream &os) {
    os << "usage: race [-a <n>] [-e <f>] [-f] [-h <n>] [-p <f>] [-s <n>] <file>"
       << endl << endl
       << "  -a <n>    Algorithm bitmask: 1=vi, 2=lrtdp, 4=lao*, 8=hdp, 16=ldfs+, 32=fvi."
       << endl
       << "  -e <f>    Epsilon. Default: 0."
       << endl
       << "  -h <n>    Heuristics: 0=zero, 1=minmin. Default: 0."
       << endl
       << "  -p <f>    Parameter p in [0,1]. Default: 0.9"
       << endl
       << "  -s <n>    Random seed. Default: 0."
       << endl
       << "  <file>    Racetrack file."
       << endl << endl;
}

int main(int argc, const char **argv) {
    FILE *is = 0;
    float p = 1.0;
    unsigned bitmap = 0;
    int h = 0;
    bool formatted = true;
    string base_name;
    string policy_type;
    Online::Evaluation::parameters_t eval_pars;

    cout << fixed;
    Algorithm::parameters_t alg_pars;

    cout << "Arguments:";
    for( int i = 0; i < argc; ++i ) {
        cout << " " << argv[i];
    }
    cout << endl;

    // parse arguments
    ++argv;
    --argc;
    while( argc > 1 ) {
        if( **argv != '-' ) break;
        switch( (*argv)[1] ) {
            case 'a':
                bitmap = strtoul(argv[1], 0, 0);
                argv += 2;
                argc -= 2;
                break;
            case 'e':
                alg_pars.epsilon_ = strtod(argv[1], 0);
                argv += 2;
                argc -= 2;
                break;
            case 'f':
                formatted = true;
                ++argv;
                --argc;
                break;
            case 'h':
                h = strtol(argv[1], 0, 0);
                argv += 2;
                argc -= 2;
                break;
            case 'l':
                eval_pars.labeling_ = true;
                ++argv;
                --argc;
                break;
            case 'p':
                p = strtod(argv[1], 0);
                argv += 2;
                argc -= 2;
                break;
            case 's':
                alg_pars.seed_ = strtoul(argv[1], 0, 0);
                argv += 2;
                argc -= 2;
                break;
            case 't':
                eval_pars.evaluation_trials_ = strtoul(argv[1], 0, 0);
                argv += 2;
                argc -= 2;
                break;
            default:
                usage(cout);
                exit(-1);
        }
    }

    if( argc == 1 ) {
        is = fopen(argv[0], "r");
    } else {
        usage(cout);
        exit(-1);
    }

	std::cout<<"build problem instances..."<<std::endl;
    cout << "seed=" << alg_pars.seed_ << endl;
    Random::set_seed(alg_pars.seed_);
    grid_t grid;
    grid.parse(cout, is);
    problem_t problem(grid, p);
    fclose(is);

	std::cout<<"create heuristic...."<<std::endl; 
    vector<pair<const Heuristic::heuristic_t<state_t>*, string> > heuristics;
    heuristics.push_back(make_pair(new Heuristic::zero_heuristic_t<state_t>, "zero"));
    //heuristics.push_back(make_pair(new Heuristic::min_min_heuristic_t<state_t>(problem, divisor), "min-min"));

    Heuristic::heuristic_t<state_t> *heuristic = 0;
    if( h == 0 ) {
        heuristic = new Heuristic::zero_heuristic_t<state_t>;
    }
	else if( h == 1 ) {
        heuristic = new Heuristic::min_min_heuristic_t<state_t>(problem, 1);
    }

	std::cout<<"solve problem...."<<std::endl;
    vector<Dispatcher::result_t<state_t> > results;
    Dispatcher::solve(problem, heuristic, problem.init(), bitmap, alg_pars, results);

    // print results
    if( !results.empty() ) {
        if( formatted ) Dispatcher::print_result<state_t>(cout, 0);
        for( unsigned i = 0; i < results.size(); ++i ) {
            Dispatcher::print_result(cout, &results[i]);
        }
    }
    // free resources
    for( unsigned i = 0; i < results.size(); ++i ) {
        delete results[i].hash_;
    }
 
    delete heuristic;

    exit(0);
}

