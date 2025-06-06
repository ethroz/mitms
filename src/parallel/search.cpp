/*--------------------------------------------------------------------
MITMS - Meet-in-the-middle quantum circuit synthesis
Copyright (C) 2013  Matthew Amy and The University of Waterloo,
Institute for Quantum Computing, Quantum Circuits Group

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Author: Matthew Amy
---------------------------------------------------------------------*/

#include "search.h"
#include <string>
#include <iomanip>
#include <pthread.h>
#include "vptree.h"
#include <vector>
#include <atomic>
#include <unistd.h>
#include <thread>
#include <chrono>

/* There is an insane amount of code repetition here. Don't judge me */

/* ------------- For use with approximate searching, for some reason I don't remember */
class map_value_iter : public std::iterator<bidirectional_iterator_tag, Circuit> {
	private:
		map_iter it;
	public:
		map_value_iter(map_iter x) : it(x) { }
		map_value_iter(const map_value_iter & mit) : it(mit.it) { }
		map_value_iter & operator++()    { ++it; return *this; }
		map_value_iter   operator++(int) { map_value_iter tmp(*this); operator++(); return tmp; }
		bool operator==(const map_value_iter & rhs) { return it == rhs.it; }
		bool operator!=(const map_value_iter & rhs) { return it != rhs.it; }
		Circuit & operator*() { return it->second; }
};

/*-------- threading stuff */
pthread_cond_t data_ready;
pthread_cond_t thrd_ready;
pthread_mutex_t data_lock;
pthread_mutex_t prnt_lock;
pthread_mutex_t result_mutex;

map_t            * data_map = NULL;
ord_circuit_list * data_res = NULL;
Circuit            data_circ;
Rmatrix            data_mat;
int                data_k = 0;
bool               data_left = false;
bool               data_avail = false;
int                data_num = 0;

void             * data_nntree = NULL;
void             * data_mat_approx = NULL;
int                data_l = 0;
list<pair<double, Circuit> > * data_res_approx;

/* -------------------------------- */

bool check_it(const Circuit & x, const Circuit & y, const Rmatrix & target) {
  Rmatrix tmp1(dim, dim), tmp2(dim, dim);
  x.to_Rmatrix(tmp1);
  y.to_Rmatrix(tmp2);

  tmp1 *= tmp2;
  if (dim == dim_proj) {
    return target.phase_eq(tmp1);
  } else {
    tmp1.submatrix(0, 0, dim_proj, dim_proj, tmp2);
    return target.phase_eq(tmp2);
  }

}

// Structure for passing search parameters to a worker thread.
struct WorkerArg {
  const map_t* mp;                    // pointer to the map to search through
  size_t start_index;                 // starting index (inclusive) in map
  size_t end_index;                   // ending index (exclusive) in map
  const Rmatrix* U;                   // pointer to U (the target matrix)
  bool left_multiply;                 // parameter used in canonicalization and check_it calls
  map_t* current_map;                 // current search map (e.g. mp[i] in original code)
  ord_circuit_list* res_list;         // shared result list
  std::atomic<size_t>* entries_processed;  // shared progress counter for map entries
};

// Worker thread function that processes its section of the map
void* worker_section(void* arg) {
  WorkerArg* warg = (WorkerArg*)arg;
  int pe = config::mod_perms ? num_perms : 1;
  int in = config::mod_invs ? 1 : 2;
  
  // Get iterator to start position
  auto it = warg->mp->begin();
  std::advance(it, warg->start_index);
  
  // Iterate over assigned map entries
  for (size_t idx = warg->start_index; idx < warg->end_index && it != warg->mp->end(); ++idx, ++it) {
    Circuit candidate = it->second;
    Rmatrix candidateMatrix(dim, dim);
    candidate.to_Rmatrix(candidateMatrix);
    
    // Try each permutation/inverse combination
    for (int k = 0; k < 2 * pe; k += in) {
      Rmatrix V(dim, dim), W(dim_proj, dim);
      // Depending on parity of k, perform the appropriate permutation
      if (k % 2 == 0) {
        candidateMatrix.permute_adj(V, k/2);
      } else {
        candidateMatrix.permute(V, k/2);
      }
      
      // Use the (possibly reduced) matrix for canonicalization
      Canon* canon_form;
      if (dim != dim_proj) {
        V.submatrix(0, 0, dim_proj, dim, W);
        canon_form = warg->left_multiply
                     ? canonicalize(W * (*(warg->U)), false, false, false)
                     : canonicalize((*(warg->U)) * W, false, false, false);
      } else {
        canon_form = warg->left_multiply
                     ? canonicalize(V * (*(warg->U)))
                     : canonicalize((*(warg->U)) * V);
      }
      
      auto trip = &(canon_form->front());
      Circuit ans = find_unitary(trip->key, trip->mat, *(warg->current_map)).second;
      if (!ans.empty()) {
        // Generate modified candidate circuits based on k
        Circuit tmp_circ = (k == 0)
                           ? candidate
                           : candidate.transform(k/2, k % 2 == 1);
        Circuit tmp_circ2 = ans.transform(-(trip->permutation), trip->adjoint);
        if (warg->left_multiply
             ? check_it(tmp_circ, tmp_circ2, *(warg->U))
             : check_it(tmp_circ2, tmp_circ, *(warg->U))) {
          Circuit tmp_circ3 = warg->left_multiply
                              ? tmp_circ.append(tmp_circ2)
                              : tmp_circ2.append(tmp_circ);
          int cst = tmp_circ3.cost();
          // Lock and insert into shared results
          pthread_mutex_lock(&result_mutex);
          warg->res_list->insert(ord_circuit_pair(cst, tmp_circ3));
          pthread_mutex_unlock(&result_mutex);
        }
        if (k != 0)
          delete_circuit(tmp_circ);
        delete_circuit(tmp_circ2);
      }
      canon_form->clear();
      delete canon_form;
    }
    
    // Update progress after processing each map entry
    (*warg->entries_processed)++;
  }
  
  pthread_exit(NULL);
}

// Reworked exact_search function.
void exact_search(Rmatrix & U) {
  int pe = config::mod_perms ? num_perms : 1;
  int in = config::mod_invs  ? 1 : 2;
  map_t* circ_table = new map_t[config::max_seq];
  map_t* left_table = (config::ancilla == 0) ? NULL : new map_t[config::max_seq];
  map_t* mp = (config::ancilla == 0) ? circ_table : left_table;
  circuit_list* base_list;
  ord_circuit_list* res_list = data_res = new ord_circuit_list;
  struct timespec start, end;
  data_mat = Rmatrix(dim, dim);

  // Generate base circuits.
  base_list = generate_base_circuits();

  // Initialize the circuit table for sequence 0.
  load_sequences(0, base_list, circ_table, NULL);

  // For each sequence length i.
  for (int i = 1; i < config::max_seq; i++) {
    if (config::ancilla == 0) {
      load_sequences(i, base_list, circ_table, NULL);
    } else {
      load_sequences(i, base_list, left_table, circ_table);
    }
    // Set the current search map (data_map in the old code).
    map_t* current_map = mp + i;

    for (int j = max(i-1, 1); j <= i; j++) {
      cout << "Looking for circuits with depth " << 2*i - (i - j) << "...\n";
      cout << "|";

      int total_map_entries = mp[j].size();
      std::atomic<size_t> entries_processed(0);
      int num_threads = config::num_threads;
      if (total_map_entries < num_threads)
        num_threads = total_map_entries;
      if (num_threads <= 0)
        num_threads = 1;

      std::vector<pthread_t> threads(num_threads);
      std::vector<WorkerArg> args(num_threads);
      pthread_mutex_init(&result_mutex, NULL);

      // Partition the map entries evenly among the worker threads
      int entries_per_thread = (total_map_entries + num_threads - 1) / num_threads;
      clock_gettime(CLOCK_MONOTONIC, &start);
      for (int t = 0; t < num_threads; t++) {
        int start_index = t * entries_per_thread;
        int end_index = std::min(start_index + entries_per_thread, total_map_entries);
        args[t].mp = &mp[j];
        args[t].start_index = start_index;
        args[t].end_index = end_index;
        args[t].U = &U;
        args[t].left_multiply = false;
        args[t].current_map = current_map;
        args[t].res_list = res_list;
        args[t].entries_processed = &entries_processed;
        pthread_create(&threads[t], NULL, worker_section, &args[t]);
      }

      // Log progress while the workers are busy
      while (entries_processed < total_map_entries) {
        int progress = (entries_processed * 37) / total_map_entries;
        cout << "\r|";
        for (int p = 0; p < progress; p++) {
          cout << "=";
        }
        cout << flush;
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
      }
      cout << "|\n";
      clock_gettime(CLOCK_MONOTONIC, &end);
      double elapsed = (end.tv_sec + (double)end.tv_nsec / 1e9) -
                       (start.tv_sec + (double)start.tv_nsec / 1e9);
      cout << fixed << setprecision(3);
      cout << "Time: " << elapsed << " s\n";

      // Join all worker threads.
      for (int t = 0; t < num_threads; t++) {
        pthread_join(threads[t], NULL);
      }

      // Print any results found.
      if (!res_list->empty()) {
        for (auto ti = res_list->begin(); ti != res_list->end(); ++ti) {
          (ti->second).print();
          cout << "Cost " << ti->first << "\n\n" << flush;
          delete_circuit(ti->second);
        }
        res_list->clear();
        if (config::early_stop) {
          pthread_mutex_destroy(&result_mutex);
          delete [] circ_table;
          if (left_table != NULL)
            delete [] left_table;
          delete base_list;
          return;
        }
      }
      cout << "----------------------------------------\n" << flush;
      pthread_mutex_destroy(&result_mutex);
    }
  }
  delete [] circ_table;
  if (left_table != NULL)
    delete [] left_table;
  delete base_list;
}

void * worker_thrd_tdepth(void * arg) {
  Circuit circ, tmp_circ, tmp_circ2, tmp_circ3, ans;
  int k, i, j, cst;
  bool left_multiply;
  map_t * mp;
  Rmatrix V(dim, dim), W(dim_proj, dim), U = *((Rmatrix *)arg);
  Canon * canon_form;
  struct triple * trip;

  pthread_mutex_lock(&data_lock);
  while(1) {
    if (data_avail) {
      // Copy our data
      circ = data_circ;
      k = data_k;
      left_multiply = data_left;
      if (k % 2 == 0) {
        data_mat.permute_adj(V, k/2);
      } else {
        data_mat.permute(V, k/2);
      }
      mp = data_map;
      data_avail = false;
      data_num++;
      // Signal the mamma thread
      pthread_cond_signal(&thrd_ready);
      // Unlock the lock
      pthread_mutex_unlock(&data_lock);

      // Perform the search
      if (dim != dim_proj) {
        V.submatrix(0, 0, dim_proj, dim, W);
        canon_form = left_multiply ? canonicalize(W*U, false, false, false) : 
                                     canonicalize(U*W, false, false, false);
      } else {
        canon_form = left_multiply ? canonicalize(V*U) : canonicalize(U*V);
      }

      trip = &(canon_form->front());
      ans = find_unitary(trip->key, trip->mat, *mp).second;
      if (!ans.empty()) {
        // Generate the two circuit halves
        tmp_circ = (k == 0) ? circ : circ.transform(k/2, k % 2 == 1);
        tmp_circ2 = ans.transform(-(trip->permutation), trip->adjoint);

        // Check that the circuit is correct
        if (left_multiply ? check_it(tmp_circ, tmp_circ2, U) : 
                            check_it(tmp_circ2, tmp_circ, U)) {
          tmp_circ3 = left_multiply ? tmp_circ.append(tmp_circ2) : 
                                      tmp_circ2.append(tmp_circ);
          cst = tmp_circ3.cost();
          pthread_mutex_lock(&prnt_lock);
          data_res->insert(ord_circuit_pair(cst, tmp_circ3));
          pthread_mutex_unlock(&prnt_lock);
        }

        if (k != 0) delete_circuit(tmp_circ);
        delete_circuit(tmp_circ2);
      }
      canon_form->clear();
      delete canon_form;

      // Lock before we reenter the loop
      pthread_mutex_lock(&data_lock);
      data_num--;
    } else {
      // Wait until data is ready
      pthread_cond_wait(&data_ready, &data_lock);
    }
  }
  pthread_exit(NULL);
}

void exact_search_tdepth(Rmatrix & U) {
  int i, j, k, ind;
  int pe = config::mod_perms ? num_perms : 1;
  int in = config::mod_invs  ?         1 : 2;
  map_t * circ_table = new map_t[config::max_seq];
  map_t * left_table = (config::ancilla == 0) ? NULL : new map_t[config::max_seq];
  map_t * mp         = (config::ancilla == 0) ? circ_table : left_table;
  map_iter it;
  circuit_list * base_list;
  ord_circuit_list * res_list = data_res = new ord_circuit_list;
  ord_circuit_iter ti;
  struct timespec start, end;

  circuit_iter c;
  Rmatrix tmp(dim, dim);

  // Do this first so that the threads don't conflict
  base_list = generate_base_circuits();

  //Initialize the circuit tables
  //Threading stuff
  pthread_mutex_init(&data_lock, NULL);
  pthread_mutex_init(&prnt_lock, NULL);
  pthread_cond_init(&data_ready, NULL);
  pthread_cond_init(&thrd_ready, NULL);

  pthread_t * thrds = new pthread_t[config::num_threads];
  for (i = 0; i < config::num_threads; i++) {
    pthread_create(thrds + i, NULL, &(worker_thrd_tdepth), (void *)(&U));
  }
  //---------------------

  load_sequences(0, base_list, circ_table, NULL);
  pthread_mutex_lock(&data_lock);
  for (i = 1; i < config::max_seq; i++) {
    if (i % 2 == 1) {
      load_sequences(i, cliff_list, circ_table, NULL);
    } else {
      data_left = true;
      load_sequences(i, base_list, circ_table, NULL);
    }

    cout << "Looking for circuits...\n";
    clock_gettime(CLOCK_MONOTONIC, &start);

    if (i <= 1) {
      // Search the cliffords
      data_map = mp + i;
      data_left = false;
      data_circ = Circuit(1);
      data_circ.to_Rmatrix(data_mat);
      data_k = 0;
      data_avail = true;
      pthread_cond_signal(&data_ready);
      pthread_cond_wait(&thrd_ready, &data_lock);
      delete_circuit(data_circ);
    } else {
      for (i == 2 ? j = 1 : j = 0; j < 2; j++) {
        if (i % 2 == 1) {
          // mp[i-1] = TCTC..., search for mp[i-1]mp[i-1] and mp[i]mp[i-1] 
          //  to get t-depth i-1
          data_left = false;
          data_map = mp + i + j - 1;
          ind = i - 1;
        } else {
          // mp[i-1] = CTCTC..., search for mp[i-2]mp[i] and mp[i-1]mp[i] 
          //  to get t-depth i-1
          data_left = true;
          data_map = mp + i;
          ind = i + j - 2;
        }
        for (it = mp[ind].begin(); it != mp[ind].end(); it++) {
          data_circ = it->second;
          (it->second).to_Rmatrix(data_mat);
          for (k = 0; k < 2*pe; k += in) {
            // Write data
            data_k = k;
            data_avail = true;
            // Signal workers that data is ready
            pthread_cond_signal(&data_ready);
            pthread_cond_wait(&thrd_ready, &data_lock);
          }
        }
      }
    }

    /* Wait for all threads to finish */
    while (data_num > 0) {
      pthread_mutex_unlock(&data_lock);
      pthread_mutex_lock(&data_lock);
    }
    for (ti = res_list->begin(); ti != res_list->end(); ++ti) {
      (ti->second).print();
      cout << "Cost " << ti->first << "\n\n" << flush;
      delete_circuit(ti->second);
      res_list->erase(ti);
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    cout << fixed << setprecision(3);
    cout << "Time: " << (end.tv_sec + (double)end.tv_nsec/1000000000) - (start.tv_sec + (double)start.tv_nsec/1000000000) << " s\n";
    cout << "----------------------------------------\n" << flush;
  }

  delete [] circ_table;
  delete [] left_table;
  delete [] base_list;
  delete [] thrds;
}

//------------------------------------------------Approximations

double circuit_dist(const Circuit & a, const Circuit & b) {
	Rmatrix U(dim, dim), V(dim, dim);
	a.to_Rmatrix(U);
	b.to_Rmatrix(V);
	return dist(U, V);
}

// Create a projector onto [0,1] that closes over the circuit's unitary
template<typename T>
class circuit_closure {
	private:
		T p;
	public:
		circuit_closure(const Circuit & circ) : p(dim, dim) { circ.to_matrix(p); }
		circuit_closure(const T & mat) : p(mat) { }
		double operator()(const Circuit & circ) const { 
			T q(dim, dim);
			circ.to_matrix(q);
			return dist(p, q);
		}
};


template<typename T, class NNtree>
void * worker_thrd_approx(void * arg) {
  Circuit circ, tmp_circ, tmp_circ2, tmp_circ3;
	const Circuit * ans;
  int k, l;
  bool left_multiply;
  double epsilon = config::precision;
  NNtree * nntable;
  T V(dim, dim), U = *((T *)arg);

  pthread_mutex_lock(&data_lock);
  while(1) {
    if (data_avail) {
      // Copy our data
      circ = data_circ;
      k = data_k;
      l = data_l;
      left_multiply = data_left;

		  if (l % 2 == 1) {
			  permute_adj(*((T *)data_mat_approx), V, l/2);
		  } else {
			  permute(*((T *)data_mat_approx), V, l/2);
		  }

      nntable = (NNtree *)data_nntree;
      data_avail = false;
      data_num++;
      // Signal the mamma thread
      pthread_cond_signal(&thrd_ready);
      // Unlock the lock
      pthread_mutex_unlock(&data_lock);

      // Perform the search
      pthread_mutex_lock(&prnt_lock);
      epsilon = config::precision;
      pthread_mutex_unlock(&prnt_lock);

      ans = nntable->nearest_neighbour(circuit_closure<T>(V), &epsilon);
      if (ans != NULL) {
        // Generate the two circuit halves
        if (left_multiply) {
          tmp_circ = circ.transform(0, true);
          tmp_circ2 = tmp_circ.append(ans->transform(-(l/2), l % 2 == 1));
        } else {
          tmp_circ = ans->transform(-(l/2), l % 2 == 1);
          tmp_circ2 = tmp_circ.append(circ.transform(0, true));
        }
        tmp_circ3 = (k == 0) ? tmp_circ2 : tmp_circ2.transform(-(k/2), k % 2 == 1);

        pthread_mutex_lock(&prnt_lock);
        data_res_approx->push_back(make_pair(epsilon, tmp_circ3));
        if (epsilon < config::precision) config::precision = epsilon;
        pthread_mutex_unlock(&prnt_lock);

        if (k != 0) delete_circuit(tmp_circ2);
        delete_circuit(tmp_circ);
      }

      // Lock before we reenter the loop
      pthread_mutex_lock(&data_lock);
      data_num--;
    } else {
      // Wait until data is ready
      pthread_cond_wait(&data_ready, &data_lock);
    }
  }
  pthread_exit(NULL);
}
/*
template<typename T>
void approx_search_gen(T & U) {
  typedef VPTree<Circuit, circuit_dist, circuit_closure<T> > NNtree;
  typedef typename NNtree::iterator             NN_iter;
  typedef typename NNtree::const_iterator const_NN_iter;

  int num = 0, p = 0;
  int i, j, k, l;
  int pe = config::mod_perms ? num_perms : 1;
  int in = config::mod_invs  ?         1 : 2;
  map_t  * circ_table = new map_t [config::max_seq];
	NNtree * NN_table   = new NNtree[config::max_seq];
  NN_iter it;
  circuit_list * base_list;
  struct timespec start, end;
	T mat(dim, dim), V(dim, dim), W(dim, dim);
  T equivs[2*pe];
	const Circuit * ans;
	Circuit tmp_circ, tmp_circ2, tmp_circ3;
	double epsilon = config::precision;

  list<pair<double, Circuit> > res_list;
  list<pair<double, Circuit> >::iterator ti;

  // Store all permutations and inversions
  for (k = 0; k < 2*pe; k += in) {
		if (k % 2 == 1) {
			permute_adj(U, equivs[k], k/2);
		} else {
			permute(U, equivs[k], k/2);
		}
  }

  // Do this first so that the threads don't conflict
  base_list = generate_base_circuits();

  load_sequences(0, base_list, circ_table, NULL);
  for (i = 1; i < config::max_seq; i++) {
   	load_sequences(i, base_list, circ_table, NULL);
		NN_table[i].build_tree(
				map_value_iter(circ_table[i].begin()),
			 	map_value_iter(circ_table[i].end()),
			 	circ_table[i].size()
		);

    // Meet in the middle - Sequences of length 2i + {0, 1}
    for (j = max(i-1, 1); j <= i; j++) {
    	cout << "Looking for circuits with depth " << 2*i - (i - j) << "...\n";
			cout << "|";
			num = 0;
			p = 0;
    	clock_gettime(CLOCK_MONOTONIC, &start);

			// Look for circuits
      for (it = NN_table[j].begin(); it != NN_table[j].end(); it++) {
        it->to_matrix(mat);
        for (k = 0; k < 2*pe; k += in) {
          V = equivs[k]*mat;
          for (l = 0; l < 2*pe; l += in) {

					  if (l % 2 == 1) {
						  permute_adj(V, W, l/2);
					  } else {
						  permute(V, W, l/2);
					  }

            ans = NN_table[i].nearest_neighbour(circuit_closure<T>(W), &epsilon);
            if (ans != NULL) {
              // Generate the two circuit halves
              tmp_circ = ans->transform(-(l/2), l % 2 == 1);
              tmp_circ2 = tmp_circ.append(it->transform(0, true));
              tmp_circ3 = (k == 0) ? tmp_circ2 : tmp_circ2.transform(-(k/2), k % 2 == 1);

              res_list.push_back(make_pair(epsilon, tmp_circ3));

              if (k != 0) delete_circuit(tmp_circ2);
              delete_circuit(tmp_circ);
            }
          }
        }
				num++;
        int xxx = num*37 / (NN_table[j].size());
				if (xxx >= p) {
          for (int yyy = p; yyy <= ceil(xxx); yyy++) {
					  cout << "=" << flush;
					  p++;
          }
				}
      }
			cout << "|\n";

			for (ti = res_list.begin(); ti != res_list.end(); ++ti) {
				(ti->second).print();
				cout << "Distance " << scientific << ti->first << "\n\n" << flush;
				delete_circuit(ti->second);
			}
      res_list.clear();

			clock_gettime(CLOCK_MONOTONIC, &end);
			cout << fixed << setprecision(3);
			cout << "Time: " << (end.tv_sec + (double)end.tv_nsec/1000000000) - (start.tv_sec + (double)start.tv_nsec/1000000000) << " s\n";
			cout << "----------------------------------------\n" << flush;

    }
  }

  delete [] circ_table;
  delete base_list;
}
*/
template<typename T>
void approx_search_gen(T & U) {
  typedef VPTree<Circuit, circuit_dist, circuit_closure<T> > NNtree;
  typedef typename NNtree::iterator             NN_iter;
  typedef typename NNtree::const_iterator const_NN_iter;

  int num = 0, p = 0;
  int i, j, k, l, ind, len;
  int pe = config::mod_perms ? num_perms : 1;
  int in = config::mod_invs  ?         1 : 2;
  map_t  * circ_table = new map_t [config::max_seq];
	NNtree * NN_table   = new NNtree[config::max_seq];
  NN_iter it;
  circuit_list * base_list;
  struct timespec start, end;
	T mat(dim, dim), V(dim, dim), W(dim, dim);
  T equivs[2*pe];
	const Circuit * ans;
	Circuit tmp_circ, tmp_circ2, tmp_circ3;
	double epsilon = config::precision;

  list<pair<double, Circuit> > res_list;
  list<pair<double, Circuit> >::iterator ti;

  data_mat_approx = &V;
  data_res_approx = &res_list;

  // Store all permutations and inversions
  for (k = 0; k < 2*pe; k += in) {
		if (k % 2 == 1) {
			permute_adj(U, equivs[k], k/2);
		} else {
			permute(U, equivs[k], k/2);
		}
  }

  //Threading stuff
  pthread_mutex_init(&data_lock, NULL);
  pthread_mutex_init(&prnt_lock, NULL);
  pthread_cond_init(&data_ready, NULL);
  pthread_cond_init(&thrd_ready, NULL);

  pthread_t * thrds = new pthread_t[config::num_threads];
  for (i = 0; i < config::num_threads; i++) {
    pthread_create(thrds + i, NULL, &(worker_thrd_approx<T, NNtree>), (void *)(&U));
  }

  // Do this first so that the threads don't conflict
  base_list = generate_base_circuits();


  load_sequences(0, base_list, circ_table, NULL);
  pthread_mutex_lock(&data_lock);
  for (i = 1; i < config::max_seq; i++) {
    if (config::tdepth && i % 2 == 1) {
      load_sequences(i, cliff_list, circ_table, NULL);
    } else {
   	  load_sequences(i, base_list, circ_table, NULL);
    }
    cout << "----------------------------------------\n" << flush;
    cout << "Building VP tree...\n";
    clock_gettime(CLOCK_MONOTONIC, &start);
		NN_table[i].build_tree(
				map_value_iter(circ_table[i].begin()),
			 	map_value_iter(circ_table[i].end()),
			 	circ_table[i].size()
		);
    clock_gettime(CLOCK_MONOTONIC, &end);
    cout << fixed << setprecision(3);
    cout << "Time: " << (end.tv_sec + (double)end.tv_nsec/1000000000) - (start.tv_sec + (double)start.tv_nsec/1000000000) << " s\n";
    cout << "----------------------------------------\n" << flush;

    if (!config::tdepth) {
      // Meet in the middle - Sequences of length 2i + {0, 1}
      data_nntree = NN_table + i;
      for (j = max(i-1, 1); j <= i; j++) {
        cout << "Looking for circuits with depth " << 2*i - (i - j) << "...\n";
        cout << "|";
        num = 0;
        p = 0;
        clock_gettime(CLOCK_MONOTONIC, &start);

        // Look for circuits
        for (it = NN_table[j].begin(); it != NN_table[j].end(); it++) {
          it->to_matrix(mat);
          data_circ = *it;
          for (k = 0; k < 2*pe; k += in) {
            V = equivs[k]*mat;
            for (l = 0; l < 2*pe; l += in) {
              // Write data
              data_l = l;
              data_k = k;
              data_avail = true;
              // Signal workers that data is ready
              pthread_cond_signal(&data_ready);
              pthread_cond_wait(&thrd_ready, &data_lock);
            }
          }
          num++;
          int xxx = num*37 / (NN_table[j].size());
          if (xxx >= p) {
            for (int yyy = p; yyy <= ceil(xxx); yyy++) {
              cout << "=" << flush;
              p++;
            }
          }
        }
        cout << "|\n";

        // Wait for all threads to finish 
        while (data_num > 0) {
          pthread_mutex_unlock(&data_lock);
          pthread_mutex_lock(&data_lock);
        }
        for (ti = res_list.begin(); ti != res_list.end(); ++ti) {
          (ti->second).print();
          cout << "Distance " << scientific << ti->first << "\n\n" << flush;
          delete_circuit(ti->second);
        }
        res_list.clear();

        clock_gettime(CLOCK_MONOTONIC, &end);
        cout << fixed << setprecision(3);
        cout << "Time: " << (end.tv_sec + (double)end.tv_nsec/1000000000) - (start.tv_sec + (double)start.tv_nsec/1000000000) << " s\n";
        cout << "----------------------------------------\n" << flush;
      }
    } else {
      // Meet in the middle - Sequences of T-depth i-1
      cout << "Looking for circuits with T-depth " << i - 1 << "...\n";
      cout << "|";
      num = 0;
      p = 0;
      clock_gettime(CLOCK_MONOTONIC, &start);
      if (i <= 1) {
        // Search the cliffords
        data_nntree = NN_table + i;
        data_circ = Circuit(1);
        for (k = 0; k < 2*pe; k += in) {
          V = equivs[k]*mat;
          data_l = 0;
          data_k = k;
          data_avail = true;
          // Signal workers that data is ready
          pthread_cond_signal(&data_ready);
          pthread_cond_wait(&thrd_ready, &data_lock);
        }
        delete_circuit(data_circ);
      } else {
        for (i == 2 ? j = 1 : j = 0; j < 2; j++) {
          if (i % 2 == 1) {
            // mp[i-1] = TCTC..., search for mp[i-1]mp[i-1] and mp[i]mp[i-1] 
            //  to get t-depth i-1
            data_left = false;
            data_nntree = NN_table + i + j - 1;
            ind = i - 1;
          } else {
            // mp[i-1] = CTCTC..., search for mp[i-2]mp[i] and mp[i-1]mp[i] 
            //  to get t-depth i-1
            data_left = true;
            data_nntree = NN_table + i;
            ind = i + j - 2;
          }
          for (it = NN_table[ind].begin(); it != NN_table[ind].end(); it++) {
            it->to_matrix(mat);
            data_circ = *it;
            for (k = 0; k < 2*pe; k += in) {
              V = data_left ? mat*equivs[k] : equivs[k]*mat;
              for (l = 0; l < 2*pe; l += in) {
                // Write data
                data_l = l;
                data_k = k;
                data_avail = true;
                // Signal workers that data is ready
                pthread_cond_signal(&data_ready);
                pthread_cond_wait(&thrd_ready, &data_lock);
              }
            }
            num++;
            int xxx = num*(18 + j) / (NN_table[ind].size());
            if (xxx >= p) {
              for (int yyy = p; yyy <= ceil(xxx); yyy++) {
                cout << "=" << flush;
                p++;
              }
            }
          }
        }
        cout << "|\n";

        // Wait for all threads to finish 
        while (data_num > 0) {
          pthread_mutex_unlock(&data_lock);
          pthread_mutex_lock(&data_lock);
        }
        for (ti = res_list.begin(); ti != res_list.end(); ++ti) {
          (ti->second).print();
          cout << "Distance " << scientific << ti->first << "\n\n" << flush;
          delete_circuit(ti->second);
        }
        res_list.clear();

        clock_gettime(CLOCK_MONOTONIC, &end);
        cout << fixed << setprecision(3);
        cout << "Time: " << (end.tv_sec + (double)end.tv_nsec/1000000000) - (start.tv_sec + (double)start.tv_nsec/1000000000) << " s\n";
        cout << "----------------------------------------\n" << flush;
      }
    }
  }

  delete [] circ_table;
  delete base_list;
}

void approx_search(Rmatrix & U) { approx_search_gen<Rmatrix>(U); }
void approx_search(Unitary & U) { approx_search_gen<Unitary>(U); }

void mem_test(int n) {
  map_t * circ_table = new map_t[config::max_seq];
  circuit_list * base_list;
  base_list = generate_base_circuits();
  for (int i = 0; i <= n; i++) {
    load_sequences(i, base_list, circ_table, NULL);
  }
  delete [] circ_table;
  delete base_list;
}
