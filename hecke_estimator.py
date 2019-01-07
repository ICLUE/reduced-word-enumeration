from random import randint, shuffle
import random
import time
import multiprocessing
import math

TRIAL_BLOCK_SIZE = 100000

def main():
    permutation = input('Permutation (comma separated): ')
    permutation = list(map(lambda x: int(x.strip()), permutation.split(',')))
    e = int(input('e: '))
    num_trials = int(input('Number of trials T: '))

    start = time.time()
    est_hecke_words = str(estimate_num_hecke_words(permutation, e, num_trials))
    end = time.time()
    print('Executed in {} seconds'.format(int(end - start)))

    print('\nEstimated number of Hecke words:\n{}.{} x 10^{}'.format(
        est_hecke_words[0],
        est_hecke_words[1:8],
        len(est_hecke_words) - 1)
    )

success_count = 0

def estimate_num_hecke_words(permutation, e, num_trials, seed=None):
    global success_count

    if seed is not None:
        random.seed(seed)

    # Run T trials, sum success values (S)
    # Run trials in large blocks, take advantage of multiple threads
    num_workers = multiprocessing.cpu_count();
    trials_left = num_trials
    while trials_left > TRIAL_BLOCK_SIZE:
        pool = multiprocessing.Pool(num_workers)
        for i in range(TRIAL_BLOCK_SIZE):
            a = pool.apply_async(run_trial, args=(permutation[:], e), callback=process_result)
        pool.close()
        pool.join()
        trials_left -= TRIAL_BLOCK_SIZE
    pool = multiprocessing.Pool(num_workers)
    for i in range(trials_left):
        a = pool.apply_async(run_trial, args=(permutation[:], e), callback=process_result)
    pool.close()
    pool.join()

    # S / T = #red(w)
    est_hecke_words = (success_count + num_trials // 2) // num_trials
    return int(est_hecke_words)

def process_result(val):
    global success_count
    success_count += val

def run_trial(current_perm, e):
    perm_size = len(current_perm)

    # Determine initial permutation length
    perm_len = 0
    for i in range(perm_size - 1):
        for j in range(i + 1, perm_size):
            if current_perm[i] > current_perm[j]:
                perm_len += 1

    current_prob = 1

    j_range = reversed(range(1, e + 1))
    for j in j_range:
        descent_list = []
        for i in range(perm_size - 1):
            if current_perm[i] > current_perm[i + 1]:
                descent_list.append(i)
        if len(descent_list) == 0:
            if e > 0:
                current_prob = 0
            if e == 0:
                current_prob = 1
        else:
            exchange_idx = descent_list[randint(0, len(descent_list) - 1)]
            if j == perm_len:
                # Perform exchange and reduce permutation length
                current_perm[exchange_idx], current_perm[exchange_idx + 1] = current_perm[exchange_idx + 1], current_perm[exchange_idx]
                current_prob *= len(descent_list)
                perm_len -= 1
            else:
                if perm_len == 1 and e > 1:
                    current_prob *= len(descent_list)
                else:
                    if randint(0, 1):
                        current_prob *= 2 * len(descent_list)
                    else:
                        # Perform exchange and reduce permutation length
                        current_perm[exchange_idx], current_perm[exchange_idx + 1] = current_perm[exchange_idx + 1], current_perm[exchange_idx]
                        current_prob *= 2 * len(descent_list)
                        perm_len -= 1

    return current_prob

if __name__ == '__main__':
    main()
