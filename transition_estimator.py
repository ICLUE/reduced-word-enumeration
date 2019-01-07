from random import randint, shuffle
import random
import time
import multiprocessing
import math

TRIAL_BLOCK_SIZE = 100000

def main():
    permutation = input('Permutation (comma separated): ')
    permutation = list(map(lambda x: int(x.strip()), permutation.split(',')))
    num_trials = int(input('Number of trials T: '))

    start = time.time()
    est_reduced_words = str(estimate_num_reduced_words(permutation, num_trials))
    end = time.time()
    print('Executed in {} seconds'.format(int(end - start)))

    print('\nEstimated number of reduced words:\n{}.{} x 10^{}'.format(
        est_reduced_words[0],
        est_reduced_words[1:8],
        len(est_reduced_words) - 1)
    )

success_count = 0

def process_result(val):
    global success_count
    success_count += val

def estimate_num_reduced_words(permutation, num_trials, seed=None):
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
            a = pool.apply_async(transition_importance, args=([permutation[:]]), callback=process_result)
        pool.close()
        pool.join()
        trials_left -= TRIAL_BLOCK_SIZE
    pool = multiprocessing.Pool(num_workers)
    for i in range(trials_left):
        a = pool.apply_async(transition_importance, args=([permutation[:]]), callback=process_result)
    pool.close()
    pool.join()

    # S / T = #red(w)
    est_reduced_words = success_count // num_trials
    return int(est_reduced_words)


def makeSwap(w, i, j):
    ans = w[:]
    ans[i - 1] = w[j - 1]
    ans[j - 1] = w[i - 1]
    return ans

def getDescents(w):
    ans  = []
    for i in range(len(w) - 1):
        if w[i] > w[i + 1]: ans.append(i)
    return ans

def getInv(w):
    ans = [w.index(i) + 1 for i in range(1, len(w) + 1)]
    return ans

def getCode(w):
    ans = [0] * len(w)
    for i in range(len(w)):
        for j in range(i + 1, len(w)):
            if w[j] < w[i]:
                ans[i] += 1
    return ans

# Returns conjugate permutation given Lehmer code
def conjugate(shape):
    ans = [0] * shape[0]
    for i in range(len(shape)):
        for j in range(shape[i]):
            ans[j] += 1
    return ans

# Calculate hook length of permutation specified by Lehmer code
def hookLength(shape):
    ans = math.factorial(sum(shape))
    conj = conjugate(shape)

    for i in range(len(shape)):
        for j in range(shape[i]):
            ans /= (shape[i] - i + conj[j] - j - 1)
    return ans

def findLastDescent(w):
    n = len(w)
    descents = getDescents(w)

    last_descent = -1
    for i in range(len(descents) - 1, -1, -1):
        d = descents[i] + 1

        w2 = []
        for j in range(1, d):
            if w[j - 1] < w[d - 1]:
                w2.append(w[j - 1])
        if len(w2) == 0: continue
        w2.sort()
        check = False
        for i in range(1, len(w2)):
            if w2[i] > w2[i - 1] + 1:
                check = True
                break
        if w[d - 1] > w2[ - 1] + 1: check = True

        if len(w2) > 0 and check:
            last_descent = d
            break

    if last_descent == -1:
        return -1, -1
    boxes_in_ld = []

    winv = getInv(w)

    maxCol = 0
    for j in range(1, w[last_descent - 1]):
        if winv[j - 1] > last_descent:
            boxes_in_ld.append((last_descent, j))
            if j > maxCol: maxCol = j

    bps = maxCol
    while (last_descent, bps) in boxes_in_ld:
        bps -= 1


    to_switch = last_descent + 1

    for i in range(last_descent + 1, n + 1):
        if w[i - 1] > bps and w[i - 1] < w[last_descent - 1]: to_switch = i

    return last_descent, to_switch


def transition(w):
    ans = []

    last_descent, to_switch = findLastDescent(w)

    if last_descent == -1:
        return []
    ans.append(last_descent)

    wprime = makeSwap(w, last_descent, to_switch)
    ans.append(wprime)

    w_dprime_list = []
    rows_list = []

    for a in range(1, last_descent):
        if wprime[a - 1] < wprime[last_descent - 1]:
            check = True
            for c in range(a + 1, last_descent):
                if wprime[a - 1] < wprime[c - 1] and wprime[c - 1] < wprime[last_descent - 1]:
                    check = False
                    break
            if check:
                w_dprime_list.append(makeSwap(wprime, a, last_descent))
                rows_list.append(a)

    ans.append(w_dprime_list)
    ans.append(rows_list)
    return ans

def transition_importance(w):
    ans = 1
    cur = w
    while True:
        # Generate all children accessible via one transition step
        # and choose one
        children = transition(cur)
        # Vexillary case
        if len(children) == 0:
            break
        ans *= len(children[2])
        cur = random.choice(children[2])

    # Multiply by hook length of vexillary permutation
    c = getCode(cur)
    c.sort(reverse=True)
    ans *= hookLength(c)

    return ans


if __name__ == '__main__':
    main()
